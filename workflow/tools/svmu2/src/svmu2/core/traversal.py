'''
Docstring for core.traversal
Produce a traversal through the graph representation of the alignment
'''    
import networkx as nx
from math import sqrt
from sortedcontainers import SortedList
from dataclasses import dataclass

@dataclass
class TrendResult:
    trend: str
    p_value: float
    significant: bool
    monotonic: bool

@dataclass(frozen=True)
class VirtualNode:
    reference_start: float
    query_start: float
    reference_end: float
    query_end: float

def is_coordinate_tuple(obj):
    ''' Check for malformed tuples when coords are used for source and sink instead of true nodes '''
    return (
        isinstance(obj, (tuple, list)) and
        len(obj) == 2 and
        all(isinstance(x, (int, float)) for x in obj)
    )

def find_source_sink_nodes(alignment, trend_result):
    """
    Determine source and sink blocks for graph traversal.
    These should be the most extreme blocks in the traversal of a monotonically increasing / decreasing plot

    returns:
        (source, sink)

    Poised for deprecation under new paradigm
    """
    # Choose corners based on trend
    if trend_result.trend == "increasing":
        start = alignment.bounding_box.bottom_left
        end   = alignment.bounding_box.top_right
    elif trend_result.trend == "decreasing":
        start = alignment.bounding_box.top_left
        end   = alignment.bounding_box.bottom_right
    else:
        raise RuntimeError(
            "Cannot determine source/sink for non-monotonic alignment"
        )

    current_min_start = None
    current_min_end = None
    min_dist_start = float("inf")
    min_dist_end = float("inf")

    for block in alignment.alignment_blocks:
        d_start = alignment.euclidean_distance(start, (block.reference_start, block.query_start))
        d_end = alignment.euclidean_distance(end, (block.reference_end, block.query_end))

        if d_start < min_dist_start:
            min_dist_start = d_start
            current_min_start = block

        if d_end < min_dist_end:
            min_dist_end = d_end
            current_min_end = block

    if current_min_start is None or current_min_end is None:
        raise RuntimeError("Failed to identify source/sink blocks (do you have an alignment with no alignment blocks?)")

    if current_min_start is current_min_end:
        raise RuntimeError("Degenerate traversal: source == sink. Not necessarily fatal but erroring out")

    source = current_min_start
    sink = current_min_end

    return (source, sink)

def evaluate_alignment_trend(aln):
    
    import numpy as np
    import pymannkendall as mk
    aln.trace = [b.weighted_reference_length for b in aln.alignment_blocks]
    aln.cumsum = np.cumsum(aln.trace)

    if len(aln.cumsum) < 3:
        trend_result = TrendResult(
            trend="undetermined",
            p_value=float("nan"),
            significant=False,
            monotonic=False,
        )

        return trend_result

    result = mk.hamed_rao_modification_test(aln.cumsum, alpha=0.2)

    monotonic = result.trend in ("increasing", "decreasing")
    trend_result = TrendResult(
        trend=result.trend,
        p_value=result.p,
        significant=result.h,
        monotonic=monotonic
    )

    return trend_result

# --- 1. Distance & Spatial Helpers ---

def _xdistance(a, b):
    return abs(b.reference_start - a.reference_end)

def _ydistance(a, b):
    return abs(b.query_start - a.query_end)

def _euclidean_distance(x1, y1, x2, y2):
    return sqrt((x2 - x1)**2 + (y2 - y1)**2)

def _orientation_aware_euclidean_distance(a, b):
    """
    Find the min distance between two alignment block objects, allowing for 
    reflections at any point defined on the plane.
    """
    scores = [
        _euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.query_start),
        _euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.query_end),
        _euclidean_distance(a.reference_end, a.query_start, b.reference_start, b.query_start),
        _euclidean_distance(a.reference_end, a.query_end, b.reference_end, b.query_start)
    ]

    if getattr(a, "reflected", False):
        scores.append(_euclidean_distance(a.reference_end, a.y2_reflection, b.reference_start, b.query_start))
    
    if getattr(b, "reflected", False):
        scores.append(_euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.y1_reflection))
        
    return min(scores)

def _get_reference_neighbors(ref_index, end_coord, window):
    left = ref_index.bisect_left((end_coord - window, -1))
    right = ref_index.bisect_right((end_coord + window, float('inf')))

    return ref_index[left:right]

# --- 2. Graph Construction ---

def build_alignment_graph(blocks, source, sink, k_nearest = 50, max_jump=100_000, weights='euclidean'):
    """
    Builds a directed graph of alignment blocks with edges weighted by spatial distance.
    
    Returns:
        tuple: (nx.DiGraph, resolved_source, resolved_sink)

    Note: I've found that using _xdistance to/from source / sink node works best around edges. 
    """
    reference_index = SortedList([(block.reference_start, i, block) for i, block in enumerate(blocks)])

    G = nx.DiGraph()
    for i, block in enumerate(blocks):
        G.add_node(i, block=block)

    weight_fn = {
        'euclidean': _orientation_aware_euclidean_distance,
        'X': _xdistance,
        'Y': _ydistance
    }.get(weights)

    if weight_fn is None:
        raise ValueError(f"Invalid weight designation for Dijkstra traversal: {weights}")

    # Build edges between existing blocks
    for i, a in enumerate(blocks):
        neighbors = _get_reference_neighbors(reference_index, a.reference_end, max_jump)
        for _, j, b in neighbors:
            if i == j:
                continue
            distance = weight_fn(a, b)
            G.add_edge(i, j, weight=distance)

    # Resolve Source
    if is_coordinate_tuple(source):
        source_node = VirtualNode(
            reference_start=source[0],
            reference_end=source[0],
            query_start=source[1],
            query_end=source[1]
        )
        G.add_node(source_node)
        
        # Grab the K nearest blocks downstream (reference_start >= source reference_end)
        start_idx = reference_index.bisect_left((source_node.reference_end, -1))
        neighbors = reference_index[start_idx : start_idx + k_nearest] # source not in ref idx, dont ever add it!
        
        for _, j, b in neighbors:
            distance = weight_fn(source_node, b) ## this is extremely sensitive
            G.add_edge(source_node, j, weight=distance)
        resolved_source = source_node
        
    elif isinstance(source, int):
        resolved_source = source
    else:
        raise TypeError("Source must be an integer block index or a (x, y) coordinate tuple.")

    # Resolve Sink
    if is_coordinate_tuple(sink):
        sink_node = VirtualNode(
            reference_start=sink[0],
            reference_end=sink[0],
            query_start=sink[1],
            query_end=sink[1]
        )
        G.add_node(sink_node)
        
        # Grab the K nearest blocks upstream (reference_start <= sink reference_start)
        end_idx = reference_index.bisect_left((sink_node.reference_start, -1))
        start_idx = max(0, end_idx - k_nearest)
        neighbors = reference_index[start_idx : end_idx] # sink not in ref idx, dont ever add it!
        
        for _, j, b in neighbors:
            distance = weight_fn(b, sink_node) ## this is extremely sensitive
            G.add_edge(j, sink_node, weight=distance)
        resolved_sink = sink_node
        
    elif isinstance(sink, int):
        resolved_sink = sink
    else:
        raise TypeError("Sink must be an integer block index or a (x, y) coordinate tuple.")

    return G, resolved_source, resolved_sink

# --- 3. Graph Traversal ---

def execute_dijkstra_traversal(G, source, sink, return_graph=False):
    ## traverse pre built graph
    try:
        # Use the variable arguments, not string literals
        path = nx.dijkstra_path(G, source, sink, weight='weight')
        total_distance = sum(G[path[i]][path[i + 1]]['weight'] for i in range(len(path) - 1))

    except nx.NetworkXNoPath:
        path = []
        total_distance = float('inf')

    # Return path_blocks without leaking virtual nodes
    path_blocks = [G.nodes[i]['block'] for i in path if 'block' in G.nodes[i] and isinstance(i, int)]

    result = {
        'source': source,
        'sink': sink,
        'path': path,
        'total_distance': total_distance,
        'path_blocks': path_blocks
    } 

    if return_graph:
        result['graph'] = G

    return result

# --- 4. Orchestration Wrapper (Backward Compatibility) ---

def dijkstra_traversal(blocks, source, sink, max_jump=100_000, weights='euclidean', return_graph=False):
    """
    Perform a Dijkstra traversal on alignment blocks with optional orientation-aware weights.
    """
    G, resolved_source, resolved_sink = build_alignment_graph(
        blocks=blocks, 
        source=source, 
        sink=sink, 
        max_jump=max_jump, 
        weights=weights
    )

    return execute_dijkstra_traversal(
        G=G, 
        source=resolved_source, 
        sink=resolved_sink, 
        return_graph=return_graph
    )

def compute_primary_synteny(alignment, max_jump=100_000, weights="euclidean"):
    if getattr(alignment, "source", None) is None or getattr(alignment, "sink", None) is None:
        raise RuntimeError(
            f"compute_primary_synteny called on non-traversable alignment "
            f"{alignment.reference} vs {alignment.query}"
        )

    blocks = alignment.alignment_blocks

    try:
        source_idx = blocks.index(alignment.source)
        sink_idx   = blocks.index(alignment.sink)
    except ValueError:
        raise RuntimeError("Source or sink block not found in alignment blocks list")

    traversal = dijkstra_traversal(
        blocks,
        source_idx,
        sink_idx,
        max_jump=max_jump,
        weights=weights
    )

    alignment.set_primary_synteny_blocks(traversal["path_blocks"])

def reflect_alignment_blocks(alignments, trend):
    """
    Reflects alignment blocks across their assigned pivot (x-axis reflection).
    Assumes that each block needing reflection has a `.pivot` attribute set.
    Updates the reference coordinates in-place.
    """
    for aln in alignments:
        for block in aln.alignment_blocks:
            if aln.slope != trend:
                point = find_xeqc_point_intersection(block)
                print(point)
                reflect_line_segment(block, point)

def find_xeqc_point_intersection(block):
    y2 = ((block.pivot - block.reference_start) * block.slope) + block.query_start
    return (block.pivot, y2)

def reflect_line_segment(block, point):
    slope = -1/block.slope
    x1, y1 = point
    block.y2_reflection = ((block.reference_end - x1) * slope) + y1
    block.y1_reflection = ((block.reference_start - x1) * slope) + y1
    block.reflected = True

def find_pivots(alignments, trend):
    for aln in alignments:
        if aln.slope != trend:  # This is an inversion
            minx = min(block.reference_start for block in aln.alignment_blocks)
            maxx = max(block.reference_start for block in aln.alignment_blocks)
            midx = (maxx + minx) / 2

            for block in aln.alignment_blocks:
                block.pivot = midx
