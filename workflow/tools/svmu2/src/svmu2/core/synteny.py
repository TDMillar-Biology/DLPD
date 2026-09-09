from svmu2.core.segmentation import (
    generate_alignment_partitions, 
    calculate_partition_bounds, 
    resolve_partition_orientations, 
    resolve_global_synteny, 
    unique_main_elements
)

from intervaltree import IntervalTree

def determine_global_slope(primary_synteny_blocks):
    """
    Determines the global orientation of the alignment by comparing the 
    query coordinates of the first and last blocks in the resolved path.
    
    Assumes blocks are pre-sorted by reference_start.

    This is a DUCT TAPE fix and needs to be made more robust
    """
    if not primary_synteny_blocks:
        return 1  # Safe default if path is empty

    first_block = primary_synteny_blocks[0]
    last_block = primary_synteny_blocks[-1]

    if last_block.query_start >= first_block.query_start:
        return 1
    else:
        return -1


def apply_final_path_state(aln, final_path_segments, primary_blocks):
    """
    Normalize alignment block state from the globally selected path only.

    This prevents provisional local partition winners from leaking into
    downstream call logic.
    """
    primary_indices = {block.index for block in primary_blocks}
    reflected_indices = {
        block.index
        for segment in final_path_segments
        if segment.orientation_sign == -1
        for block in segment.traversal_data.get("path_blocks", [])
    }

    for block in aln.alignment_blocks:
        block.part_of_primary_synteny = block.index in primary_indices
        block.reflected = block.index in reflected_indices

def resolve_alignment_synteny(aln, expected_breakpoints=0):
    """
    Computes the primary syntenic path and structural variants for a single alignment.
    """
    dim = expected_breakpoints + 1

    # 1. Partition the alignment space
    hlines, vlines = calculate_partition_bounds(aln, breakpoints=expected_breakpoints)

    aln.vlines = vlines
    aln.hlines = hlines
    
    partitions = generate_alignment_partitions(aln, vlines, hlines)

    # 2. Traverse and resolve the optimal path
    local_optimal_traversals = resolve_partition_orientations(aln, partitions)
    final_path = resolve_global_synteny(local_optimal_traversals, rows=dim, cols=dim)

    aln.local_optimal_traversals = local_optimal_traversals
    aln.final_path_segments = final_path

    # 3. Extract the unique main blocks
    main_elements = unique_main_elements(final_path)
    aln.set_primary_synteny_blocks(main_elements)
    apply_final_path_state(aln, final_path, main_elements)
    
    # Calculate the global slope dynamically from the optimal path
    aln.slope = determine_global_slope(main_elements)
    return aln

def build_reference_synteny_tree(alignment):
    """REFERENCE SPACE SYNTENY TREE"""
    if alignment.primary_synteny_blocks is None:
        raise RuntimeError("primary synteny blocks not set, use dijkstra traversal before this tool")
    if len(alignment.primary_synteny_blocks) == 0:
        #raise RuntimeError("primary synteny block list is empty")  
        return IntervalTree() # handling case of dijkstra failed, but we dont want to crash

    tree = IntervalTree()
    for block in alignment.primary_synteny_blocks:
        tree[block.left_most:block.right_most] = block

    alignment.reference_synteny_tree = tree
    return tree

def build_query_synteny_tree(alignment):
    """QUERY SPACE SYNTENY TREE"""
    if alignment.primary_synteny_blocks is None:
        raise RuntimeError("primary synteny blocks not set, use dijkstra traversal before this tool")
    if len(alignment.primary_synteny_blocks) == 0:
        #raise RuntimeError("primary synteny block list is empty")  
        return IntervalTree() # handling case of dijkstra failed, but we dont want to crash

    tree = IntervalTree()
    for block in alignment.primary_synteny_blocks:

        start = min(block.top, block.bottom)
        end = max(block.top, block.bottom)

        # IntervalTree uses half-open intervals [start, end).
        # Check to see if ends are inclusive
        
        tree[start:end] = block

    alignment.query_synteny_tree = tree
    return tree

def build_primary_synteny_trees(alignment):
    build_reference_synteny_tree(alignment)
    build_query_synteny_tree(alignment)
    