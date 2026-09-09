'''
Trevor Millar
Contains data structures needed to represent alignments as collections of line segments on 
a cartesian plane

These are frozen representations of state, no algos found here
'''

from dataclasses import dataclass
from math import sqrt

from intervaltree import IntervalTree


@dataclass(frozen=True)
class LineSegment:
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int

@dataclass(frozen=True)
class BoundingBox:
    ''' 
    Holds the most extreme coordinates (x,y) of rectangle containg the minimal 
    space required to hold every alignment block in a 2D plane. Is immutable
    '''
    bottom_left: tuple
    top_left: tuple
    top_right: tuple
    bottom_right: tuple

class Alignment:
    ''' 
    This is the full alignment reported by the aligner for a single chromosome
    It is composed of multiple individual line segments denoting synteny
    '''
    def __init__(self, reference, query, reference_length, query_length, parent=None):
        self.reference = reference
        self.query = query
        self.reference_length = int(reference_length)
        self.query_length = int(query_length)
        self.parent = parent
        self.alignment_blocks = None
        self.zero_coverage_sites = None
        self.covered_sites = None
        self.coverage_by_position = []
        self.traversal_len = None
        self.origin = (0,0) ## poised for removal
        self.slope = None ## poised for removal
        self.primary_synteny_blocks = None # Djikstra traversal defines these
        self.primary_synteny_tree = None # main diagonal / orthologous
        self.exhaustive_tree = None   # all homologous blocks
        self.source = None
        self.sink = None
        self.score = None
        self.trace = None
        self.segment_start = None
        self.segment_end = None
        self.x_cumsum = None
        self.y_cumsum = None
        self.vlines = []
        self.hlines = []
        self.local_optimal_traversals = None
        self.final_path_segments = None

    def __str__(self):
        return f'{self.reference} by {self.query} alignment, consisting of {0 if self.alignment_blocks is None else len(self.alignment_blocks)} alignment blocks'

    def add_alignment_block(self, line, indel_map, index, source_hash=None):
        if self.alignment_blocks is None:
            self.alignment_blocks = []
        self.alignment_blocks.append(AlignmentBlock(self.reference, *line, self.query, indel_map, index, source_hash))
        self.alignment_blocks.sort(key=lambda x: x.left_most) ## keep the alignment blocks sorted by leftmost position on reference

    def set_primary_synteny_blocks(self, blocks):
        """
        Register primary syntenic alignment blocks
        (e.g. result of Dijkstra traversal).
        """
        if blocks is None or len(blocks) == 0:
            print(f"Primary synteny block list is empty for {self.reference} {self.query}")
            #raise ValueError("Primary synteny block list is empty")

        self.primary_synteny_blocks = blocks
        self.primary_synteny_tree = None  # invalidate any stale tree

    def clone_with_blocks(self, blocks, segment_start=None, segment_end=None):
        """Return a shallow copy of this alignment limited to a block subset."""
        clone = Alignment(
            self.reference,
            self.query,
            self.reference_length,
            self.query_length,
            parent=self,
        )
        clone.alignment_blocks = list(blocks)
        clone.primary_synteny_blocks = [
            block for block in (self.primary_synteny_blocks or [])
            if block in blocks
        ]
        clone.primary_synteny_tree = None
        clone.exhaustive_tree = None
        clone.slope = self.slope
        clone.source = self.source if self.source in blocks else None
        clone.sink = self.sink if self.sink in blocks else None
        clone.origin = self.origin
        clone.segment_start = segment_start
        clone.segment_end = segment_end
        clone.local_optimal_traversals = None
        clone.final_path_segments = None
        return clone

    def calculate_coverage(self):
        ''' 
        Treat start,end of alignment as events with value (+1,-1)
        Find segments along ref length that have no alignments
        '''
        events = []
        for block in self.alignment_blocks:
            events.append((block.left_most, 1))
            events.append((block.right_most + 1, -1))
        
        events.sort()
        coverage = 0
        last_pos = 0
        self.zero_coverage_sites = 0
        
        for pos, change in events: # go through events, add to zero_coverage_sites the uncovered region
            if pos > last_pos and coverage == 0:
                self.zero_coverage_sites += pos - last_pos  # count uncovered region
            coverage += change
            last_pos = pos

        # Count trailing uncovered region
        if last_pos < self.reference_length:
            self.zero_coverage_sites += self.reference_length - last_pos

        self.covered_sites = self.reference_length - self.zero_coverage_sites

    def calculate_range_coverage(self):
            ''' 
            Treat start,end of alignment as events with value (+1,-1)
            Find segments along qry length that have no alignments
            '''
            events = []
            for block in self.alignment_blocks:
                events.append((block.bottom, 1))
                events.append((block.top + 1, -1))
            
            events.sort()
            coverage = 0
            last_pos = 0
            self.zero_coverage_sites = 0
            
            for pos, change in events: # go through events, add to zero_coverage_sites the uncovered region
                if pos > last_pos and coverage == 0:
                    self.zero_coverage_sites += pos - last_pos  # count uncovered region
                coverage += change
                last_pos = pos

            # Count trailing uncovered region
            if last_pos < self.query_length:
                self.zero_coverage_sites += self.query_length - last_pos

            self.covered_sites = self.query_length - self.zero_coverage_sites

    def compute_bounding_box(self):
        """
        Computes an immutable BoundingBox to the alignment object (see BoundingBox class definition)
        """
        if not self.alignment_blocks:
            raise ValueError("No alignment blocks provided.")

        minx = min(b.left_most for b in self.alignment_blocks)
        maxx = max(b.right_most for b in self.alignment_blocks)
        miny = min(b.bottom for b in self.alignment_blocks)
        maxy = max(b.top for b in self.alignment_blocks)

        bottom_left  = (minx, miny)
        top_left     = (minx, maxy)
        top_right    = (maxx, maxy)
        bottom_right = (maxx, miny)

        self.bounding_box = BoundingBox(
            bottom_left=bottom_left,
            top_left=top_left,
            top_right=top_right,
            bottom_right=bottom_right
        )
    @staticmethod
    def euclidean_distance(p1, p2):
        '''Need to be careful here about namespace, would be good to check, maybe define globally if possible'''
        if len(p1) != len(p2):
            raise ValueError("Points must have the same dimensionality")
        return sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))
    
    @staticmethod
    def x_distance(p1, p2):
        """
        Reference-space (x-axis) distance only.
        Used when reference order should dominate (primary synteny endpoints).
        """
        if len(p1) < 1 or len(p2) < 1:
            raise ValueError("Points must have at least one dimension")

        return abs(p2[0] - p1[0])

class AlignmentBlock:
    def __init__(self, reference, reference_start, reference_end, query_start, query_end, errors, similarity_errors, stop_codons, query , indel_map, index, source_hash = None):
        self.reference = reference
        self.query = query
        self.reference_start = int(reference_start)
        self.reference_end = int(reference_end)
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.errors = int(errors)
        self.similarity_errors = int(similarity_errors)
        self.stop_codons = int(stop_codons)
        self.slope = (self.query_end - self.query_start) / (self.reference_end - self.reference_start)
        self.rounded_slope = round(self.slope, 0)
        self.length = self.reference_end - self.reference_start # consider removal
        self.left_most = min(self.reference_start, self.reference_end)
        self.right_most = max(self.reference_start, self.reference_end)
        self.bottom = min(self.query_end, self.query_start)
        self.top = max(self.query_end, self.query_start)
        self.domain = (self.left_most, self.right_most)
        self.length_as_sole_coverage = 0
        self.y_intercept = self.query_start - (self.slope * self.reference_start) ## b = y - mx
        self.domain_length = self.reference_end - self.reference_start
        self.x_mid = (self.reference_start + self.reference_end) / 2
        self.y_mid = (self.query_start + self.query_end) / 2
        self.z_scaled_y_intercept = None
        self.start = (self.reference_start, self.query_start)
        self.end = (self.reference_end, self.query_end)
        self.left_neighbor = None
        self.right_neighbor = None
        self.shares_domain = []
        self.shares_range = []
        self.part_of_primary_synteny = False
        self.is_repeat = False
        self.repeat_content = 0
        self.index = index
        self.weighted_reference_length = self.slope * (self.right_most - self.left_most)
        self.weighted_query_length = self.slope * (self.top - self.bottom)
        self.y1_reflection = None
        self.y2_reflection = None
        self.pivot = self.y_mid
        self.reflected = False
        self.indel_map = indel_map
        self.segments = None
        self.segment_tree = None
        self.source_hash = source_hash

    def __str__(self):
        return f'{self.reference_start} {self.reference_end} {self.query_start} {self.query_end} {self.errors} {self.similarity_errors} {self.stop_codons} {self.slope}'

    def build_segments(self):
        """
        Derive ungapped syntenic segments from delta indel map.
        Effectively deprecated, kept in case I explore constituent line segments in the future
        """
        if self.segments is not None:
            return

        ref_pos = self.reference_start
        qry_pos = self.query_start
        self.segments = []

        for d in self.indel_map:
            if d == 0:
                break

            seg_len = abs(d)
            if seg_len > 1:
                self.segments.append(
                    LineSegment(
                        ref_start=ref_pos,
                        ref_end=ref_pos + seg_len,
                        query_start=qry_pos,
                        query_end=qry_pos + seg_len
                    )
                )

            ref_pos += seg_len
            qry_pos += seg_len

            if d > 0:
                ref_pos += 1
            else:
                qry_pos += 1

    def build_interval_tree(self):
        """Build reference-space interval tree over segments."""
        if self.segment_tree is not None:
            return

        if self.segments is None:
            raise RuntimeError("build_segments() must be called first")

        self.segment_tree = IntervalTree()
        for seg in self.segments:
            self.segment_tree[seg.ref_start:seg.ref_end] = seg
