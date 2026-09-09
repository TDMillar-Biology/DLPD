'''
Docstring for orchestration.synteny
Operate on alignment graph space to determine elements of primary synteny
'''

from svmu2.orchestration.parse import run_parse
from svmu2.core.synteny import resolve_alignment_synteny, build_primary_synteny_trees
from svmu2.IO.breakpoints import load_breakpoint_mapping

def resolve_synteny(primary_alignments, breakpoint_map=None):
    """
    Core synteny resolution logic. Resolves the primary synteny for given alignments.
    
    Args:
        primary_alignments: Dict of alignment objects to process
        breakpoint_map: Optional dict mapping (reference, query) pairs to breakpoint counts
    
    Returns:
        The input primary_alignments dict with synteny fields populated
    """
    if breakpoint_map is None:
        breakpoint_map = {}

    for aln_id, aln in primary_alignments.items():
        # Look up the specific alignment pair using the tuple key
        pair_key = (aln.reference, aln.query)
        
        # Default to 0 if the pair is not in the breakpoint map
        bkps = breakpoint_map.get(pair_key, 0) 
        
        resolve_alignment_synteny(aln, expected_breakpoints=bkps)
        build_primary_synteny_trees(aln)

    return primary_alignments

def run_synteny(args):
    """
    CLI entry point for the synteny module. Handles argument parsing and delegates
    to resolve_synteny() for the actual synteny resolution logic.
    """
    _, primary_alignments = run_parse(args)

    # Accept both the current CLI name and the older internal name.
    breakpoints_path = getattr(args, "breakpoints", None) or getattr(args, "breakpoints_tsv", None)

    # Load the mapping if provided, otherwise default to an empty dict
    breakpoint_map = load_breakpoint_mapping(breakpoints_path) if breakpoints_path else {}

    return resolve_synteny(primary_alignments, breakpoint_map)
