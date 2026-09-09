'''
Docstring for orchestration.parse
Orchestrate parsing of delta file and selection of main alignments
Two layer structure here is purposeful -- keep CLI clean and API access
'''
from svmu2.IO.delta import parse_delta_file
from svmu2.IO.sam import parse_sam_file 
from svmu2.IO.paf import parse_paf_file 
from svmu2.core.selection import select_primary_alignments

def parse(alignment_file, fmt="delta"):
    # Dispatch to the correct parser based on the format flag
    if fmt == "delta":
        all_alns = parse_delta_file(alignment_file)
    elif fmt == "sam":
        all_alns = parse_sam_file(alignment_file)
    elif fmt == "paf":
        all_alns = parse_paf_file(alignment_file)

    else:
        raise ValueError(f"Unsupported alignment format: {fmt}")

    primary = select_primary_alignments(all_alns)

    if not primary: 
        raise ValueError("No primary alignments found")
    
    return all_alns, primary

def run_parse(args):
    # Retrieve the format flag, defaulting to 'delta' if not present
    fmt = getattr(args, 'format', 'delta')
    # Keeping args.delta for backwards compatibility with downstream modules
    return parse(args.alignment, fmt=fmt)
