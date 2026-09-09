'''
Docstring for orchestration.call
Orchestrate the necessary steps for variant calling
'''

from svmu2.orchestration.synteny import run_synteny
from svmu2.core.classify import create_domain_range_trees, extract_collinear_gap_segments_from_path, call_all_inversions
from svmu2.IO.vcf import write_vcf

def call(args):
    alns = run_synteny(args)
    SVs = []
    
    for _, alignment in alns.items():
        domain_tree, range_tree = create_domain_range_trees(alignment.alignment_blocks)
        INDELS = extract_collinear_gap_segments_from_path(
            alignment.primary_synteny_blocks,
            alignment.reference,
            domain_tree,
            range_tree,
            write_bnds=args.write_bnds,
        )
        INVERSIONS = call_all_inversions(
            alignment.final_path_segments,
            alignment.primary_synteny_blocks,
            alignment.slope,
        )
        SVs.extend(INDELS + INVERSIONS)

    return alns, SVs

def run_call(args):
    _, SVs = call(args)
    write_vcf(SVs, output_path=args.out, sample=args.sample)
