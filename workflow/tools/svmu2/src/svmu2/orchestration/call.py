'''
Docstring for orchestration.call
Orchestrate the necessary steps for variant calling
'''

from svmu2.orchestration.synteny import run_synteny
from svmu2.core.classify import create_domain_range_trees, extract_collinear_gap_segments_from_path, call_all_inversions
from svmu2.IO.vcf import write_vcf

class IntraChainVariant:
    """Lightweight object to perfectly mimic a DotPlotLineSegment for write_vcf"""
    def __init__(self, sv_dict):
        self.chrom = sv_dict.get("chrom")
        self.sv_type = sv_dict.get("svtype")
        
        self.reference_start = sv_dict.get("pos")
        self.reference_end = sv_dict.get("end")

        # write_vcf uses (query_end - query_start) to calculate INS length.
        # We can mock this by setting start to 0 and end to the actual length.
        svlen = sv_dict.get("svlen")
        self.query_start = 0
        self.query_end = svlen
        
        # Default structural attributes expected by write_vcf logic
        self.theta = 0  # Setting to 0 bypasses the 45-degree filter and inversion logic
        self.range_partners = None
        self.domain_partners = None
        self.event_ID = f"intra_{self.reference_start}_{self.sv_type}"

def call(args):
    alns = run_synteny(args)
    SVs = []
    
    for _, alignment in alns.items():
        domain_tree, range_tree = create_domain_range_trees(alignment.alignment_blocks)
        
        # 1. Call inter-chain collinear gaps (Large INDELs)
        INDELS = extract_collinear_gap_segments_from_path(
            alignment.primary_synteny_blocks,
            alignment.reference,
            domain_tree,
            range_tree,
            write_bnds=args.write_bnds,
        )
        
        # 2. Call inter-chain inversions
        INVERSIONS = call_all_inversions(
            alignment.final_path_segments,
            alignment.primary_synteny_blocks,
            alignment.slope,
        )
        SVs.extend(INDELS + INVERSIONS)

        # 3. Call intra-chain micro-indels (if flag is thrown)
        if getattr(args, 'include_intra', False):
            print('yes')
            intra_variants = []
            if alignment.primary_synteny_blocks:
                for block in alignment.primary_synteny_blocks:
                    print(block.indel_map)
                    # Extract the dictionary list from the indel_map
                    block_indels = block.call_intra_chain_indels()
                    
                    # Wrap them in the object and append
                    for sv_dict in block_indels:
                        intra_variants.append(IntraChainVariant(sv_dict))
            print(len(intra_variants))
            SVs.extend(intra_variants)

    return alns, SVs

def run_call(args):
    _, SVs = call(args)
    write_vcf(SVs, output_path=args.out, sample=args.sample)