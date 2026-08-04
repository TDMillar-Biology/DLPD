#!/usr/bin/env python

import argparse
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

def get_contiguity_data(fasta):
    """Returns sorted lengths and cumulative sum."""
    try:
        lengths = sorted(
            [len(r.seq) for r in SeqIO.parse(fasta, "fasta")],
            reverse=True
        )
        if not lengths:
            return None, None
        return np.array(lengths), np.cumsum(lengths)
    except Exception as e:
        print(f"[ERROR] Could not parse {fasta}: {e}")
        return None, None

def plot_contiguity(
    results_dir,
    output_file,
    reference_fastas=None,
    strains=None,
    assembly_pattern="curated_assembly/*.fasta",
    genome_size=180000000  # Default for Drosophila (180Mb)
):
    results_dir = Path(results_dir)
    
    # Professional styling: White background, thin grid
    plt.style.use('seaborn-v0_8-whitegrid') 
    fig, ax = plt.subplots(figsize=(9, 6))

    # 1. References (Gray lines)
    if reference_fastas:
        for ref_path in reference_fastas:
            ref_path = Path(ref_path)
            lengths, cumsum = get_contiguity_data(ref_path)
            if cumsum is None: continue
            
            # X-axis is the index of contigs, but we will focus on the curve shape
            x_pct = (cumsum / genome_size) * 100
            ax.plot(x_pct, lengths / 1e6, color='#95a5a6', alpha=0.6, 
                    linestyle='--', linewidth=1.5, label=f"{ref_path.stem} (Ref)")

    # 2. Strains (Primary Data)
    if strains is None:
        strains = sorted([p.name for p in results_dir.iterdir() 
                         if p.is_dir() and p.name != "aggregate"])

    # High-quality color palette
    colors = plt.cm.turbo(np.linspace(0.1, 0.9, len(strains)))

    for i, strain in enumerate(strains):
        fasta_files = list((results_dir / strain).glob(assembly_pattern))
        if not fasta_files: continue
        
        lengths, cumsum = get_contiguity_data(fasta_files[0])
        if cumsum is None: continue

        # The 'NG' style curve: Y = contig length, X = cumulative % of genome
        x_pct = (cumsum / genome_size) * 100
        
        # We plot the 'step' to show where each contig ends
        ax.step(x_pct, lengths / 1e6, color=colors[i], linewidth=2, 
                where='post', label=strain, alpha=0.8, zorder=3)

    # 3. Formatting to match Nature/Science standards
    ax.set_xlabel("Percentage of Genome Size (%)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Contig Length (Mb)", fontsize=12, fontweight='bold')
    ax.set_title("Contiguity (NG) Curve", fontsize=14, pad=15, fontweight='bold')

    # Add N50 reference line
    ax.axvline(50, color='black', linestyle=':', alpha=0.4, zorder=1)
    ax.text(51, ax.get_ylim()[1]*0.9, 'N50', fontsize=10, color='black', alpha=0.6)

    ax.set_xlim(0, 105) # Allow slightly over 100% for assemblies larger than ref
    ax.set_ylim(0, None)
    
    # Grid and Legend
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", frameon=False)

    plt.tight_layout()
    
    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference-fastas", nargs="*", default=[])
    parser.add_argument("--strains", nargs="*")
    parser.add_argument("--pattern", default="curated_assembly/*.fasta") 
    parser.add_argument("--genome-size", type=int, default=180000000)

    args = parser.parse_args()

    plot_contiguity(
        results_dir=args.results,
        output_file=args.output,
        reference_fastas=args.reference_fastas,
        strains=args.strains,
        assembly_pattern=args.pattern,
        genome_size=args.genome_size
    )

if __name__ == "__main__":
    main()
