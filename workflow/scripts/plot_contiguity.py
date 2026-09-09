#!/usr/init/env python

import argparse
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def get_contiguity_data(fasta):
    """Returns sorted tuples of (length, index_number) and cumulative sums."""
    try:
        raw_lengths = sorted(
            [len(r.seq) for r in SeqIO.parse(fasta, "fasta")],
            reverse=True
        )
        if not raw_lengths:
            return None, None
        
        # Create tuples of (length, index_number) where index is 1-based
        contig_tuples = [(length, i + 1) for i, length in enumerate(raw_lengths)]
        
        lengths_arr = np.array([t[0] for t in contig_tuples])
        cumsum_arr = np.cumsum(lengths_arr)
        
        return contig_tuples, cumsum_arr
    except Exception as e:
        print(f"[ERROR] Could not parse {fasta}: {e}")
        return None, None

def plot_contiguity(
    results_dir,
    output_file,
    reference_fastas=None,
    strains=None,
    assembly_pattern="curated_assembly/*.fasta",
    max_contigs=100
):
    results_dir = Path(results_dir)
    
    # Professional styling: White background, thin grid
    plt.style.use('seaborn-v0_8-whitegrid') 
    fig, ax = plt.subplots(figsize=(9, 6))

    # ==========================================
    # PLUG IT IN RIGHT HERE (References Section)
    # ==========================================
    if reference_fastas:
        # Generates distinct shades, avoiding pure white or black, scaling to any number of references
        ref_colors = plt.cm.Greys(np.linspace(0.3, 0.8, len(reference_fastas)))
        ref_styles = ['--', '-.', ':', '-']
        
        for idx, ref_path in enumerate(reference_fastas):
            ref_path = Path(ref_path)
            contig_tuples, cumsum = get_contiguity_data(ref_path)
            if cumsum is None: continue
            
            x_indices = [t[1] for t in contig_tuples]
            y_cumsum_mb = cumsum / 1e6
            
            ax.plot(x_indices, y_cumsum_mb, color=ref_colors[idx], alpha=0.9, 
                    linestyle=ref_styles[idx % len(ref_styles)], linewidth=1.8, 
                    label=f"{ref_path.stem} (Ref)")

    # 2. Strains (Primary Data)
    if strains is None:
        strains = sorted([p.name for p in results_dir.iterdir() 
                         if p.is_dir() and p.name != "aggregate"])

    # High-quality color palette
    colors = plt.cm.turbo(np.linspace(0.1, 0.9, len(strains)))

    for i, strain in enumerate(strains):
        fasta_files = list((results_dir / strain).glob(assembly_pattern))
        if not fasta_files: continue
        
        contig_tuples, cumsum = get_contiguity_data(fasta_files[0])
        if cumsum is None: continue

        x_indices = [t[1] for t in contig_tuples]
        y_cumsum_mb = cumsum / 1e6
        
        ax.plot(x_indices, y_cumsum_mb, color=colors[i], linewidth=2, 
                linestyle='-', label=strain, alpha=0.8, zorder=3)

    # 3. Formatting to match standards
    ax.set_xlabel("Contig Index (Rank, log10)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Cumulative Assembly Length (Mb)", fontsize=12, fontweight='bold')
    ax.set_title("Cumulative Assembly Footprint", fontsize=14, pad=15, fontweight='bold')

    # Apply Log Scale to X-axis
    ax.set_xscale('log')
    ax.set_xlim(1, max_contigs)

    # Custom X-axis ticks and integer labels
    desired_ticks = [1, 5, 20, 100]
    ax.set_xticks(desired_ticks)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.ticklabel_format(style='plain', axis='x')

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
    parser.add_argument("--max-contigs", type=int, default = 100)
    args = parser.parse_args()

    plot_contiguity(
        results_dir=args.results,
        output_file=args.output,
        reference_fastas=args.reference_fastas,
        strains=args.strains,
        assembly_pattern=args.pattern,
        max_contigs=args.max_contigs
    )

if __name__ == "__main__":
    main()