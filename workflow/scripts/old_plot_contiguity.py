#!/usr/bin/env python

import argparse
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

def get_cumsum_curve(fasta):
    """Extracts lengths, sorts descending, and returns cumulative sum in bp."""
    try:
        lengths = sorted(
            [len(r.seq) for r in SeqIO.parse(fasta, "fasta")],
            reverse=True
        )
        if not lengths:
            return None
        return np.cumsum(lengths)
    except Exception as e:
        print(f"[ERROR] Could not parse {fasta}: {e}")
        return None

def plot_contiguity(
    results_dir,
    output_file,
    reference_fastas=None,
    strains=None,
    assembly_pattern="curated_assembly/*.fasta",
    max_contigs=20
):
    results_dir = Path(results_dir)
    
    # Use a built-in matplotlib style that doesn't require extra libs
    plt.style.use('bmh') 
    fig, ax = plt.subplots(figsize=(10, 7))

    # ---------------------------------------------------------
    # 1. Plot References (Grayscale / Dashed)
    # ---------------------------------------------------------
    if reference_fastas:
        for ref_path in reference_fastas:
            ref_path = Path(ref_path)
            curve = get_cumsum_curve(ref_path)
            if curve is None: continue
            
            x_vals = np.arange(1, len(curve) + 1)
            ax.step(x_vals, curve, color='#7f8c8d', alpha=0.7, 
                    linewidth=2, linestyle='--', where='post',
                    label=f"{ref_path.stem} (Ref)", zorder=2)

    # ---------------------------------------------------------
    # 2. Plot Strains (Color / Solid)
    # ---------------------------------------------------------
    if strains is None:
        strains = sorted([p.name for p in results_dir.iterdir() 
                         if p.is_dir() and p.name != "aggregate"])

    # Use a standard matplotlib colormap (Tableau Colors are high contrast)
    colors = plt.cm.tab10(np.linspace(0, 1, len(strains)))

    for i, strain in enumerate(strains):
        fasta_files = list((results_dir / strain).glob(assembly_pattern))
        if not fasta_files:
            print(f"[WARN] No FASTA found for {strain}")
            continue
        
        curve = get_cumsum_curve(fasta_files[0])
        if curve is None: continue

        x_vals = np.arange(1, len(curve) + 1)

        ax.step(x_vals, curve, color=colors[i], linewidth=2.5, 
                where='post', label=strain, alpha=1.0, zorder=3)

    # ---------------------------------------------------------
    # Formatting
    # ---------------------------------------------------------
    ax.set_xlim(1, max_contigs)
    ax.set_xticks(np.arange(1, max_contigs + 1))
    
    # Scientific notation for large BP values is usually cleaner
    ax.ticklabel_format(style='plain', axis='y') 
    
    ax.set_ylabel("Cumulative Assembly Size (bp)", fontsize=12, fontweight='bold')
    ax.set_xlabel(f"Contig Rank (Top {max_contigs})", fontsize=12, fontweight='bold')
    ax.set_title("Assembly Contiguity Curve", fontsize=14, pad=15, fontweight='bold')
    
    # Legend outside to avoid overlapping data
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", borderaxespad=0, frameon=True)

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
    parser.add_argument("--max-contigs", type=int, default=20)

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
