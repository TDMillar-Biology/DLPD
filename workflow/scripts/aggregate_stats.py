#!/usr/bin/env python

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


# -----------------------------
# Assembly stats
# -----------------------------
def compute_stats(fasta):
    lengths = sorted([len(r.seq) for r in SeqIO.parse(fasta, "fasta")], reverse=True)

    total = sum(lengths)
    n = len(lengths)

    cumsum = np.cumsum(lengths)
    idx = np.where(cumsum >= total / 2)[0][0]

    n50 = lengths[idx]
    l50 = idx + 1

    return total, n50, l50, n


# -----------------------------
# BUSCO parser
# -----------------------------
def parse_busco(summary_file):
    with open(summary_file) as f:
        for line in f:
            if line.strip().startswith("C:"):
                line = line.strip().replace(",", "")

                return (
                    float(line.split("C:")[1].split("%")[0]),
                    float(line.split("S:")[1].split("%")[0]),
                    float(line.split("D:")[1].split("%")[0]),
                    float(line.split("F:")[1].split("%")[0]),
                    float(line.split("M:")[1].split("%")[0]),
                )

    return [None] * 5


# -----------------------------
# QV parser
# -----------------------------
def parse_qv(qv_file):
    if not os.path.exists(qv_file):
        return None

    df = pd.read_csv(qv_file, sep="\t")
    return float(df["qv"].iloc[0])


# -----------------------------
# Core aggregation
# -----------------------------
def aggregate(results_dir, output_file, strains=None, busco_dirname="busco"):
    results_dir = Path(results_dir)

    # discover strains
    if strains is None:
        strains = [
            p.name for p in results_dir.iterdir()
            if p.is_dir() and p.name != "aggregate"
        ]

    rows = []

    for strain in strains:
        base = results_dir / strain

        if not base.exists():
            print(f"[WARN] Skipping missing strain: {strain}")
            continue

        # -------------------------
        # Contig-level stats
        # -------------------------
        contig_fastas = list(base.glob("assembly/*.fasta"))
        if not contig_fastas:
            print(f"[WARN] No assembly found for {strain}")
            continue

        contig_fasta = contig_fastas[0]
        total, cN50, cL50, n_contigs = compute_stats(contig_fasta)

        # -------------------------
        # Scaffold-level stats (future-proof)
        # -------------------------
        scaffold_fastas = list(base.glob("curated_assembly/*.fasta"))

        if scaffold_fastas:
            _, sN50, sL50, n_scaffolds = compute_stats(scaffold_fastas[0])
        else:
            sN50, sL50, n_scaffolds = [None] * 3

        # -------------------------
        # BUSCO
        # -------------------------
        busco_files = list(base.glob(f"{busco_dirname}/run_*/short_summary*.txt"))

        if busco_files:
            c, s, d, f_, m = parse_busco(busco_files[0])
        else:
            print(f"[WARN] No BUSCO summary for {strain}")
            c, s, d, f_, m = [None] * 5

        # -------------------------
        # QV
        # -------------------------
        qv_file = base / "qc" / f"{strain}_qc.tsv"
        qv = parse_qv(qv_file)

        rows.append({
            "strain": strain,
            "assembly_size": total,
            "contig_N50": cN50,
            "contig_L50": cL50,
            "scaffold_N50": sN50,
            "scaffold_L50": sL50,
            "n_contigs": n_contigs,
            "n_scaffolds": n_scaffolds,
            "QV": qv,
            "BUSCO_complete": c,
            "BUSCO_single": s,
            "BUSCO_dup": d,
            "BUSCO_frag": f_,
            "BUSCO_missing": m,
        })

    df = pd.DataFrame(rows)

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_file, sep="\t", index=False)

    print(f"[INFO] Wrote: {output_file}")


# -----------------------------
# CLI
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Aggregate assembly QC metrics across strains"
    )

    parser.add_argument(
        "--results",
        required=True,
        help="Results directory (e.g. results/)"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file"
    )

    parser.add_argument(
        "--strains",
        nargs="+",
        help="Optional subset of strains to include"
    )

    parser.add_argument(
        "--busco-dirname",
        default="busco",
        help="BUSCO directory name inside each strain (default: busco)"
    )

    args = parser.parse_args()

    aggregate(
        results_dir=args.results,
        output_file=args.output,
        strains=args.strains,
        busco_dirname=args.busco_dirname,
    )


if __name__ == "__main__":
    main()