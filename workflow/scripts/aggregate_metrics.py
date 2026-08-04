#!/usr/bin/env python

import argparse
import os
from pathlib import Path

import pandas as pd


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
# seqkit parser
# -----------------------------
def parse_seqkit_stats(stats_file):
    if not os.path.exists(stats_file):
        print(f"[WARN] Missing stats file: {stats_file}")
        return None, None, None

    df = pd.read_csv(stats_file, sep=r'\s+', thousands=",")
    df.columns = df.columns.str.strip()
    print("this is running")
    print(stats_file)
    print(df)
    row = df.iloc[0]

    return (
        int(row["sum_len"]),
        int(row["N50"]),
        int(row["num_seqs"]),
    )
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
        # Contig-level stats (seqkit)
        # -------------------------
        contig_stats_file = base / "qc" / f"{strain}.seqkit.stats.tsv"

        assembly_size, contig_N50, n_contigs = parse_seqkit_stats(contig_stats_file)

        if assembly_size is None:
            print(f"[WARN] Missing seqkit contig stats for {strain}")

        # -------------------------
        # Scaffold-level stats (future-proof)
        # -------------------------
        scaffold_stats_file = base / "qc" / f"{strain}.scaffold.stats.tsv"

        scaffold_size, scaffold_N50, n_scaffolds = parse_seqkit_stats(scaffold_stats_file)

        if scaffold_size is None:
            scaffold_N50, n_scaffolds = None, None

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

        # -------------------------
        # Collect row
        # -------------------------
        rows.append({
            "strain": strain,
            "assembly_size": assembly_size,
            "contig_N50": contig_N50,
            "scaffold_N50": scaffold_N50,
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
        description="Aggregate assembly QC metrics (seqkit-based)"
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
        help="Optional subset of strains"
    )

    parser.add_argument(
        "--busco-dirname",
        default="busco",
        help="BUSCO directory name (default: busco)"
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
