#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from Bio import SeqIO


def fasta_ids(fasta):
    return [record.id for record in SeqIO.parse(fasta, "fasta")]


def infer_original_contig(curated_id, original_ids):
    """
    Map a curated contig back to its original parent.

    Exact match is preferred.

    Otherwise, if curated contig ends in _<integer>,
    repeatedly strip the terminal suffix and check against
    original contig IDs.

    Example:
        ptg000001l_1 -> ptg000001l
    """
    original_set = set(original_ids)

    if curated_id in original_set:
        return curated_id

    candidate = curated_id

    while True:
        match = re.match(r"^(.*)_([0-9]+)$", candidate)

        if not match:
            break

        candidate = match.group(1)

        if candidate in original_set:
            return candidate

    return None


def parse_agp(agp):
    """
    Parse standard AGP and return:

        component_id -> scaffold/object name

    Gap rows (N/U) are ignored.

    AGP:
      col 1 = object/scaffold
      col 5 = component type
      col 6 = component ID
    """
    component_to_scaffold = {}

    with open(agp) as fh:
        for line in fh:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 6:
                continue

            scaffold = fields[0]
            component_type = fields[4]

            if component_type in {"N", "U"}:
                continue

            component = fields[5]

            if component in component_to_scaffold:
                previous = component_to_scaffold[component]

                if previous != scaffold:
                    raise ValueError(
                        f"AGP component {component} occurs in multiple "
                        f"scaffolds: {previous}, {scaffold}"
                    )

            component_to_scaffold[component] = scaffold

    return component_to_scaffold


def build_mapping(original_fasta, curated_fasta, agp):
    original_ids = fasta_ids(original_fasta)
    curated_ids = fasta_ids(curated_fasta)

    original_set = set(original_ids)
    component_to_scaffold = parse_agp(agp)

    rows = []

    # ------------------------------------------------------
    # Curated contigs/fragments
    # ------------------------------------------------------

    for curated in curated_ids:
        original = infer_original_contig(curated, original_ids)

        if curated in component_to_scaffold:
            status = "scaffolded"
            scaffold = component_to_scaffold[curated]
        else:
            status = "unplaced"
            scaffold = "Unplaced"

        rows.append({
            "original_contig": original,
            "curated_contig": curated,
            "status": status,
            "scaffold": scaffold,
        })

    df = pd.DataFrame(rows)

    # ------------------------------------------------------
    # Detect original contigs that disappeared completely
    # ------------------------------------------------------

    represented_originals = set(
        x for x in df["original_contig"]
        if pd.notna(x)
    )

    missing_originals = original_set - represented_originals

    return df, missing_originals


def draw_ribbon(
    ax,
    x0,
    x1,
    y0_low,
    y0_high,
    y1_low,
    y1_high,
    alpha=0.45,
):
    dx = (x1 - x0) * 0.45

    verts = [
        (x0, y0_low),

        (x0 + dx, y0_low),
        (x1 - dx, y1_low),
        (x1, y1_low),

        (x1, y1_high),

        (x1 - dx, y1_high),
        (x0 + dx, y0_high),
        (x0, y0_high),

        (x0, y0_low),
    ]

    codes = [
        Path.MOVETO,

        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,

        Path.LINETO,

        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,

        Path.CLOSEPOLY,
    ]

    patch = PathPatch(
        Path(verts, codes),
        alpha=alpha,
        linewidth=0,
    )

    ax.add_patch(patch)


def calculate_positions(groups, gap=0.5):
    """
    Given:
        {
            node_name: number_of_items,
            ...
        }

    return:
        {
            node_name: (y0, y1)
        }
    """
    positions = {}

    y = 0

    for name, count in groups.items():
        positions[name] = (y, y + count)
        y += count + gap

    return positions

def is_derivative(curated_id, original_id):
    """
    True if curated_id is a split/derived descendant rather than an unchanged original.
    """
    if original_id is None:
        return False
    return curated_id != original_id

def make_alluvial(df, output, title=None):
    """
    Four layers:

        Original
        Curated
        Placement
        Scaffold

    Labeling behavior:
      - Original: label all, right-justified to the left of the node
      - Curated: label only derived contigs (e.g. ptg000001l_1)
      - Right side: label final curated contig endpoints
    """

    df = df.copy()

    df["is_derivative"] = [
        is_derivative(c, o)
        for c, o in zip(df["curated_contig"], df["original_contig"])
    ]

    df = df.sort_values(
        ["original_contig", "curated_contig", "status", "scaffold"],
        na_position="last"
    ).reset_index(drop=True)

    # ------------------------------------------------------
    # Node definitions
    # ------------------------------------------------------

    original_counts = (
        df.dropna(subset=["original_contig"])
          .groupby("original_contig")
          .size()
          .to_dict()
    )

    curated_nodes = {x: 1 for x in df["curated_contig"]}

    placement_nodes = {
        "scaffolded": int((df["status"] == "scaffolded").sum()),
        "unplaced": int((df["status"] == "unplaced").sum()),
    }

    scaffold_order = ["X", "2L", "2R", "3L", "3R", "4"]

    observed_scaffolds = list(df["scaffold"].unique())

    final_scaffolds = [x for x in scaffold_order if x in observed_scaffolds]
    extra_scaffolds = sorted(
        x for x in observed_scaffolds
        if x not in scaffold_order and x != "Unplaced"
    )
    final_scaffolds.extend(extra_scaffolds)
    if "Unplaced" in observed_scaffolds:
        final_scaffolds.append("Unplaced")

    scaffold_nodes = {
        scaffold: int((df["scaffold"] == scaffold).sum())
        for scaffold in final_scaffolds
    }

    gap = 0.35

    positions = [
        calculate_positions(original_counts, gap),
        calculate_positions(curated_nodes, gap),
        calculate_positions(placement_nodes, gap),
        calculate_positions(scaffold_nodes, gap),
    ]

    x_positions = [0, 1.5, 3.0, 4.5]
    node_width = 0.10

    fig_height = max(8, len(df) * 0.18)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    stage_names = ["Original", "Curated", "Placement", "Scaffold"]

    # ------------------------------------------------------
    # Draw nodes
    # ------------------------------------------------------

    for stage_idx, stage_positions in enumerate(positions):
        x = x_positions[stage_idx]

        for name, (y0, y1) in stage_positions.items():
            rect = Rectangle(
                (x - node_width / 2, y0),
                node_width,
                y1 - y0,
                alpha=0.8,
            )
            ax.add_patch(rect)

            # Node labels only for placement + scaffold stages
            if stage_idx >= 2:
                ax.text(
                    x,
                    (y0 + y1) / 2,
                    str(name),
                    ha="center",
                    va="center",
                    fontsize=8,
                )

    # ------------------------------------------------------
    # Offsets for ribbon packing
    # ------------------------------------------------------

    offsets = []
    for stage_positions in positions:
        offsets.append({name: y0 for name, (y0, y1) in stage_positions.items()})

    # Keep track of centers for labels
    original_label_pos = {}
    curated_label_pos = {}
    final_label_pos = {}

    # ------------------------------------------------------
    # Original -> Curated
    # ------------------------------------------------------

    for _, row in df.iterrows():
        original = row["original_contig"]
        curated = row["curated_contig"]

        if pd.isna(original):
            continue

        y0 = offsets[0][original]
        y1 = offsets[1][curated]

        draw_ribbon(
            ax,
            x_positions[0] + node_width / 2,
            x_positions[1] - node_width / 2,
            y0, y0 + 1,
            y1, y1 + 1,
        )

        # Save label positions
        original_label_pos.setdefault(original, []).append(y0 + 0.5)
        curated_label_pos[curated] = y1 + 0.5

        offsets[0][original] += 1
        offsets[1][curated] += 1

    # Reset curated offsets for next stage
    offsets[1] = {name: y0 for name, (y0, y1) in positions[1].items()}

    # ------------------------------------------------------
    # Curated -> Placement
    # ------------------------------------------------------

    for _, row in df.iterrows():
        curated = row["curated_contig"]
        placement = row["status"]

        y0 = offsets[1][curated]
        y1 = offsets[2][placement]

        draw_ribbon(
            ax,
            x_positions[1] + node_width / 2,
            x_positions[2] - node_width / 2,
            y0, y0 + 1,
            y1, y1 + 1,
        )

        offsets[1][curated] += 1
        offsets[2][placement] += 1

    # Reset placement offsets
    offsets[2] = {name: y0 for name, (y0, y1) in positions[2].items()}

    # ------------------------------------------------------
    # Placement -> Scaffold
    # ------------------------------------------------------

    for _, row in df.iterrows():
        placement = row["status"]
        scaffold = row["scaffold"]
        curated = row["curated_contig"]

        y0 = offsets[2][placement]
        y1 = offsets[3][scaffold]

        draw_ribbon(
            ax,
            x_positions[2] + node_width / 2,
            x_positions[3] - node_width / 2,
            y0, y0 + 1,
            y1, y1 + 1,
        )

        final_label_pos[curated] = y1 + 0.5

        offsets[2][placement] += 1
        offsets[3][scaffold] += 1

    # ------------------------------------------------------
    # Layer titles
    # ------------------------------------------------------

    for x, name in zip(x_positions, stage_names):
        ax.text(
            x,
            -1.5,
            name,
            ha="center",
            va="bottom",
            fontsize=13,
            fontweight="bold",
        )

    if title is None:
        title = "Contig flow through curation and scaffolding"

    ax.set_title(title, fontsize=15, pad=20)

    # ------------------------------------------------------
    # Left labels: original contigs, right-justified
    # ------------------------------------------------------

    original_x = x_positions[0] - node_width / 2 - 0.04

    for original, ys in original_label_pos.items():
        y = sum(ys) / len(ys)
        ax.text(
            original_x,
            y,
            original,
            ha="right",
            va="center",
            fontsize=7,
        )

    # ------------------------------------------------------
    # Curated labels: derivatives only
    # ------------------------------------------------------

    curated_x = x_positions[1] - node_width / 2 - 0.03

    for _, row in df.iterrows():
        if not row["is_derivative"]:
            continue

        curated = row["curated_contig"]
        y = curated_label_pos[curated]

        ax.text(
            curated_x,
            y,
            curated,
            ha="right",
            va="center",
            fontsize=6,
        )

    # ------------------------------------------------------
    # Far-right labels: final curated contig placement
    # ------------------------------------------------------

    final_x = x_positions[3] + node_width / 2 + 0.04

    for _, row in df.iterrows():
        curated = row["curated_contig"]
        y = final_label_pos[curated]

        ax.text(
            final_x,
            y,
            curated,
            ha="left",
            va="center",
            fontsize=6,
        )

    max_y = max(
        max(y1 for y0, y1 in stage.values())
        for stage in positions if stage
    )

    ax.set_xlim(-0.9, 5.3)
    ax.set_ylim(max_y + 1, -2)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Trace contigs from original assembly through "
            "curation and scaffolding."
        )
    )

    parser.add_argument(
        "--original",
        required=True,
    )

    parser.add_argument(
        "--curated",
        required=True,
    )

    parser.add_argument(
        "--agp",
        required=True,
    )

    parser.add_argument(
        "--out-plot",
        required=True,
    )

    parser.add_argument(
        "--out-table",
        required=True,
    )

    parser.add_argument(
        "--title",
        default=None,
    )

    args = parser.parse_args()

    df, missing_originals = build_mapping(
        args.original,
        args.curated,
        args.agp,
    )

    df.to_csv(
        args.out_table,
        sep="\t",
        index=False,
    )

    make_alluvial(
        df,
        args.out_plot,
        args.title,
    )

    print(f"Original contigs represented: {df['original_contig'].nunique()}")
    print(f"Curated contigs:              {len(df)}")
    print(
        f"Scaffolded contigs:           "
        f"{(df['status'] == 'scaffolded').sum()}"
    )
    print(
        f"Unplaced contigs:             "
        f"{(df['status'] == 'unplaced').sum()}"
    )

    if missing_originals:
        print()
        print("WARNING: original contigs with no curated descendant:")

        for contig in sorted(missing_originals):
            print(f"  {contig}")

    print()
    print(f"Wrote table: {args.out_table}")
    print(f"Wrote plot:  {args.out_plot}")


if __name__ == "__main__":
    main()