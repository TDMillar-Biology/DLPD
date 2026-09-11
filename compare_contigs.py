#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def contig_ids(fasta):
    return {record.id for record in SeqIO.parse(fasta, "fasta")}


def agp_components(agp):
    """
    Return the set of component contig IDs used in an AGP file.

    Standard AGP columns:
      1 object
      2 object_beg
      3 object_end
      4 part_number
      5 component_type
      6 component_id / gap_length
      ...
    """
    components = set()

    with open(agp) as fh:
        for line in fh:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 6:
                continue

            component_type = fields[4]

            # N/U are gaps. Everything else is a sequence component.
            if component_type not in {"N", "U"}:
                components.add(fields[5])

    return components


def print_group(name, values):
    print(f"## {name}: {len(values)}")
    for value in sorted(values):
        print(value)
    print()


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare contig IDs across original, curated, and scaffolded "
            "assemblies and inspect contig usage in an AGP file."
        )
    )

    parser.add_argument("--original", required=True)
    parser.add_argument("--curated", required=True)
    parser.add_argument("--scaffolded", required=True)
    parser.add_argument("--agp", required=True)

    args = parser.parse_args()

    original = contig_ids(args.original)
    curated = contig_ids(args.curated)
    scaffolded = contig_ids(args.scaffolded)
    agp = agp_components(args.agp)

    groups = {
        "original_only":
            original - curated - scaffolded,

        "curated_only":
            curated - original - scaffolded,

        "scaffolded_only":
            scaffolded - original - curated,

        "original_and_curated_only":
            (original & curated) - scaffolded,

        "original_and_scaffolded_only":
            (original & scaffolded) - curated,

        "curated_and_scaffolded_only":
            (curated & scaffolded) - original,

        "present_in_all_three":
            original & curated & scaffolded,
    }

    print("## Assembly counts")
    print(f"Original FASTA:        {len(original)}")
    print(f"Curated FASTA:         {len(curated)}")
    print(f"Scaffolded FASTA:      {len(scaffolded)}")
    print(f"AGP source components: {len(agp)}")
    print()

    print("## FASTA set comparison")
    print()

    for name, values in groups.items():
        print_group(name, values)

    print("## Scaffolding-specific comparison")
    print()

    print_group(
        "curated_contigs_used_in_scaffolds",
        curated & agp
    )

    print_group(
        "curated_contigs_not_used_in_scaffolds",
        curated - agp
    )

    print_group(
        "agp_components_not_found_in_curated_fasta",
        agp - curated
    )

    print_group(
        "agp_components_also_present_as_scaffolded_fasta_records",
        agp & scaffolded
    )


if __name__ == "__main__":
    main()
