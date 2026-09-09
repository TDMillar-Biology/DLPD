'''
Docstring for models.line
primitive line dataclass for plotting
'''

from dataclasses import dataclass

@dataclass
class LinePrimitive:
    x: tuple
    y: tuple
    color: str
    hover_text: str | None = None

def build_alignment_primitives(aln, SVs=None, x_marker=None):
    primitives = []

    for block in aln.alignment_blocks:
        primitives.append(
            LinePrimitive(
                x=(block.reference_start, block.reference_end),
                y=(block.query_start, block.query_end),
                color="black",
                hover_text = f"({block.reference_start}, {block.query_start}) - ({block.reference_end}, {block.query_end})"
            )   
        )
    if aln.primary_synteny_blocks:
        for block in aln.primary_synteny_blocks:
            primitives.append(
                LinePrimitive(
                    x=(block.reference_start, block.reference_end),
                    y=(block.query_start, block.query_end),
                    color="blue",
                    hover_text = f"({block.reference_start}, {block.query_start}) - ({block.reference_end}, {block.query_end})"
                )   
            )
    if SVs:
        for SV in SVs:
            primitives.append(
                LinePrimitive(
                    x=(SV.reference_start, SV.reference_end),
                    y=(SV.query_start, SV.query_end),
                    color="orange"
                )   
            )

    return primitives
