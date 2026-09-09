'''
Orchestrate the necessary steps for the plot command
'''

from svmu2.models.line import build_alignment_primitives
from svmu2.orchestration.parse import run_parse
from svmu2.orchestration.synteny import resolve_synteny
from svmu2.IO.breakpoints import load_breakpoint_mapping

from pathlib import Path
import matplotlib.pyplot as plt
from svmu2.visualization.renderers import render_matplotlib


def render_plot(aln, ax=None):
    """
    Generates alignment primitives and renders them.
    If 'ax' is provided, it draws directly on the existing axis, making it 
    easy to embed in composite figures.
    """
    primitives = build_alignment_primitives(aln)
    xlabel = aln.reference
    ylabel = aln.query
    
    # Pass the axis down to the rendering function
    fig, ax = render_matplotlib(primitives, xlabel, ylabel, ax=ax)
    
    return fig, ax


def run_plot(args):
    """
    CLI entry point for plotting structural variants.
    """
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    all_alns, primary = run_parse(args)

    # 1. Synteny Integration
    # If the --synteny flag is passed, resolve primary synteny blocks
    if getattr(args, "synteny", False):
        breakpoints_path = getattr(args, "breakpoints", None) or getattr(args, "breakpoints_tsv", None)
        breakpoint_map = load_breakpoint_mapping(breakpoints_path) if breakpoints_path else {}
        
        # Populate the synteny fields on the primary alignments
        primary = resolve_synteny(primary, breakpoint_map)

    # 2. Target Selection
    targets = primary.values()
    if getattr(args, "reftarget", None):
        targets = [aln for aln in all_alns if aln.reference == args.reftarget]
        if not targets:
            raise ValueError(f"{args.reftarget} not a reference sequence in {args.alignment}")
    if getattr(args, "qrytarget", None):
        targets = [aln for aln in all_alns if aln.query == args.qrytarget]
        if not targets:
            raise ValueError(f"{args.qrytarget} not a query sequence in {args.alignment}")
    if getattr(args, "all", None):
        targets = all_alns

    # 3. Render and Save
    for aln in targets:
        # Exclusively render static matplotlib figures
        fig, ax = render_plot(aln)
        
        img_format = getattr(args, "img_format", "pdf")
        if img_format == "png":
            outpath = out_dir / f"{aln.reference}.{aln.query}.png"
        else:
            outpath = out_dir / f"{aln.reference}.{aln.query}.pdf"
        
        if getattr(args, "popup", None):
            plt.show()
        else:
            fig.savefig(outpath)
            
        plt.close(fig)