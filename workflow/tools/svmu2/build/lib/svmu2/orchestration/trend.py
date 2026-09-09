"""
svmu2.orchestration.trend

Debugging tool to visualize the cumulative sum of alignment blocks and 
verify the partition boundaries for structural variant calling.
"""

from svmu2.orchestration.synteny import run_synteny
from svmu2.visualization.trend import plot_alignment_dotplot, plot_signals

def run_trend(args):
    # 1. Parse alignments, load breakpoints, and resolve the global synteny 
    primary_alignments = run_synteny(args)

    # 2. Iterate through the resolved alignments and generate the debug plots
    for aln in primary_alignments.values():
        
        # Plot the 1D cumulative sum signals
        plot_signals(aln)
        
        # Plot the 2D dotplot with the partition boundaries overlaid
        plot_alignment_dotplot(
            aln, 
            ax=None, 
            vlines=aln.vlines, 
            hlines=aln.hlines
        )
