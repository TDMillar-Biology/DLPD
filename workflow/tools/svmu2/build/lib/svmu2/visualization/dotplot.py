'''
Compatibility module for dotplot rendering helpers.

Rendering implementations live in svmu2.visualization.renderers.
'''

from svmu2.visualization.renderers import (
    plot_bounding_box,
    plot_interactive_dotplot,
    render_alignment_blocks,
    render_aln_block,
    render_sv,
)

__all__ = [
    "plot_bounding_box",
    "plot_interactive_dotplot",
    "render_alignment_blocks",
    "render_aln_block",
    "render_sv",
]
