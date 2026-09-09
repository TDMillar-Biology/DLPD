'''
Docstring for visualization.renderers
hold matplotlib and plotly renderers working on primitive line objects from
svmu.models.line:PrimitiveLine
'''

def render_matplotlib(primitives, xlabel, ylabel, ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))
    else:
        fig = ax.figure

    for p in primitives:
        ax.plot(p.x, p.y, color=p.color)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig, ax

def render_plotly(primitives, title):
    import plotly.graph_objects as go
    traces = []

    for p in primitives:
        traces.append(go.Scatter(
            x=p.x,
            y=p.y,
            mode='lines',
            line=dict(color=p.color),
            hovertext=p.hover_text,
            hoverinfo='text' if p.hover_text else 'skip'
        ))

    fig = go.Figure(data=traces)
    fig.update_layout(title=title)
    return fig


def render_sv(sv, ax):
    ax.plot(
        [sv.reference_start, sv.reference_start],
        [sv.query_start, sv.query_end],
        color="yellow",
    )
    ax.plot(
        [sv.reference_start, sv.reference_end],
        [sv.query_start, sv.query_start],
        color="yellow",
    )


def render_aln_block(block, ax, color="black"):
    return ax.plot(
        [block.reference_start, block.reference_end],
        [block.query_start, block.query_end],
        color=color,
    )


def render_alignment_blocks(aln, xlabel, ylabel):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 6))

    for block in aln.alignment_blocks:
        render_aln_block(block, ax)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig, ax


def plot_bounding_box(aln, ax, color="red", linestyle="--", linewidth=1.5):
    """
    Overlay the alignment bounding box on an existing dotplot axis.
    """
    if not hasattr(aln, "bounding_box"):
        raise RuntimeError("Bounding box not computed; call compute_bounding_box() first")

    bb = aln.bounding_box

    xs = [
        bb.bottom_left[0],
        bb.top_left[0],
        bb.top_right[0],
        bb.bottom_right[0],
        bb.bottom_left[0],
    ]
    ys = [
        bb.bottom_left[1],
        bb.top_left[1],
        bb.top_right[1],
        bb.bottom_right[1],
        bb.bottom_left[1],
    ]

    ax.plot(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth)


def plot_interactive_dotplot(best_delta, new_SVs=None, x_marker=None, output_path=None, auto_open=True):
    import plotly.graph_objects as go
    import plotly.io as pio

    plotly_data = []

    if new_SVs:
        for sv in new_SVs:
            color = "yellow" if sv.sv_type == "BND" else "orange"
            plotly_data.append(
                go.Scatter(
                    x=[sv.reference_start, sv.reference_end],
                    y=[sv.query_start, sv.query_end],
                    mode="lines",
                    line=dict(color=color),
                    hovertext=f"SV type: {sv.sv_type}",
                    hoverinfo="text",
                )
            )

    for block in best_delta.alignment_blocks:
        block_color = (
            "grey" if (not block.is_repeat and not block.part_of_primary_synteny) else
            "blue" if block.part_of_primary_synteny else
            "black"
        )
        hover_text = (
            f"Index: {block.index}<br>"
            f"Ref: {block.reference_start}-{block.reference_end}<br>"
            f"Query: {block.query_start}-{block.query_end}<br>"
            f"Repeat: {block.is_repeat}<br>"
            f"Primary Synteny: {block.part_of_primary_synteny}<br>"
            f"Reflected: {getattr(block, 'reflected', False)}"
        )
        plotly_data.append(
            go.Scatter(
                x=[block.reference_start, block.reference_end],
                y=[block.query_start, block.query_end],
                mode="lines",
                line=dict(color=block_color),
                hovertext=hover_text,
                hoverinfo="text",
            )
        )

        if block.reflected:
            reflect_hover_text = (
                f"Index: {block.index}<br>"
                f"Ref: {block.reference_start}-{block.reference_end}<br>"
                f"Query: {block.y1_reflection}-{block.y2_reflection}<br>"
                f"Repeat: {block.is_repeat}<br>"
                f"Primary Synteny: {block.part_of_primary_synteny}<br>"
                f"Reflected: {getattr(block, 'reflected', False)}"
            )
            plotly_data.append(
                go.Scatter(
                    x=[block.reference_start, block.reference_end],
                    y=[block.y1_reflection, block.y2_reflection],
                    mode="lines",
                    line=dict(color="green"),
                    hovertext=reflect_hover_text,
                    hoverinfo="text",
                )
            )

    if x_marker is not None:
        miny = min(b.query_start for b in best_delta.alignment_blocks)
        maxy = max(b.query_end for b in best_delta.alignment_blocks)
        plotly_data.append(
            go.Scatter(
                x=[x_marker, x_marker],
                y=[miny, maxy],
                mode="lines",
                line=dict(color="red", dash="dash"),
                name="X Marker",
                hoverinfo="skip",
            )
        )

    layout = go.Layout(
        title=f"Interactive Dotplot: {best_delta.reference}",
        xaxis=dict(title="Reference"),
        yaxis=dict(title="Query"),
        showlegend=False,
    )

    plotly_fig = go.Figure(data=plotly_data, layout=layout)
    if output_path is None:
        output_path = f"{best_delta.reference}_interactive.html"
    pio.write_html(plotly_fig, file=output_path, auto_open=auto_open)
    return plotly_fig
