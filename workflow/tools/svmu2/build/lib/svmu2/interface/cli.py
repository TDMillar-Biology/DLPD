import argparse

def main():
    ########## ARGUMENT PARSING ##########
    parser = argparse.ArgumentParser(
        prog="svmu2",
        description="Structural Variants from MUmmer and associated tools"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    ##### PLOT #####
    plot = subparsers.add_parser("plot", help="Static Dotplot of genome to genome alignment")
    plot.add_argument("--alignment", '-d', required=True, help='Path to alignment file (.delta or .sam)')
    plot.add_argument("--format", choices=("delta", "sam", "paf"), default="delta", help="Format of the input alignment file")
    plot.add_argument("--out_dir", '-o', required=True, help='Out dir to write output image to')
    plot.add_argument("--strain", required=False)
    plot.add_argument("--contig", default=None)
    plot.add_argument("--reftarget", default=None, required=False, help = "plot all alignments involving reftarget")
    plot.add_argument("--qrytarget", default=None, required=False, help = "plot all alignments involving qrytarget")
    plot.add_argument("--all", action = 'store_true',default = None, required=False, help="Force plotting of all alignments")
    plot.add_argument("--popup", action = 'store_true', required = False, help="Pop up the figure as opposed to writing to file")
    plot.add_argument("--img-format", dest="img_format", choices=("pdf", "png"), default="pdf", help="Output image format")
    plot.add_argument('--synteny', action='store_true', help="Run synteny resolution to highlight the main diagonal.")
    plot.add_argument('--breakpoints', type=str, metavar="FILE", help="Path to a TSV mapping breakpoints for synteny blocks.")

    ##### PLOT-SV #####
    plot_sv = subparsers.add_parser("plot-sv", help="Dotplot with primary synteny path and SV break segments")
    plot_sv.add_argument("--alignment", '-d', required=True, help='Path to alignment file (.delta or .sam)')
    plot_sv.add_argument("--format", choices=("delta", "sam", "paf"), default="delta", help="Format of the input alignment file")
    plot_sv.add_argument("--out_dir", '-o', required=True, help='Out dir to write output image to')
    plot_sv.add_argument("--reftarget", default=None, required=False, help="Plot only the selected reference sequence")
    plot_sv.add_argument("--qrytarget", default=None, required=False, help="Plot only the selected query sequence")
    plot_sv.add_argument("--img-format", dest="img_format", choices=("pdf", "png"), default="pdf", help="Output image format")
    plot_sv.add_argument("--interactive", action='store_true', help="Write an interactive HTML plot instead of a static image")
    plot_sv.add_argument("--popup", action='store_true', required=False, help="Pop up the figure as opposed to writing to file")
    plot_sv.add_argument('--write-bnds', '-wb', action='store_true', help='Also plot BND breakend records for each SV')
    plot_sv.add_argument("--debug-trend", action='store_true', help="Write cumulative trend debug plots")
    plot_sv.add_argument("--debug-trend-dir", default=None, help="Output directory for cumulative trend debug plots")
    plot_sv.add_argument("--show-debug-trend", action='store_true', help="Display cumulative trend debug plots interactively")
    plot_sv.add_argument("--breakpoints", "-b", required=False, help="Path to curated breakpoint mapping TSV")

    ##### CALL #####
    call = subparsers.add_parser("call", help="Call SVs from genome to genome alignment. Write vcf")
    call.add_argument("--alignment", '-d', required=True, help='Path to alignment file (.delta or .sam)')
    call.add_argument("--format", choices=("delta", "sam", "paf"), default="delta", help="Format of the input alignment file")
    call.add_argument("--out", '-o', required=True, help='Path to write output vcf')
    call.add_argument('--sample', '-s', default='SAMPLE', help='Sample name for VCF output')
    call.add_argument('--write-bnds', '-wb', action='store_true', help='Also output BND breakends for each SV (default: False, niche use case, not recommended for most users)')
    call.add_argument("--debug-trend", action='store_true', help="Write cumulative trend debug plots")
    call.add_argument("--debug-trend-dir", default=None, help="Output directory for cumulative trend debug plots")
    call.add_argument("--show-debug-trend", action='store_true', help="Display cumulative trend debug plots interactively")
    call.add_argument("--breakpoints", "-b", required=False, help="Path to curated breakpoint mapping TSV")

    ##### INTERACTIVE #####
    interactive = subparsers.add_parser("interactive", help="Interactive Dotplot of genome to genome alignment")
    interactive.add_argument("--alignment", '-d', required=True, help='Path to alignment file (.delta or .sam)')
    interactive.add_argument("--format", choices=("delta", "sam", "paf"), default="delta", help="Format of the input alignment file")
    interactive.add_argument("--out_dir", '-o', required=True, help='Out dir to write output html to')
    interactive.add_argument("--popup", action='store_true')

    ##### TREND #####
    trend = subparsers.add_parser("trend", help = "Debugging function to visualize cumulative sum of weighted values of alignment blocks")
    trend.add_argument("--alignment", '-d', required=True, help='Path to alignment file (.delta or .sam)')
    trend.add_argument("--format", choices=("delta", "sam", "paf"), default="delta", help="Format of the input alignment file")
    trend.add_argument("--breakpoints", "-b", required=False, help="Path to curated breakpoint mapping TSV")
    
    ##### INSPECT #####
    inspect = subparsers.add_parser("inspect", help="Parse delta file headers to generate a breakpoint mapping TSV template")
    inspect.add_argument("--alignment", '-d', required=True, help='Path to alignment file (.delta or .sam)')
    inspect.add_argument("--format", choices=("delta", "sam", "paf"), default="delta", help="Format of the input alignment file")
    inspect.add_argument("--out", '-o', required=True, help='Path to write the output TSV template')
    inspect.add_argument("--force", '-f', action='store_true', help='Force overwrite if the output file already exists')

    ##### FILTER ##### 
    filter_cmd = subparsers.add_parser(
        "filter", help=("Filter SAM/BAM alignment to retain only primary syntenic 'main' diagonal elements (currently supported for SAM only)"),
    )
    filter_cmd.add_argument( "--alignment", "-a", dest="alignment", required=True, help="Path to alignment file (.sam or .bam)",)
    filter_cmd.add_argument("--format", choices=("sam", "paf",), default="paf", help=("Format of the input alignment file (paf / sam only for now)"),)
    filter_cmd.add_argument("--out", "-o", default=None, help=("Explicit path to write output filtered alignment. If omitted, defaults to {out_dir}/{prefix}{alignment_filename}"),)
    filter_cmd.add_argument("--prefix","-p",default="filtered_",help="Prefix for output file name if --out is not provided (default:'filtered_')",)
    filter_cmd.add_argument("--out_dir",default=None,help="Directory to write output file (defaults to input file's directory)",)
    filter_cmd.add_argument("--breakpoints", "-b", required=False, help="Path to curated breakpoint mapping TSV",)

    args = parser.parse_args()
    ########## COMMAND DISPATCH ##########

    if args.command == "plot":
        from svmu2.orchestration.plot import run_plot
        run_plot(args)

    elif args.command == "plot-sv":
        from svmu2.orchestration.plot_sv import run_plot_sv
        run_plot_sv(args)

    elif args.command == "call":
        from svmu2.orchestration.call import run_call
        run_call(args)

    elif args.command == "interactive":
        from svmu2.orchestration.plot import run_plot
        run_plot(args)

    elif args.command == "trend":
        from svmu2.orchestration.trend import run_trend
        run_trend(args)

    elif args.command == "inspect":
        from svmu2.orchestration.inspect import run_inspect
        run_inspect(args)

    elif args.command == "filter":
        from svmu2.orchestration.filter import run_filter
        run_filter(args)

    else:
        raise RuntimeError(f"Unknown command: {args.command}")

if __name__ == '__main__':
    main()