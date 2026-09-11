rule assemble:
    input:
        reads=get_reads
    output:
        hap1="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.gfa",
        hap2="results/{strain}/assembly/{strain}.bp.hap2.p_ctg.gfa"
    threads: 16
    container:
        "workflow/containers/images/assembly.sif"
    resources:
        mem_mb=64000,
        runtime=600,
        tasks=1
    shell:
        """
        mkdir -p results/{wildcards.strain}/assembly

        hifiasm -t {threads} \
            -o results/{wildcards.strain}/assembly/{wildcards.strain} \
            {input.reads} \
            > results/{wildcards.strain}/assembly/hifiasm.log 2>&1
        """

rule gfa2fa:
    input:
        hap1="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.gfa",
        hap2="results/{strain}/assembly/{strain}.bp.hap2.p_ctg.gfa",
        primary="results/{strain}/assembly/{strain}.bp.p_ctg.gfa"
    output:
        hap1fasta="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.fasta",
        hap2fasta="results/{strain}/assembly/{strain}.bp.hap2.p_ctg.fasta",
        primaryfasta="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    threads: 1
    resources:
        mem_mb=7500,
        runtime=10,
        tasks=1
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' {input.hap1} > {output.hap1fasta}
        awk '/^S/{{print ">"$2"\\n"$3}}' {input.hap2} > {output.hap2fasta}
        awk '/^S/{{print ">"$2"\\n"$3}}' {input.primary} > {output.primaryfasta}
        """

rule mummer_for_curation:
    input:
        reference=config["references"]["ISO1"]["main"],
        query="results/{strain}/assembly/{strain}.bp.p_ctg.fasta",
    output:
        delta="results/{strain}/mummer/{strain}_r6_main.delta"
    params:
        prefix="results/{strain}/mummer/{strain}_r6_main"
    resources:
        mem_mb=16000,
        runtime=120,
        tasks=1
    threads: 8
    container:
        "workflow/containers/images/python_mummer.sif"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mummer

        nucmer -p {params.prefix} \
            {input.reference} {input.query} -t {threads}
        """

rule dotplot:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main.delta"
    output:
        outdir=directory("results/{strain}/dotplots/for_curation")
    threads: 8
    resources:
        mem_mb=1600,
        runtime=60,
        ntasks=1
    container:
        "workflow/containers/images/svmu2.sif"
    log:
        "logs/dotplot/{strain}.log"
    shell:
        """
        mkdir -p {output.outdir}

        svmu2 plot \
            --alignment {input.delta} \
            --format delta \
            --out_dir {output.outdir} \
            --strain {wildcards.strain}_for_curation
        """

rule curate_assembly:
    input:
        fasta="results/{strain}/assembly/{strain}.bp.p_ctg.fasta",
        breaks="config/curation/{strain}.tsv"
    output:
        fasta="results/{strain}/curated_assembly/{strain}.curated.fasta",
        break_log="results/{strain}/curated_assembly/{strain}.curated.log"
    container:
        "workflow/containers/images/python_mummer.sif"
    resources:
        mem_mb=16000,
        runtime=60,
        tasks=1
    log: 
        "logs/asm_curation_{strain}.errorlog"
    shell:
        """
        mkdir -p results/{wildcards.strain}/curated_assembly/
    
        python3 workflow/scripts/breakasm.py \
            {input.fasta} \
            --breakfile {input.breaks} \
            --out-fasta {output.fasta} \
            --out-log {output.break_log} \
            > {log} 2>&1
        """

rule mummer_for_scaffolding:
    input:
        reference=config["references"]["ISO1"]["main"],
        query="results/{strain}/curated_assembly/{strain}.curated.fasta",
    output:
        delta="results/{strain}/mummer/{strain}_r6_main_curated.delta"
    params:
        prefix="results/{strain}/mummer/{strain}_r6_main_curated"
    resources:
        mem_mb=16000,
        runtime=120,
        tasks=1
    threads: 8
    container:
        "workflow/containers/images/python_mummer.sif"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mummer

        nucmer -p {params.prefix} \
            {input.reference} {input.query} -t {threads}
        """

rule scaffold_with_daedalus:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main_curated.delta",
        fasta="results/{strain}/curated_assembly/{strain}.curated.fasta"
    output:
        coords="results/{strain}/scaffold/{strain}.coords",
        fasta="results/{strain}/scaffold/{strain}.scaffolded.fasta",
        agp="results/{strain}/scaffold/{strain}.scaffolded.agp"
    params:
        prefix="results/{strain}/scaffold/{strain}"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
        ntasks=1
    container:
        "workflow/containers/images/daedalus.sif"
    log:
        "logs/scaffold/{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/scaffold

        # Step 1: delta → coords
        show-coords -rclT {input.delta} > {output.coords}
        
        # Step 2: run daedalus
        daedalus scaffold \
            --coords {output.coords} \
            --fasta {input.fasta} \
            --prefix {wildcards.strain}.scaffolded \
            --out_dir results/{wildcards.strain}/scaffold \
            --chromosomes X 2L 2R 3L 3R 4 \
            > {log} 2>&1
      """

rule dotplot_curated_all:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main_curated.delta"
    output:
        outdir=directory("results/{strain}/dotplots/curated_all")
    resources:
        mem_mb=1600,
        runtime=60,
        tasks=1
    container:
        "workflow/containers/images/svmu2.sif"
    log:
        "logs/dotplot/{strain}_curated.log"
    shell:
        """
        mkdir -p {output.outdir}

        svmu2 plot \
            --alignment {input.delta} \
            --format delta \
            --out_dir {output.outdir} \
            --strain {wildcards.strain}_curated \
            --img-format png \
            --all \
            > {log} 2>&1
        """

rule dotplot_curated:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main_curated.delta"
    output:
        outdir=directory("results/{strain}/dotplots/curated/")
    resources:
        mem_mb=1600,
        runtime=60,
        tasks=1
    container:
        "workflow/containers/images/svmu2.sif"
    log:
        "logs/dotplot/{strain}_curated.log"
    shell:
        """
        mkdir -p {output.outdir}

        svmu2 plot \
            --alignment {input.delta} \
            --format delta \
            --out_dir {output.outdir} \
            --strain {wildcards.strain}_curated \
            --img-format png \
            > {log} 2>&1
        """

rule mummer_scaffolded:
    input:
        reference=config["references"]["ISO1"]["main"],
        query="results/{strain}/scaffold/{strain}.scaffolded.fasta"
    output:
        delta="results/{strain}/mummer/{strain}_r6_main_scaffolded.delta"
    params:
        prefix="results/{strain}/mummer/{strain}_r6_main_scaffolded"
    resources:
        mem_mb=16000,
        runtime=120,
        tasks=1
    threads: 8
    container:
        "workflow/containers/images/python_mummer.sif"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mummer

        nucmer -p {params.prefix} \
            {input.reference} {input.query} -t {threads}
        """


rule dotplot_scaffolded:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main_scaffolded.delta"
    output:
        outdir=directory("results/{strain}/dotplots/scaffolded")
    resources:
        mem_mb=1600,
        runtime=60,
        tasks=1
    container:
        "workflow/containers/images/svmu2.sif"
    log:
        "logs/dotplot/{strain}_scaffolded.log"
    shell:
        """
        mkdir -p {output.outdir}

        svmu2 plot \
            --alignment {input.delta} \
            --format delta \
            --out_dir {output.outdir} \
            --strain {wildcards.strain}_scaffolded \
            > {log} 2>&1
        """

rule plot_contig_flow:
    input:
        original="results/{strain}/assembly/{strain}.bp.p_ctg.fasta",
        curated="results/{strain}/curated_assembly/{strain}.curated.fasta",
        agp="results/{strain}/scaffold/{strain}.scaffolded.agp"
    output:
        plot="results/{strain}/qc/{strain}.contig_flow.png",
        table="results/{strain}/qc/{strain}.contig_flow.tsv"
    resources:
        mem_mb=1600,
        runtime=60,
        ntasks=1
    log:
        "logs/contig_flow/{strain}.log"
    container:
        "workflow/containers/images/python_mummer.sif"
    shell:
        """
        mkdir -p results/{wildcards.strain}/qc logs/contig_flow

        python3 workflow/scripts/plot_contig_flow.py \
            --original {input.original} \
            --curated {input.curated} \
            --agp {input.agp} \
            --out-plot {output.plot} \
            --out-table {output.table} \
            --title "{wildcards.strain} contig flow" \
            > {log} 2>&1
        """