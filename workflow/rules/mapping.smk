rule map_reads:
    input:
        reads=get_reads,
        assembly="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    output:
        bam="results/{strain}/mapping/{strain}.bam",
        bai="results/{strain}/mapping/{strain}.bam.bai"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600,
        tasks=1
    container:
        "workflow/containers/images/mapping_qc.sif"
    log:
        "logs/mapping/{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mapping

        minimap2 -ax map-pb -t {threads} {input.assembly} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam} -O bam

        samtools index {output.bam} {output.bai} \
        > {log} 2>&1
        """

rule call_variants:
    input:
        bam="results/{strain}/mapping/{strain}.bam",
        assembly="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    output:
        outdir=directory("results/{strain}/variants")
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600,
        ntasks=1
    log:
        "logs/variants/{strain}.log"
    singularity:
        "containers/pmdv.sif"
    shell:
        """
        mkdir -p {output.outdir}

        run_pepper_margin_deepvariant call_variant \
            -b {input.bam} \
            -f {input.assembly} \
            -p "{wildcards.strain}_pmdv" \
            -o {output.outdir} \
            -t {threads} \
            --hifi \
            > {log} 2>&1
        """
