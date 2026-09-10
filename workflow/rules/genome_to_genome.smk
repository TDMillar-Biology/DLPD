rule align_iso1_to_assembly:
    input:
        target_assembly="results/{strain}/curated_assembly/{strain}.curated.fasta",
        query_iso1=config["references"]["ISO1"]["full"]
    output:
        paf="results/{strain}/mapping/ISO1_to_{strain}.raw.paf"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        ntasks=1,
        slurm_partition="medium"
    container:
        "workflow/containers/images/mapping_qc.sif"
    log:
        "logs/mapping/ISO1_to_{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mapping logs/mapping

        # Direct PAF emission with --cs for paftools.js compatibility
        minimap2 -cx asm5 --cs --secondary=no -t {threads} \
            {input.target_assembly} {input.query_iso1} \
            > {output.paf} 2> {log}
        """

rule svmu_filter_main_diagonal:
    input:
        paf="results/{strain}/mapping/ISO1_to_{strain}.raw.paf"
    output:
        filtered_paf=temp("results/{strain}/mapping/ISO1_to_{strain}.filtered.paf"),
        sorted_paf="results/{strain}/mapping/ISO1_to_{strain}.filtered.sorted.paf"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=30,
        ntasks=1,
        slurm_partition="short"
    container:
        "workflow/containers/images/svmu2.sif"
    log:
        "logs/mapping/svmu_filter_{strain}.log"
    shell:
        """
        # 1. Filter PAF to the main diagonal using native PAF format
        svmu2 filter \
            -a {input.paf} \
            --format paf \
            -o {output.filtered_paf} \
            > {log} 2>&1

        # 2. Sort directly by Target Name (k6) and Target Start (k8) for paftools.js
        sort -k6,6 -k8,8n {output.filtered_paf} > {output.sorted_paf} 2>> {log}
        """

rule paftools_call_micro_variants:
    """
    Calls small indels and SNPs from the purified, sorted syntenic PAF.
    """
    input:
        target_assembly="results/{strain}/curated_assembly/{strain}.curated.fasta",
        paf="results/{strain}/mapping/ISO1_to_{strain}.filtered.sorted.paf"
    output:
        vcf="results/{strain}/mappable_variants/ISO1_to_{strain}.paftools_micro.vcf"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=30,
        ntasks=1,
        slurm_partition="short"
    container:
        "workflow/containers/images/mapping_qc.sif"
    log:
        "logs/variants/paftools_call_{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/variants logs/variants
        paftools.js call -f {input.target_assembly} {input.paf} > {output.vcf} 2> {log}
        """