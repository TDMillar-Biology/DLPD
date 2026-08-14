rule align_iso1_to_assembly:
    input:
        target_assembly="results/{strain}/curated_assembly/{strain}.curated.fasta",
        query_iso1=config["references"]["ISO1"]["full"]
    output:
        bam="results/{strain}/mapping/ISO1_to_{strain}.bam",
        bai="results/{strain}/mapping/ISO1_to_{strain}.bam.bai"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        ntasks=1,
        slurm_partition="medium" 
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/mapping/ISO1_to_{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mapping logs/mapping

        # 1. minimap2: --eqx for extended CIGAR, --secondary=no for primary only
        # 2. samtools view: -F 2304 drops any lingering secondary/supplementary alignments
        # 3. samtools sort: organizes and compresses into the final BAM
        minimap2 -ax asm5 --eqx --secondary=no -t {threads} {input.target_assembly} {input.query_iso1} | \
        samtools view -F 2304 -u | \
        samtools sort -@ {threads} -o {output.bam} -O bam \
        > {log} 2>&1

        samtools index {output.bam} {output.bai} \
        >> {log} 2>&1
        """