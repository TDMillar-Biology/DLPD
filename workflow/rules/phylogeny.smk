# rules/phylogeny.smk

rule map_reads_for_phylogeny:
    input:
        reads=get_reads,
        reference=config["references"]["ISO1"]["full"]
    output:
        bam="results/{strain}/phylogeny/mapping/{strain}_ref.bam",
        bai="results/{strain}/phylogeny/mapping/{strain}_ref.bam.bai"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600,
        tasks=1
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/phylogeny/mapping/{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/phylogeny/mapping logs/phylogeny/mapping

        minimap2 -ax map-hifi -t {threads} {input.reference} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam} -O bam

        samtools index {output.bam} {output.bai} \
        > {log} 2>&1
        """

rule call_variants_for_phylogeny:
    input:
        bam="results/{strain}/phylogeny/mapping/{strain}_ref.bam",
        reference=config["references"]["ISO1"]["full"]
    output:
        vcf="results/{strain}/phylogeny/variants/{strain}_pmdv.vcf.gz",
        tbi="results/{strain}/phylogeny/variants/{strain}_pmdv.vcf.gz.tbi"
    params:
        outdir="results/{strain}/phylogeny/variants"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600,
        ntasks=1
    log:
        "logs/phylogeny/variants/{strain}.log"
    singularity:
        "containers/pmdv.sif"
    shell:
        """
        mkdir -p {params.outdir} logs/phylogeny/variants

        run_pepper_margin_deepvariant call_variant \
            -b {input.bam} \
            -f {input.reference} \
            -p "{wildcards.strain}_pmdv" \
            -o {params.outdir} \
            -t {threads} \
            --hifi \
            > {log} 2>&1
        """

rule merge_vcfs:
    input:
        vcfs = expand("results/{sample}/phylogeny/variants/{sample}_pmdv.renamed.vcf.gz", sample=PHYLOGENY_STRAINS),
        tbis = expand("results/{sample}/phylogeny/variants/{sample}_pmdv.renamed.vcf.gz.tbi", sample=PHYLOGENY_STRAINS)
    output:
        filtered_vcf = "results/aggregate/phylogeny/core_genome.filtered.vcf.gz"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
        ntasks=1
    conda:
        "../envs/phylogeny.yaml"
    log:
        "logs/phylogeny/merge_vcfs.log"
    shell:
        """
        mkdir -p results/aggregate/phylogeny logs/phylogeny

        # 1. Merge and automatically backfill missing genotypes to 0/0
        # 2. Force bcftools to recalculate accurate allele counts across all 12 strains
        # 3. Filter using the freshly updated, mathematically accurate MAC and AC values
        bcftools merge -0 --threads {threads} {input.vcfs} | \
        bcftools +fill-tags -O u | \
        bcftools view -v snps -m2 -M2 -i 'QUAL>30 && F_MISSING<0.1 && MAC>0 && AC!=AN' -O z -o {output.filtered_vcf} \
        > {log} 2>&1
        """

rule vcf_to_fasta:
    input:
        vcf="results/aggregate/phylogeny/core_genome.filtered.vcf.gz"
    output:
        fasta="results/aggregate/phylogeny/core_genome.filtered.min4.fasta"
    params:
        script="workflow/scripts/vcf2phylip.py",
        outdir="results/aggregate/phylogeny" # Moved outdir to params
    resources:
        mem_mb=8000,
        runtime=30,
        ntasks=1
    conda:
        "../envs/python.yaml"
    log:
        "logs/phylogeny/vcf_to_fasta.log"
    shell:
        """
        # Ensure the output directory exists before the script runs
        mkdir -p {params.outdir}
        
        # Pass the params.outdir to the script
        python3 {params.script} -i {input.vcf} --output-folder {params.outdir} -f > {log} 2>&1
        """

rule build_ml_tree:
    input:
        fasta="results/aggregate/phylogeny/core_genome.filtered.varsites.fasta"
    output:
        tree="results/aggregate/phylogeny/core_genome.filtered.varsites.fasta.treefile"
    threads: 16
    resources:
        mem_mb=32000,
        runtime=600,
        ntasks=1
    conda:
        "../envs/phylogeny.yaml"
    log:
        "logs/phylogeny/iqtree.log"
    shell:
        """
        iqtree2 -s {input.fasta} \
            -pre results/aggregate/phylogeny/core_genome.filtered.varsites.fasta \
            -m GTR+ASC -B 1000 -T {threads} \
            > {log} 2>&1
        """

rule reheader_vcf:
    input:
        vcf = "results/{sample}/phylogeny/variants/{sample}_pmdv.vcf.gz"
    output:
        vcf = "results/{sample}/phylogeny/variants/{sample}_pmdv.renamed.vcf.gz",
        tbi = "results/{sample}/phylogeny/variants/{sample}_pmdv.renamed.vcf.gz.tbi"
    log:
        "logs/phylogeny/reheader_{sample}.log"
    conda: 
        "../envs/phylogeny.yaml" # Make sure your bcftools/tabix env is referenced here
    threads: 1
    resources:
        mem_mb=32000,
        runtime=600,
        ntasks=1
    shell:
        """
        # Create a temporary file with the correct sample name
        echo "{wildcards.sample}" > {wildcards.sample}_name.txt
        
        # Apply the new sample name to the VCF header
        bcftools reheader -s {wildcards.sample}_name.txt {input.vcf} -o {output.vcf} 2> {log}
        
        # Index the new VCF (required for bcftools merge)
        tabix -p vcf {output.vcf} 2>> {log}
        
        # Clean up the temporary text file
        rm {wildcards.sample}_name.txt
        """

rule extract_variable_sites:
    input:
        fasta="results/aggregate/phylogeny/core_genome.filtered.min4.fasta"
    output:
        varsites="results/aggregate/phylogeny/core_genome.filtered.varsites.fasta"
    conda:
        "../envs/phylogeny.yaml"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=600,
        ntasks=1
    log:
        "logs/phylogeny/extract_variable_sites.log"
    shell:
        """
        # The -c flag keeps only variable sites, -o specifies the output file
        snp-sites -c -o {output.varsites} {input.fasta} > {log} 2>&1
        """