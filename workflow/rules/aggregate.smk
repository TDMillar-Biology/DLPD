rule plot_contiguity:
    input:
        # Ensures the plot only runs after ALL curated assemblies are finished
        assemblies = expand("results/{strain}/curated_assembly/{strain}.curated.fasta", strain=STRAINS),
        # Uses the list generated at the top of your Snakefile
        refs = REF_FASTAS
    output:
        plot = "results/aggregate/plots/contiguity.png"
    resources:
        mem_mb=1600,
        runtime=60,
        ntasks=1
    log:
        "logs/contiguity_plot/contiguity.log"
    container:
        "workflow/containers/images/python_mummer.sif"
    params:
        results_dir = "results",
        # Pattern to find the curated assemblies within each strain folder
        pattern = "curated_assembly/*.curated.fasta",
        max_contigs = 20
    shell:
        """
        # Ensure the output directory exists
        mkdir -p $(dirname {output.plot})
        
        python3 workflow/scripts/plot_contiguity.py \
            --results {params.results_dir} \
            --reference-fastas {input.refs} \
            --pattern "{params.pattern}" \
            --output {output.plot} > {log} 2>&1
        """

rule aggregate_metrics:
    input:
        seqkit=expand("results/{strain}/qc/{strain}.seqkit.stats.tsv", strain=STRAINS),
        qv=expand("results/{strain}/qc/{strain}_qc.tsv", strain=STRAINS)
    output:
        "results/aggregate/tables/assembly_metrics.tsv"
    resources:
        mem_mb=1600,
        runtime=60,
        ntasks=1
    log:
        "logs/aggregate/aggregate.log"
    container:
        "workflow/containers/images/python_mummer.sif"
    shell:
        """
        mkdir -p logs/aggregate results/aggregate/tables

        python3 workflow/scripts/aggregate_metrics.py \
            --results results \
            --output {output} \
            > {log} 2>&1
        """
