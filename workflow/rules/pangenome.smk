# rules/pangenome.smk

PANGENOME_REF = config['pangenome_reference']
PANGENOME_REF_PATH = config['pangenome_reference_path']

rule generate_cactus_seqfile:
    """
    Cactus requires a seqfile mapping strains to their assemblies.
    """
    input:
        fastas = expand("results/{strain}/scaffold/{strain}.scaffolded.fasta", strain=STRAINS)
    output:
        seqfile = "results/pangenome/seqfile.txt"
    resources:
        mem_mb=8000,
        runtime=30,
        ntasks=1
    run:
        with open(output.seqfile, "w") as out:
            out.write(f"{PANGENOME_REF}\t{PANGENOME_REF_PATH}\n")
            for strain, fasta in zip(STRAINS, input.fastas):
                out.write(f"{strain}\t{fasta}\n")

rule build_pangenome_graph:
    input:
        seqfile = "results/pangenome/seqfile.txt"
    output:
        outdir = directory("results/pangenome/cactus_out"),
        gfa = "results/pangenome/cactus_out/PANGENOME_REF.gfa.gz"
    log:
        "logs/pangenome/cactus_pangenome.log"
    threads: 32
    resources:
        mem_mb = 128000,
        runtime = 1440, # 1 day, matching the Grace medium partition limit
        ntasks = 1,
        slurm_partition = "long" # Overrides the default 'medium'
    container:
        "containers/cactus.sif"
    shell:
        """
        rm -rf ./jobstore_pangenome

        cactus-pangenome \
            ./jobstore_pangenome \
            {input.seqfile} \
            --outDir {output.outdir} \
            --outName PANGENOME_REF \
            --reference {PANGENOME_REF} \
            --maxCores {threads} \
            > {log} 2>&1
        """