snakemake --dag > dag.tmp
cat dag.tmp | dot -Tpng > figures/dag.png
