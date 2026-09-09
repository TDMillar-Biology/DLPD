#!/usr/bin/env bash
set -euo pipefail

sudo apptainer build \
    workflow/containers/images/python_mummer.sif \
    workflow/containers/defs/python_mummer.def

sudo apptainer build \
    workflow/containers/images/mapping_qc.sif \
    workflow/containers/defs/mapping_qc.def

sudo apptainer build \
    workflow/containers/images/busco.sif \
    workflow/containers/defs/busco.def

apptainer pull \
    workflow/containers/images/compleasm.sif \
    docker://quay.io/biocontainers/compleasm:0.2.9--pyhdfd78af_0

sudo apptainer build \
    workflow/containers/images/phylogeny.sif \
    workflow/containers/defs/phylogeny.def

sudo apptainer build \
    workflow/containers/images/assembly.sif \
    workflow/containers/defs/assembly.def

sudo apptainer build \
    workflow/containers/images/seqkit.sif \
    workflow/containers/defs/seqkit.def

sudo apptainer build \
    workflow/containers/images/daedalus.sif \
    workflow/containers/defs/daedalus.def

sudo apptainer build \
    workflow/containers/images/svmu2.sif \
    workflow/containers/defs/svmu2.def

