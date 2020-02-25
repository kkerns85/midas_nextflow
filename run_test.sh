#!/bin/bash

set -e

echo "Processing paired-end data"
NXF_VER=20.01.0 \
    nextflow \
    run \
    -c nextflow.config \
    -profile testing \
    midas_workflow.nf \
    --manifest test_data/manifest.paired.csv \
    --db db/ \
    --output_folder test_output/paired/ \
    --species_cov 0.1 \
    -with-docker ubuntu:18.04 \
    -w work/ \
    -process.executor local \
    -with-report \
    -with-trace \
    -resume


echo "Processing single-end data"
NXF_VER=20.01.0 \
    nextflow \
    run \
    -c nextflow.config \
    -profile testing \
    midas_workflow.nf \
    --manifest test_data/manifest.single.csv \
    --db db/ \
    --output_folder test_output/single/ \
    --species_cov 0.1 \
    --single \
    -with-docker ubuntu:18.04 \
    -w work/ \
    -process.executor local \
    -with-report \
    -with-trace \
    -resume
