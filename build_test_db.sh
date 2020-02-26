#!/bin/bash

set -e

NXF_VER=20.01.0 \
    nextflow \
    run \
    build_db.nf \
    -c nextflow.config \
    -profile testing \
    --genome_folder test_data/genomes \
    --mapfile test_data/genomes.mapfile \
    --output_folder db/ \
    -with-docker ubuntu:18.04 \
    -w work/ \
    -process.executor local \
    -with-report report.build_db.html \
    -with-trace trace.build_db.txt \
    -resume
