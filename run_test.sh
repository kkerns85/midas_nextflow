#!/bin/bash

set -e

NXF_VER=20.01.0 \
    nextflow \
    run \
    -c nextflow.config \
    -profile testing \
    midas_workflow.nf \
    --manifest test_data/manifest.txt \
    -with-docker ubuntu:18.04 \
    -w work/ \
    -process.executor local \
    -with-report \
    -with-trace \
    -resume
