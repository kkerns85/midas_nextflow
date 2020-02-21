#!/bin/bash

set -e

nextflow \
    run \
    /Users/kriskerns_home/midas_workflow.nf \
    --manifest /Users/kriskerns_home/manifest.txt\
    --input_folder s3://mcleanlab-nextflow-midas/scratch-delete-30/nextflow/input/ \
    --output_folder s3://mcleanlab-nextflow-midas/scratch-delete-30/nextflow/output/ \
    -work-dir s3://mcleanlab-nextflow-midas/scratch-delete-30/nextflow/work/ \
    -with-docker "ubuntu:16.04" \
