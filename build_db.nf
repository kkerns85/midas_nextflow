#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Build a microbial pan-genome database using MIDAS

    Usage:

    nextflow run kkerns85/midas_nextflow/build_db.nf <ARGUMENTS>
    
    Required Arguments:
      --genome_folder       Folder with input genomes
      --mapfile             Mapping file used by MIDAS to organize genomes
      --output_folder       Destination for database

    See the MIDAS documentation for the formatting requirements for MIDAS.

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime, or omits the --manifest
if (params.help || params.genome_folder == null || params.mapfile == null || params.output_folder == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

process buildMIDASdb {
    container "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    label 'mem_medium'
    publishDir "${params.output_folder}"
    
    input:
    path "GENOMES" from file(params.genome_folder)
    path "mapfile" from file(params.mapfile)
    
    output:
    path "*"

    """
#!/bin/bash

set -e

build_midas_db.py \
    GENOMES \
    mapfile \
    ./ \
    --threads ${task.cpus} \
    --compress
"""
}
