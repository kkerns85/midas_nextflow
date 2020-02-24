#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Analyze microbial pan-genomes using MIDAS

    Usage:

    nextflow run kkerns85/midas_nextflow <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output_folder       Folder to place analysis outputs (default ./midas)
      --output_prefix       Text used as a prefix for output files (default: midas)      

    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This value can not be repeated.
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime, or omits the --manifest
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}


// Set default options, which are overridden by the user with, e.g., --output_folder OTHER_VALUE
params.output_folder =  'midas'
params.output_prefix =  'midas'


// Parse the manifest CSV
// Along the way, make sure that the appropriate columns were provided
fastq_ch = Channel.from(
    file(
        params.manifest
    ).splitCsv(
        header: true,
        sep: ","
    )
).filter { 
    r -> (r.specimen != null)
}.ifEmpty { 
    exit 1, "Cannot find values in the 'specimen' column: ${params.manifest}"
}.filter { 
    r -> (r.R1 != null)
}.ifEmpty { 
    exit 1, "Cannot find values in the 'R1' column: ${params.manifest}"
}.filter { 
    r -> (r.R2 != null)
}.ifEmpty { 
    exit 1, "Cannot find values in the 'R2' column: ${params.manifest}"
}.map {
    r -> [r["specimen"], file(r["R1"]), file(r["R2"])]
}


params.input_type = "fastq"
params.input_type_knead = "*.kneaddata.trimmed.fastq"

process kneaddata {
    container "biobakery/kneaddata:0.7.3_cloud"
    label 'mem_medium'
    
    input:
    tuple val(specimen), file(R1), file(R2) from fastq_ch
    
    output:
    tuple val(specimen), file("*_kneaddata.trimmed.1.fastq"), file("*_kneaddata.trimmed.2.fastq") into trimmed_fastq_ch

    """
#!/bin/bash

set -e

kneaddata --input ${R1} --input ${R2} --output ./
"""
}

process midas {
    container "https://github.com/FredHutch/docker-midas.git"
    cpus 16
    memory "256 GB"
    publishDir "${params.output_folder}"
    input:
    file input_fastq from (***Should be from kneaddata output?***)
    val input_type from params.input_type_knead
    output:
    file "${input_fastq}.midas.species.tsv",
    file "${input_fastq}.midas.genes.tsv",
    file "${input_fastq}.midas.snp.tsv",
    file "${input_fastq}.midas.merged.tsv"
    """
    run.py --input_type ${input_type} --tmp_dir ./ -o ${input_fastq}.midas.tsv ${input_fastq}
    """
}



