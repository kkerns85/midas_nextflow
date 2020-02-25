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
      --db                  Folder containing the MIDAS database

    Options:
      --output_folder       Folder to place analysis outputs (default ./midas)
      --output_prefix       Text used as a prefix for output files (default: midas)
      --species_cov         Coverage (depth) threshold for species inclusion (default: 3.0)

    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This value can not be repeated.
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime, or omits the --manifest
if (params.help || params.manifest == null || params.db == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make sure that the manifest file can be found
if (file(params.manifest).isEmpty()){

    // Print a helpful log message
    log.info"""
    Cannot find the file specified by --manifest ${params.manifest}
    """.stripIndent()

    // Exit out and do not run anything else
    exit 0
}

// Make sure that the database file can be found
if (file(params.db).isEmpty()){

    // Print a helpful log message
    log.info"""
    Cannot find the file specified by --db ${params.db}
    """.stripIndent()

    // Exit out and do not run anything else
    exit 0
}

// Set default options, which are overridden by the user with, e.g., --output_folder OTHER_VALUE
params.output_folder =  'midas'
params.output_prefix =  'midas'
params.species_cov = 3.0


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
    container "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    label "mem_veryhigh"
    publishDir "${params.output_folder}"

    input:
    tuple val(specimen), file(R1), file(R2) from trimmed_fastq_ch
    file DB from file(params.db)

    output:
    file "${specimen}.midas.species.txt.gz"
    file "${specimen}.midas.species.log.gz"
    file "${specimen}.midas.genes.log.gz"
    file "${specimen}.midas.genes.summary.txt.gz"
    file "${specimen}.midas.genes.tar"
    file "${specimen}.midas.snps.log.gz"
    file "${specimen}.midas.snps.summary.txt.gz"
    file "${specimen}.midas.snps.tar"
"""
#!/bin/bash

set -e

echo "Running species summary"
    
# Run the species abundance summary
run_midas.py \
    species \
    OUTPUT \
    -1 ${R1} \
    -2 ${R2} \
    -t ${task.cpus} \
    -d ${DB}

# Run the gene abundance summary
echo "Running gene summary"
run_midas.py \
    genes \
    OUTPUT \
    -1 ${R1} \
    -2 ${R2} \
    -t ${task.cpus} \
    -d ${DB} \
    --species_cov ${params.species_cov}

# Run the SNP summary
echo "Running SNP summary"
run_midas.py \
    snps \
    OUTPUT \
    -1 ${R1} \
    -2 ${R2} \
    -t ${task.cpus} \
    -d ${DB} \
    --species_cov ${params.species_cov}

echo "Gathering output files"

# Species-level results
mv OUTPUT/species/species_profile.txt ${specimen}.midas.species.txt
mv OUTPUT/species/log.txt ${specimen}.midas.species.log

# Gene-level results
mv OUTPUT/genes/log.txt ${specimen}.midas.genes.log
mv OUTPUT/genes/summary.txt ${specimen}.midas.genes.summary.txt
tar cvf ${specimen}.midas.genes.tar OUTPUT/genes/output/*

# SNP-level results
mv OUTPUT/snps/log.txt ${specimen}.midas.snps.log
mv OUTPUT/snps/summary.txt ${specimen}.midas.snps.summary.txt
tar cvf ${specimen}.midas.snps.tar OUTPUT/snps/output/*

# Compress output files
gzip ${specimen}.midas.*log
gzip ${specimen}.midas.*txt

echo "Done"

"""
}



