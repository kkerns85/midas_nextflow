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
      --db_midas            Folder containing the MIDAS database
      --db_knead            Folder containing the Kneaddata database
    Options:
      --output_folder       Folder to place analysis outputs (default ./midas)
      --output_prefix       Text used as a prefix for output files (default: midas)
      --species_cov         Coverage (depth) threshold for species inclusion (default: 3.0)
      --single              Input data is single-end (default: treat as paired-end)
      --merge_sample_depth  Corresponds to the --sample_depth parameter in the merge_midas.py command (default: 1.0)
      --no_knead            Skip Kneaddata (default: false)
      
    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This value can not be repeated.
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If you specify --single, then only data from `R1` will be used
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime, or omits the --manifest
if (params.help || params.manifest == null || params.db_midas == null || params.db_knead == null){
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

// Make sure that the Midas database file can be found
if (file(params.db_midas).isEmpty()){

    // Print a helpful log message
    log.info"""
    Cannot find the file specified by --db_midas ${params.db_midas}
    """.stripIndent()

    // Exit out and do not run anything else
    exit 0
}

// Make sure that the Kneaddata database file can be found
if (file(params.db_midas).isEmpty()){

    // Print a helpful log message
    log.info"""
    Cannot find the file specified by --db_knead ${params.db_knead}
    """.stripIndent()

    // Exit out and do not run anything else
    exit 0
}

// Set default options, which are overridden by the user with, e.g., --output_folder OTHER_VALUE
params.output_folder =  'midas'
params.output_prefix =  'midas'
params.species_cov = 3.0
params.merge_sample_depth = 1.0
params.single = false
params.no_knead = false


// Parse the manifest CSV
// Along the way, make sure that the appropriate columns were provided
if (params.single || params.no_knead){
    trimmed_fastq_ch = Channel.from(
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
    }.map {
        r -> [r["specimen"], [file(r["R1"])]]
    }
} else if (params.single) {
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
    }.map {
        r -> [r["specimen"], [file(r["R1"])]]
    }
} else if (params.no_knead){
    trimmed_fastq_ch = Channel.from(
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
        r -> [r["specimen"], [file(r["R1"]), file(r["R2"])]]
    }
} else {
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
        r -> [r["specimen"], [file(r["R1"]), file(r["R2"])]]
    }
}

process kneaddata {
    container "biobakery/kneaddata:0.7.5_cloud"
    label 'mem_medium'

    input:
    tuple val(specimen), file("${specimen}.R*.fastq.gz") from fastq_ch

    output:
    tuple val(specimen), file("${specimen}.R1_kneaddata.trimmed.*.fastq.gz") into trimmed_fastq_ch

    """
#!/bin/bash
set -e
if [[ -s ${specimen}.R2.fastq.gz ]]; then
    kneaddata --input ${specimen}.R1.fastq.gz --input ${specimen}.R2.fastq.gz --output ./ -t ${task.cpus}
else
    mv ${specimen}.R.fastq.gz ${specimen}.R1.fastq.gz
    kneaddata --input ${specimen}.R1.fastq.gz --output ./ -t ${task.cpus}
    mv ${specimen}.R1_kneaddata.trimmed.fastq ${specimen}.R1_kneaddata.trimmed.1.fastq
fi
gzip ${specimen}.R1_kneaddata.trimmed.[12].fastq
"""
}

process midas {
    container "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    label "mem_veryhigh"
    publishDir "${params.output_folder}/${specimen}"

    input:
    tuple val(specimen), file("${specimen}.R*.fastq.gz") from trimmed_fastq_ch
    file DB from file(params.db_midas)

    output:
    file "${specimen}.species.tar.gz" into species_ch
    file "${specimen}.genes.tar.gz" into gene_ch
    file "${specimen}.snps.tar.gz" into snps_ch

"""
#!/bin/bash
set -e
echo "Running species summary"
# If the input is single-end, change the filename to match the pattern used for paired-end
if [[ ! -s ${specimen}.R2.fastq.gz ]]; then
    mv ${specimen}.R.fastq.gz ${specimen}.R1.fastq.gz
fi
# Run the same command differently, depending on whether the input is single- or paired-end
if [[ -s ${specimen}.R2.fastq.gz ]]; then
    # Run the species abundance summary
    run_midas.py \
        species \
        ${specimen} \
        -1 ${specimen}.R1.fastq.gz \
        -2 ${specimen}.R2.fastq.gz \
        -t ${task.cpus} \
        -d ${DB}
else
    # Run the species abundance summary
    run_midas.py \
        species \
        ${specimen} \
        -1 ${specimen}.R1.fastq.gz \
        -t ${task.cpus} \
        -d ${DB}
fi
# Run the gene abundance summary
if [[ -s ${specimen}.R2.fastq.gz ]]; then
    echo "Running gene summary"
    run_midas.py \
        genes \
        ${specimen} \
        -1 ${specimen}.R1.fastq.gz \
        -2 ${specimen}.R2.fastq.gz \
        -t ${task.cpus} \
        -d ${DB} \
        --species_cov ${params.species_cov}
else
    echo "Running gene summary"
    run_midas.py \
        genes \
        ${specimen} \
        -1 ${specimen}.R1.fastq.gz \
        -t ${task.cpus} \
        -d ${DB} \
        --species_cov ${params.species_cov}
fi
# Run the SNP summary
echo "Running SNP summary"
if [[ -s ${specimen}.R2.fastq.gz ]]; then
    run_midas.py \
        snps \
        ${specimen} \
        -1 ${specimen}.R1.fastq.gz \
        -2 ${specimen}.R2.fastq.gz \
        -t ${task.cpus} \
        -d ${DB} \
        --species_cov ${params.species_cov}
else
    run_midas.py \
        snps \
        ${specimen} \
        -1 ${specimen}.R1.fastq.gz \
        -t ${task.cpus} \
        -d ${DB} \
        --species_cov ${params.species_cov}
fi
echo "Gathering output files"
# Species-level results
echo "Tarring up species results"
tar cvf ${specimen}.species.tar ${specimen}/species/*
gzip ${specimen}.species.tar
# Gene-level results
echo "Tarring up gene results"
tar cvf ${specimen}.genes.tar ${specimen}/genes/*
gzip ${specimen}.genes.tar
# SNP-level results
echo "Tarring up SNP results"
tar cvf ${specimen}.snps.tar ${specimen}/snps/*
gzip ${specimen}.snps.tar
echo "Done"
"""
}

process midas_merge_species {
    container "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    label "mem_veryhigh"
    publishDir "${params.output_folder}"

    input:
    file species_tar_list from species_ch.toSortedList()
    file DB from file(params.db_midas)

    output:
    file "SPECIES/*"

"""
#!/bin/bash
set -e
ls -lahtr
# Keep track of the folders created while unpacking input files
input_string=""
echo "Unpacking all of the input files"
for tarfile in ${species_tar_list}; do
    echo "Making sure that \$tarfile was downloaded correctly"
    [[ -s \$tarfile ]]
    echo "Unpacking \$tarfile"
    tar xzvf \$tarfile
    # Add this folder to the input string
    input_string="\$input_string,\$( echo \$tarfile | sed 's/.species.tar.gz//' )"
    echo "Updated input string: \$input_string"
done
# Remove the leading comma from the input string
input_string=\$( echo \$input_string | sed 's/^,//' )
echo "Merging species results"
merge_midas.py \
    species \
    SPECIES \
    -i \$input_string \
    -t list \
    -d ${DB} \
    --sample_depth ${params.merge_sample_depth}
echo "Done merging data"
ls -lahtr SPECIES
echo "Compressing output files"
find SPECIES -type f | xargs gzip
echo "Done"
"""
}


process midas_merge_genes {
    container "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    label "mem_veryhigh"
    publishDir "${params.output_folder}"

    input:
    file genes_tar_list from gene_ch.toSortedList()
    file DB from file(params.db_midas)

    output:
    file "GENES/*"

"""
#!/bin/bash
set -e
ls -lahtr
# Keep track of the folders created while unpacking input files
input_string=""
echo "Unpacking all of the input files"
for tarfile in ${genes_tar_list}; do
    echo "Making sure that \$tarfile was downloaded correctly"
    [[ -s \$tarfile ]]
    echo "Unpacking \$tarfile"
    tar xzvf \$tarfile
    # Add this folder to the input string
    input_string="\$input_string,\$( echo \$tarfile | sed 's/.genes.tar.gz//' )"
    echo "Updated input string: \$input_string"
done
# Remove the leading comma from the input string
input_string=\$( echo \$input_string | sed 's/^,//' )
echo "Merging gene results"
merge_midas.py \
    genes \
    GENES \
    -i \$input_string \
    -t list \
    -d ${DB} \
    --sample_depth ${params.merge_sample_depth}
echo "Done merging data"
ls -lahtr GENES
echo "Compressing output files"
find GENES -type f | xargs gzip
echo "Done"
"""
}

process midas_merge_snps {
    container "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    label "mem_veryhigh"
    publishDir "${params.output_folder}"

    input:
    file snps_tar_list from snps_ch.toSortedList()
    file DB from file(params.db_midas)

    output:
    file "SNPS/*"

"""
#!/bin/bash
set -e
ls -lahtr
# Keep track of the folders created while unpacking input files
input_string=""
echo "Unpacking all of the input files"
for tarfile in ${snps_tar_list}; do
    echo "Making sure that \$tarfile was downloaded correctly"
    [[ -s \$tarfile ]]
    echo "Unpacking \$tarfile"
    tar xzvf \$tarfile
    # Add this folder to the input string
    input_string="\$input_string,\$( echo \$tarfile | sed 's/.snps.tar.gz//' )"
    echo "Updated input string: \$input_string"
done
# Remove the leading comma from the input string
input_string=\$( echo \$input_string | sed 's/^,//' )
echo "Merging snps results"
merge_midas.py \
    snps \
    SNPS \
    -i \$input_string \
    -t list \
    -d ${DB} \
    --sample_depth ${params.merge_sample_depth}
echo "Done merging data"
touch SNPS/DONE
ls -lahtr SNPS
echo "Compressing output files"
find SNPS -type f | xargs gzip
echo "Done"
"""
}
