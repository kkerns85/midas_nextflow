#!/usr/bin/env nextflow

fastq_ch = Channel.from(file(params.manifest).readLines())
                  .map {it -> file(it)}

params.input_type = "fastq"
params.input_type_knead = "kneaddata.trimmed.fastq"

process kneaddata {
    container "https://github.com/brianmorganpalmer/kneaddata.git"
    cpus 16
    memory "256 GB"
    publishDir "${params.output_folder}"
    input:
    file input_fastq from fastq_ch
    val input_type from params.input_type
    output:
    file "${input_fastq}.kneaddata.trimmed.fastq"
    reference-db:
    file "s3://mcleanlabmidas/Kneaddata_database/ribosomal_RNA_db/"
    """
    kneaddata --input $INPUT --reference-db $DATABASE --output $OUTPUT_DIR
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
    file "${input_fastq}.midas.tsv"
    """
    run.py --input_type ${input_type} --tmp_dir ./ -o ${input_fastq}.midas.tsv ${input_fastq}
    """
}
