# midas_nextflow
Running Midas using nextflow

****Test****

This workflow is designed to help streamline bacterial metagenomic and metatranscriptomic data analysis using the Nextflow workflow manager. Raw reads (fastq or fastq.gz) are initailly passed through Kneaddata (Huttenhower lab_ verison) which initially trims raw reads via trimmomatic (lab, verison) and then something via Bowtie2 (lab verision) and finally removes bacterial ribosomal RNA using the Silva 128 rRNA database. Trimmed and filtered reads are then passed on to MIDAS (version _ lab) in order to assign taxonomy using the phy-eco single marker gene set which provides species id with strain level variation, gene content, and single-nucleotide-polymorphisms (SNP's). Assigned reads are merged with genes and then annotated using Prokka: rapid prokaryotic genome annotation. 

More information is available at: 

Kneaddata: http://huttenhower.sph.harvard.edu/kneaddata

Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

MiDAS: https://github.com/snayfach/MIDAS

Prokka: https://github.com/tseemann/prokka

Nextflow: https://www.nextflow.io

