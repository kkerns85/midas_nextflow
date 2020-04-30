# midas_nextflow
Running Midas using Nextflow
Optimized for Running on AWS Batch

****Beta****

This workflow is designed to help streamline bacterial metagenomic and metatranscriptomic data analysis using the Nextflow workflow manager. Raw reads (fastq or fastq.gz) are initailly passed through Kneaddata which initially trims raw reads via Trimmomatic (v. 0.33) and then aligns them via Bowtie2 (v. >= 2.2) and finally removes bacterial ribosomal RNA using the Silva rRNA database (v. 128). Trimmed, aligned, and filtered reads are then passed on to MIDAS in order to assign taxonomy using the Phy-Eco single-copy marker gene set and pangenome database (v1.2) which provides species taxonomic level assignment, gene content, and strain level variation via single-nucleotide-polymorphism (SNP's) analysis for each metagenome.

This minimal test_data set was developed using Saccharibacteria Nanosynbacter lyticus HMT 952, formerly TM7x, a egnimatic member of the candidate phyla radiation (CPR), further emphazing the application of this workflow. 

#Install Nextflow
- https://www.nextflow.io/docs/latest/getstarted.html

#Input 
- Manifest file (Local) 
  - Path to paried or single end fastq or fastq.gz (R1 and R2) (Local or S3://)
  - Metadata additional Columns with simple header format 
- Work directory (Local or S3://)
- Path to Midas Database (Local or S3://)
- Path to Kneaddata Database (Human, Mouse, and/or Bacterial rRNA) (Local or S3://)

#Output
- Trimmed and Filtered Reads from Kneaddata (Stored Local or S3://)
- Species Analysis from MiDAS
- Gene Analysis from MiDAS
- SNP Analysis from MiDAS
- Merged files for Species, Genes, and SNP analysis

#Options
    --single            Run single end reads
    --no-knead          Skip kneaddata process
    --profile           local or AWS batch
    --output_folder     Folder to place analysis outputs (default ./midas)
    --output_prefix     Text used as a prefix for output files (default: midas)
    --species_cov       Coverage (depth) threshold for species inclusion (default: 3.0)
    --merge_sample_depth  Corresponds to the --sample_depth parameter in the merge_midas.py command (default: 1.0)

#Manifest
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This value can not be repeated.
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If you specify --single, then only data from `R1` will be used

Initial Analysis using R (Coming Soon!)


More information is available at: 

Kneaddata: http://huttenhower.sph.harvard.edu/kneaddata

Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

MiDAS: https://github.com/snayfach/MIDAS 

MiDAS Docker: https://github.com/FredHutch/docker-midas

Nextflow: https://www.nextflow.io

McLean Lab: Univeristy of Washington School of Dentistry, Department of Periodontics

# Authors
> Sam Minot and 
> Kristopher Kerns 


