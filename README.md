# midas_nextflow
Running Midas using nextflow

****Test****

This workflow is designed to help streamline bacterial metagenomic and metatranscriptomic data analysis using the Nextflow workflow manager. Raw reads (fastq or fastq.gz) are initailly passed through Kneaddata which initially trims raw reads via Trimmomatic (v. 0.33) and then aligns them via Bowtie2 (v. >= 2.2) and finally removes bacterial ribosomal RNA using the Silva rRNA database (v. 128). Trimmed, aligned, and filtered reads are then passed on to MIDAS in order to assign taxonomy using the Phy-Eco single-copy marker gene set which provides species id with strain level variation, gene content, and single-nucleotide-polymorphisms (SNP's) for each metagenome.

This minimal test_data set was developed using Saccharibacteria nanosynbacter lyticus HMT 952, formerly TM7x, a egnimatic member of the candidate phyla radiation (CPR), further emphazing the application of this workflow. 

More information is available at: 

Kneaddata: http://huttenhower.sph.harvard.edu/kneaddata

Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

MiDAS: https://github.com/snayfach/MIDAS

Maybe use FA.nf: https://github.com/guigolab/FA-nf/blob/master/flowchart.png

Nextflow: https://www.nextflow.io

McLean Lab: Univeristy of Washington School of Dentistry, Department of Periodontics

# Authors
> Sam Minot and 
> Kristopher Kerns 


