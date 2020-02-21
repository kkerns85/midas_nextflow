
Description and statistics for pan-genome files

Summary Statistics
############

Genomes: 1
Genes: 757
Gene clusters (99% identity): 753
Gene clusters (95% identity): 753
Gene clusters (90% identity): 752
Gene clusters (85% identity): 751
Gene clusters (80% identity): 748
Gene clusters (75% identity): 745
		
Output files
############
genes.ffn
  all genes from specified genomes
  
centroids.ffn
  gene sequences from 99% identity gene clusters
  used for recruiting metagenomic reads
  
gene_info.txt
  information for all genes from genes.ffn
  the fields centroid_{99,95,90,95,80,75} indicate mappings between gene_id and gene clusters
