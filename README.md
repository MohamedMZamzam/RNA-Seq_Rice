# RNA-Seq for rice
<p align="justify">

  
## QC and read mapping in Rice
To perform QC on the rawdata fasq files and remove low qality reads and adaptor sequences etc. Then to mapped the raw reads to reference genome, sort and index the bam file (mapped and unmampped reads)

### Files and datasets
1. Raw data files (R1.fastq.gz and R2.fastq.gz): These are the file from your own experiments or files downloaded from publical repositories like NCBI
2. The transcript annotation file (all.gff3) [download from here](http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir). This file is used when you are mapping with respect to rice MSU locus ID
3.  Reference genome index (mainchrs_rap). This is generated from Os-Nipponbare-Reference-IRGSP-1.0 (IRGSP-1.0_genome.fasta) as demonstrated in [Reference_genome_index_rap.R](https://github.com/MohamedMZamzam/RNA-Seq_Rice/blob/main/Reference_genome_index_rap.R) 



