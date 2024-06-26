# RNA-Seq for rice
<p align="justify">

  
## QC and read mapping in Rice
To perform QC on the rawdata fasq files and remove low qality reads and adaptor sequences etc. Then to mapped the raw reads to reference genome, sort and index the bam file (mapped and unmampped reads)

### Files and datasets
1. Raw data files (R1.fastq.gz and R2.fastq.gz): These are the file from your own experiments or files downloaded from publical repositories like NCBI
2. The transcript annotation file (transcripts.gff) [download from here](https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2024-01-11.tar.gz)
3. The transcript annotation file (all.gff3) [download from here](http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir). This file is used when you are mapping with respect to rice MSU locus ID
4.  Reference genome index (mainchrs_rap). This is generated from Os-Nipponbare-Reference-IRGSP-1.0 (IRGSP-1.0_genome.fasta) [download from here](https://rapdb.dna.affrc.go.jp/download/irgsp1.html)  as demonstrated in [Reference_genome_index_rap.R](https://github.com/MohamedMZamzam/RNA-Seq_Rice/blob/main/Reference_genome_index_rap.R) 

## Feature Counts and DEG analysis
Count the raw and normalized reads mapped to each genomic feature and compare the differential expression levels
### Files and datasets
1. we used the sorted bam files from the previous step
2. For making the SAF, we used the transcript annotation file (transcripts.gff) or the locus annotation file (locus.gff) [download from here](https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2024-01-11.tar.gz)

## Annotation
Annotate the list of DEG from the information available on database
### Files and datasets
1. Annotation files can be downloaded from different database. The example file used here is downloaded from Oryzabase databse

