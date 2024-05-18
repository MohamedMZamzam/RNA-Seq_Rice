# Set the working directory where your output files to be saved
setwd("D:/folder/subfolder")
# Load required libraries#

library(Rsubread)
library (tidyverse)
library(Rsamtools)
library (ShortRead)
library (Rfastp)
#run QC with Rfastp to remove the adapotor and he low quality reads

#We also need to provide the prefix we want for the outputs to outputFastq.

#It is recommended to give the adaptor sequnece as it is not always recognized by the rfastp

json_report <- rfastp(read1 = "R1.fastq.gz", read2 = "R2.fastq.gz", adapterSequenceRead2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", outputFastq = "unpaired" )
#to see the QCSummary
qcSummary(json_report)

#Load the transcript annotation file
Transcript <- rtracklayer::readGFF("D://IGV/RAP REF GENOME/transcripts(1).gff") %>% as.data.frame

#make the annotation file data frame
SAF <- data.frame(GeneID = Transcript$Locus_id, Chr = Transcript$seqid, Start = Transcript$start,
                  End = Transcript$end, Strand = Transcript$strand)

#Map the trimmed reads to reference genome#
MyMapped <- subjunc("mainchrs_rap",
                    readfile1 = "unpaired_R1.fastq.gz", 
                    readfile2 = "unpaired_R2.fastq.gz", 
                    output_format = "BAM",
                    output_file = "Mapped.rap.bam", 
                    useAnnotation = TRUE, 
                    annot.ext = SAF, 
                    isGTF = FALSE)

#Sort and Index the reads in the .bam file
sortBam("Mapped.rap.bam", "Sorted_Mapped.rap")
indexBam("Sorted_Mapped.rap.bam")

#This sorted file can be investigated with IGV browser,
# for featurecount to do the DGE and can be further investigated by sametools for looking at the unmapped reads which may contains your transgene and other important data