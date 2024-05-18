# Set the working directory where your output files to be saved
setwd("D:/Sources/Ref_genomes")

#Load libraries# 
library(Rsubread)

buildindex("mainchrs_rap", "IRGSP-1.0_genome.fasta", memory = 8000, indexSplit = TRUE)
