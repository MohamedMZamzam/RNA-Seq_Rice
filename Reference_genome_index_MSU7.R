#set the working directory
setwd("D:/Sources/Ref_genomes/")
#Load the ibraries
library (BSgenome.Osativa.MSU.MSU7)
library(Rsubread)

#generate ref genome index from BSgenome.Osativa.MSU.MSU7
mainChromosomes <- paste0("Chr", c(1:12, "M", "C", "Un", "Sy"))
mainChrSeq <- lapply(mainChromosomes, function(x) BSgenome.Osativa.MSU.MSU7[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)

writeXStringSet(mainChrSeqSet, "BSgenome.Osativa.MSU.MSU7.mainChrs.fa")


buildindex("mainchrs_MSU7", "BSgenome.Osativa.MSU.MSU7.mainChrs.fa", memory = 8000,
           indexSplit = TRUE)
