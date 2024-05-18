#Set working directory
setwd(":/Rawdata (1)/double/")
#Load ibraries
library(Rsubread)
library (tidyverse)
library(Rsamtools)
library (ShortRead)
library (ggplot2)
library (Rfastp)

#prepare the feature annotation file
Transcript <- rtracklayer::readGFF("locus.gff") %>% as.data.frame

SAF <- data.frame(GeneID = Transcript$ID, Chr = Transcript$seqid,
                  Start = Transcript$start,
                  End = Transcript$end, Strand = Transcript$strand)

#Count feature for each of your samples

fc_WT_B1 <- featureCounts("Sorted_WT_B1.rap.bam", allowMultiOverlap = TRUE, 
                          countMultiMappingReads = TRUE,
                          annot.ext = SAF, isPairedEnd = TRUE)

fc_WT_B2 <- featureCounts("Sorted_WT_B2.rap.bam", allowMultiOverlap = TRUE,
                          countMultiMappingReads = TRUE,
                          annot.ext = SAF, isPairedEnd = TRUE)

fc_double_B1 <- featureCounts("Sorted_mutant_B1.rap.bam", allowMultiOverlap = TRUE,
                              countMultiMappingReads = TRUE,
                              annot.ext = SAF, isPairedEnd = TRUE)

fc_double_B2 <- featureCounts("Sorted_mutant_B2.rap.bam", allowMultiOverlap = TRUE,
                              countMultiMappingReads = TRUE,
                              annot.ext = SAF, isPairedEnd = TRUE)

#prepare unnormalized raw counts#

X <- data.frame(fc_WT_B1$annotation["GeneID"], 
                fc_WT_B1$counts, fc_WT_B2$counts, fc_double_B1$counts, 
                fc_double_B2$counts, stringsAsFactors = FALSE)

X <-`colnames<-`(X, c ("GeneID", "WT_B1", "WT_B2", "mutant_B1", "mutant_B2"))

write.table(X, file = "counts.xls", quote = FALSE, sep = "\t", row.names = FALSE)


#data grouping for comparison#
Counts <- X [, 2:5]
datagroup <- c ("WT", "WT", "mutant", "mutant")
d <- DGEList(counts = Counts, group = factor(datagroup))

keep <- filterByExpr(d, group = factor(datagroup))

d <- d[keep, , keep.lib.sizes=FALSE]


#The TMM normalization is applied to account for the compositional biases#

d <- calcNormFactors(d)

#An MDS plot shows the relative similarities of the six samples#

plotMDS(d, col=rep(1:2, each=2))

#PCA polt#

plotMDS(d, gene.selection = "common")

d1 <- estimateCommonDisp(d, verbose = T)
d1 <- estimateTagwiseDisp(d1)
dim (d1)

##pseudcounts normalized counts##
pseducounts <- d1[["pseudo.counts"]]
pseducounts_table <- cbind(rownames(pseducounts), pseducounts)
write.table (pseducounts_table, file = "pseducounts_.table.xls", quote = FALSE,
             sep = "\t", row.names = FALSE)
#Extract DEG#
de1 <- exactTest(d1, pair = c (2, 1))
topTags(de1, n= 23561)

n = topTags(de1, n= 23561)

l <- n[["table"]]

Table_l <- cbind(rownames(l), l)

#extract DEG that are significant as per specified threasholds#

sig <- filter (l, abs (logFC) >=1, FDR <= 0.05)

Table_sig <- cbind(rownames(sig), sig)

write.table(Table_sig, "Table_sig.xls", quote = FALSE,
            sep = "\t", row.names = FALSE)

