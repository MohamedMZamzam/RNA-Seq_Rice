library (Rsamtools)
library (tidyverse)
library (Biostrings)
library(seqinr)
library (Rsubread)
library(VariantAnnotation)
#for dealing with Bam and Sam files#

SNP_WT__E_B1 <- Rsubread::exactSNP (readFile = "C://Rawdata (1)/OsM2/0.2-5/Sorted_WT02-05B1.rap.bam", isBAM = T,
                                 refGenomeFile = "D://IGV/RAP REF GENOME/IRGSP-1.0_genome.fasta", outputFile ="SNP_WT_B1.VCF",
                                 minAllelicFraction = 0.4)

VCF_WT_E_B1 <- VariantAnnotation::readVcf("SNP_WT_B1.VCF")

A <- as.data.frame(VCF_WT_E_B1@rowRanges)
Names_A <- rownames(A)
A <- cbind(A, Names_A)
View (A)
SNP_WT_E_B2 <- Rsubread::exactSNP (readFile = "C://Rawdata (1)/OsM2/0.2-5/Sorted_WT02-05B2.rap.bam", isBAM = T,
                                 refGenomeFile = "D://IGV/RAP REF GENOME/IRGSP-1.0_genome.fasta", outputFile ="SNP_WT_B2.VCF",
                                 minAllelicFraction = 0.4)

VCF_WT_E_B2 <- VariantAnnotation::readVcf("SNP_WT_B2.VCF")
B<-as.data.frame(VCF_WT_E_B2@rowRanges)
Names_B <- rownames(B)
B <- cbind(B, Names_B)

View (B)

#SNP in double B1##

SNP_WT_E_B3 <- Rsubread::exactSNP (readFile = "C://Rawdata (1)/OsM2/0.2-5/Sorted_WTB3-0_2-0_5.rap.bam", isBAM = T,
                                refGenomeFile = "D://IGV/RAP REF GENOME/IRGSP-1.0_genome.fasta", outputFile ="SNP_WT_B3.VCF",
                                minAllelicFraction = 0.4)

VCF_WT_E_B3 <- VariantAnnotation::readVcf("SNP_WT_B3.VCF")
C <- as.data.frame(VCF_WT_E_B3@rowRanges)
Names_C <- rownames(C)
C <- cbind(C, Names_C)
View (C) 


common_WT_B1B2 <- inner_join(A, B)
common_WT_B1B2 <- distinct(common_WT_B1B2)

common_WT_B1B2B3 <- inner_join(C, common_WT_B1B2)
common_WT_B1B2B3 <- distinct(common_WT_B1B2B3)
common_WT_B1B2B3_SNP_RES <- common_WT_B1B2B3%>% filter(common_WT_B1B2B3$width<2)


View (common_WT_B1B2B3_SNP_RES)
writexl::write_xlsx(as.data.frame (common_WT_B1B2B3_SNP_RES$Names_A), "SNP_RES_WT_E.xlsx")
writexl::write_xlsx(as.data.frame (common_WT_B1B2B3_SNP_RES[, c(7, 1:4)]), "SNP_RES_WT_E_full.xlsx")

#Variation in TP309 0.2-0.5 cm WT infl RNA-Seq but not in SNP_SEEK#

SNP_SEEK_TP309_E_WT <- read_tsv("SNP_SEEK_TP309_E_WT.txt")%>%
  tidyr::separate(col = 'Not SNP positions: 2877', into = paste0(c("Chr", "position")), sep = ",")

SNP_SEEK_TP309_E_WT$position<- gsub(")","",as.character(SNP_SEEK_TP309_E_WT$position))
SNP_SEEK_TP309_E_WT$Chr<- substring(as.character (SNP_SEEK_TP309_E_WT$Chr),2)
SNP_SEEK_TP309_E_WT$position <- as.numeric(SNP_SEEK_TP309_E_WT$position)
SNP_SEEK_TP309_E_WT <- inner_join(SNP_SEEK_TP309_E_WT, common_WT_B1B2B3_SNP_RES, by = c("position"= "start"))
SNP_SEEK_TP309_E_WT <- SNP_SEEK_TP309_E_WT [, c (1, 2, 8)]
SNP_SEEK_TP309_E_WT$SNP <-SNP_SEEK_TP309_E_WT$Names_C
SNP_SEEK_TP309_E_WT <- SNP_SEEK_TP309_E_WT [, c (1, 2, 4)]
