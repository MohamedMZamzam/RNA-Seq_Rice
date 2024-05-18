#Set working directory
setwd(":/Rawdata (1)/double/")
#Load libraries
library (tidyverse)
library(readxl)
#load the annotation file#
Anno_OBASE <- read_tsv ("C:/Rawdata (1)/annotation_O.BASE/OryzabaseGeneListAll_20221210010000.txt")
#Load the DEG list that you want to annotate#
DEG <-read_xls ("C:/Rawdata (1)/OsM2/0.2-5/DE_early_3wt_2mutant.xls")

#annotate the DEG by joining#
DEG__annotated <- left_join(DEG, Anno_OBASE, by = c ("ID" = "RAP ID")) 

