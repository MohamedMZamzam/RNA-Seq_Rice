##annot files tracking##
library (tidyverse)
Anno_OBASE_15.11.22 <- read_tsv ("C:/Rawdata (1)/annotation_O.BASE/OryzabaseGeneListAll_20221115010001.txt")
Anno_OBASE_10.12.22 <- read_tsv ("C:/Rawdata (1)/annotation_O.BASE/OryzabaseGeneListAll_20221210010000.txt")

New_info <- setdiff(Anno_OBASE_10.12.22, Anno_OBASE_15.11.22)


de_early <- readxl::read_xlsx("C:/transgenic_plants_ analyses/osm2_ko/Received from Raghavram/RNA seq/Okayed list 3vs 2 email 13.10.22/DE_early_3wt_2mut.xlsx")
de_early <- filter(de_early, abs (de_early$logFC) >= 0.5849625, de_early$FDR <= 0.05)
de_early_new_info <- inner_join(de_early, New_info, by = c ("rownames(l)" = "RAP ID"))


de_middle <- readxl::read_xlsx("C:/transgenic_plants_ analyses/osm2_ko/Received from Raghavram/RNA seq/Okayed list 3vs 2 email 13.10.22/DE_Middle_3wt_2mut.xlsx")
de_middle <- filter(de_middle, abs (de_middle$logFC) >= 0.5849625, de_middle$FDR <= 0.05)
de_middle_new_info <- inner_join(de_middle, New_info, by = c ("rownames(l)" = "RAP ID"))


de_Late <- readxl::read_xlsx("C:/transgenic_plants_ analyses/osm2_ko/Received from Raghavram/RNA seq/Okayed list 3vs 2 email 13.10.22/DE_Late_3wt_2mut.xlsx")
de_Late <- filter(de_Late, abs (de_Late$logFC) >= 0.5849625, de_Late$FDR <= 0.05)
de_Late_new_info <- inner_join(de_Late, New_info, by = c ("rownames(l)" = "RAP ID"))
---------------
  
#New annotation####

Early <- readxl::read_xlsx("C:/transgenic_plants_ analyses/osm2_ko/Received from Raghavram/RNA seq/Okayed list 3vs 2 email 13.10.22/DE_early_3wt_2mut.xlsx")
Early__2FC <- filter(Early, abs (Early$logFC) >= 1, Early$FDR <= 0.05)
Early__2FC <- left_join(Early__2FC, Anno_OBASE_10.12.22, by = c ("rownames(l)" = "RAP ID"))

writexl::write_xlsx(Early__2FC, "Dataset2.xlsx")


Middle <- readxl::read_xlsx("C:/transgenic_plants_ analyses/osm2_ko/Received from Raghavram/RNA seq/Okayed list 3vs 2 email 13.10.22/DE_middle_3wt_2mut.xlsx")
Middle__2FC <- filter(Middle, abs (Middle$logFC) >= 1, Middle$FDR <= 0.05)
Middle__2FC <- left_join(Middle__2FC, Anno_OBASE_10.12.22, by = c ("rownames(l)" = "RAP ID"))
writexl::write_xlsx(Middle__2FC, "Dataset2_1.xlsx")


Late <- readxl::read_xlsx("C:/transgenic_plants_ analyses/osm2_ko/Received from Raghavram/RNA seq/Okayed list 3vs 2 email 13.10.22/DE_Late_3wt_2mut.xlsx")
Late__2FC <- filter(Early, abs (Late$logFC) >= 1, Late$FDR <= 0.05)
Late__2FC <- left_join(Late__2FC, Anno_OBASE_10.12.22, by = c ("rownames(l)" = "RAP ID")) 
writexl::write_xlsx(Late__2FC, "Dataset2_2.xlsx")

