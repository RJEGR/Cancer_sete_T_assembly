# SEARCHING MALIGNAN, METASTASIS, AND SURVIVAL GENES
# RICARDO GOMEZ-REYES
# LOAD SWISSPROT AND FILTER
# LOAD DE RESULTS AND BIND FILTERED GENES W/ DETAILS ABOUT UP/DOWN EXPRESSION (BY CONTRAST C OR D)

# REVISE BIND_DATA_SOURCE_TO_DE.R 
# GENERATE SOME PLOTS

rm(list = ls()) # Limpiar la memoria de la sesion de R

if(!is.null(dev.list())) dev.off()


WHICH_PROTEINS <- c('PALB2_HUMAN',	'BRCA1_HUMAN',	'BRCA2_HUMAN',	'CHK2_HUMAN',	'CADH1_HUMAN',	'PTEN_HUMAN',	'STK11_HUMAN',	'P53_HUMAN',	'ATM_HUMAN',	'BARD1_HUMAN',	'FANCJ_HUMAN','BACH1_HUMAN','CASP8_HUMAN',	'CTLA4_HUMAN',	'CP19A_HUMAN',	'FGFR2_HUMAN',	'H19_HUMAN',	'LSP1_HUMAN',	'M3K1_HUMAN',	'MRE11_HUMAN',	'NBN_HUMAN',	'RAD51_HUMAN',	'TERT_HUMAN',	'LOX15_HUMAN',	'CO4A6_HUMAN',	'LMB13_HUMAN',	'MTAP_MACMU',	'PA24A_HUMAN',	'ATTY_HUMAN')


options(stringsAsFactors = FALSE) # 

library(DESeq2)

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

rds_f_l <- list.files(path, pattern = 'CONTRAST_',  full.names = TRUE)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

swiss_f <- list.files(path, pattern = 'Trinotate.xls.blastx.tsv',  full.names = TRUE)

# Database source ----

swissdf <- read_tsv(swiss_f) %>%
  mutate(orf = ifelse(is.na(protein), FALSE, TRUE)) %>%
  group_by(transcript, orf) %>%
  filter(identity == max(identity)) %>%
  ungroup() %>%
  distinct(transcript, uniprot, identity, name, genus, orf) %>%
  dplyr::rename("transcript_id"="transcript", "protein_name"="name")


# FILTER

length(WHICH_PROTEINS)

sum(keep <- WHICH_PROTEINS %in% unique(sort(swissdf$uniprot)))

WHICH_PROTEINS[!keep]

bind_ids <- function(x) {
  x <- unique(x)
  x <- paste(x, sep = ';', collapse = ';')
}

swissdf %>% 
  filter(uniprot %in% WHICH_PROTEINS) %>% 
  group_by(uniprot, protein_name) %>%
  summarise(across(transcript_id, .fns = bind_ids), .groups = "drop_last") 

ANNOT <- swissdf %>% 
  filter(uniprot %in% WHICH_PROTEINS) 


# BIND TO DE RESULTS

path <- "~/Documents/DOCTORADO/human_cancer_dataset/"

file_name <- "ALL_MULTIPLE_CONTRAST_Annot_count_down_up_genes_LOGFC_1.xlsx"

file_name <- paste0(path, file_name)

library(readxl)

print(DF1 <- read_excel(file_name) %>% filter(uniprot %in% WHICH_PROTEINS))

# Solo 5 fueron significativas en la comparacion control vs muestras con Cancer

# LOADING CONTRAST C AND D ONLY (W/O FILTERING)



CNT_C <- rds_f_l[grepl("CONTRAST_C_DDS", rds_f_l)]
CNT_D <- rds_f_l[grepl("CONTRAST_D_DDS", rds_f_l)]

CNT_C <- read_rds(CNT_C)
CNT_D <- read_rds(CNT_D)

CNT_C <- dds2res(CNT_C)
CNT_D <- dds2res(CNT_D)

write_rds(list(CNT_C, CNT_D), file = paste0(path, "CONTRAST_C_AND_D_FOR_PUB.rds"))

OUT1 <- CNT_C %>% right_join(ANNOT) %>% select(-genus)
OUT2 <- CNT_D %>% right_join(ANNOT) %>% select(-genus)


write_excel_csv(OUT1, file = paste0(path, "/CONTRAST_C_GRADOS_HISTOLOGICOS.xls"))

write_excel_csv(OUT2, file = paste0(path, "/CONTRAST_D_METASTASIS_NO_METASTASIS.xls"))

bind_samples <- function(x) {
  x <- unique(x)
  x <- paste(x, collapse = ' and ')
}

recode_to <- c(`-1` = "Up in Sample B", `1` = "Down in sample B")


OUT1 %>% 
  filter(padj < 0.05) %>%
  mutate(log2FoldChange = sign(log2FoldChange)) %>%
  dplyr::mutate(log2FoldChange = recode_factor(log2FoldChange, !!!recode_to)) %>%
  # distinct(uniprot, sampleB) %>% 
  # right_join(ANNOT) %>%
  drop_na(sampleB) %>% 
  group_by(uniprot, log2FoldChange) %>%
  summarise(n = n(), across(sampleB, .fns = bind_samples), .groups = "drop_last")  %>%
  view()

  
OUT2 %>% 
  filter(padj < 0.05) %>%
  mutate(log2FoldChange = sign(log2FoldChange)) %>%
  dplyr::mutate(log2FoldChange = recode_factor(log2FoldChange, !!!recode_to)) %>%
  # distinct(uniprot, sampleB) %>% 
  # right_join(ANNOT) %>%
  drop_na(sampleB) %>% 
  group_by(uniprot, log2FoldChange) %>%
  summarise(n = n(), across(sampleB, .fns = bind_samples), .groups = "drop_last") %>%
  view()



