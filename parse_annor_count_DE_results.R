# Parse Diff Exp results
# based on annotation and abundance parse diffExp results 
# for annotation select uniprot id

rm(list = ls()) # Limpiar la memoria de la sesion de R


options(stringsAsFactors = FALSE) # 

library(tidyverse)

url <- 'https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R'

source(url)


# Annot
path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

trino_file <- list.files(path = path, pattern = 'xls$', full.names = T)
x <- data.table::fread(trino_file, sep = '\t', na.strings = '.')

load(paste0(path, '/count_annot_multiple_contrast_vs_control_up_down_genes.Rdata'))

# data 

data %>% as_tibble(rownames = 'ids') -> data

# prepare query ids
parse_results <- function(df) {
  
  llist <- function(x) {
    x <- paste(x, sep = ';', collapse = ';')
    x <- list(x)
    x <- unlist(x)
  }
  
  df %>% 
    dplyr::select(ids, sampleB) %>% 
    group_by(ids) %>%
    summarise(Freq = length(sampleB), 
      across(sampleB, .fns = llist)) -> distinct_trans
  
  
  # prepare exclusive results 
  
  distinct_trans %>% filter(Freq == 1) %>%
    left_join(df %>% 
        dplyr::select(ids, padj, log2FoldChange), 
      by = "ids") -> distinct_only
  
  distinct_trans %>% pull(ids) %>% unique()-> query.ids
  
  str(query.ids)
  
  sum(keep <- x$transcript_id %in% query.ids)
  
  annot <- x[keep,]
  
  # length(unique(annot$transcript_id))
  # 
  # length(query.ids)
  
  # split trinotate results
  
  blastx <- split_blast(annot, "sprot_Top_BLASTX_hit")
  
  
  blastx %>% 
    mutate(orf = ifelse(is.na(protein), 
      FALSE, TRUE)) %>%
    group_by(transcript, orf) %>%
    filter(identity == max(identity)) %>%
    distinct(transcript, uniprot, 
      name, genus, identity) -> blastx_unique
  
  
  distinct_trans %>% 
    left_join(distinct_only) %>%
    left_join(blastx_unique, 
      by = c('ids'='transcript')) %>%
    left_join(data) -> savedf
  
  return(savedf)
  
  
}

savedf_down <- parse_results(down_df)
savedf_up <- parse_results(up_df)

write_excel_csv(savedf_down,
  file = paste0(path, '/Cancer_vs_Control_annot_count_down_genes.xls'))

write_excel_csv(savedf_up,
  file = paste0(path, '/Cancer_vs_Control_annot_count_up_genes.xls'))

# sblastx <- trinotateR::summary_blast(blastx)
# 
# length(unique(blastx$transcript))
# 
# table(blastx$domain)

blastp <- split_blast(annot, "sprot_Top_BLASTP_hit")
sblastp <- trinotateR::summary_blast(blastp)

length(unique(blastp$transcript))



# omit () -----
# UniProt.ws: A package for retrieving data from the UniProt web service
# BiocManager::install("UniProt.ws")
# convert 

# uniprots <- blastp %>% 
#   filter(genus %in% 'Homo') %>% 
#   mutate(uniprot = sapply(strsplit(uniprot, "_"), "[", 1)) %>%
#   pull(uniprot) %>% unique()
# 
# 
# str(uniprots)

# library(UniProt.ws)
# 
# availableUniprotSpecies("Homo sapiens")
# lookupUniprotSpeciesFromTaxId(taxId=9606)
# No data is available for the keys provided.
# up <- UniProt.ws(taxId=9606)
# select(up, head(uniprots), "ENTREZ_GENE")
