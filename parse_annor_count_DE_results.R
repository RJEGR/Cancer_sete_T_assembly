# Parse Diff Exp results
# based on annotation and abundance parse diffExp results 
# for annotation select uniprot id

rm(list = ls()) # Limpiar la memoria de la sesion de R


options(stringsAsFactors = FALSE) # 

library(tidyverse)

# url <- 'https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R'

# source(url)


# Annot
path <- "~/Downloads/CORAL/"

trino_file <- list.files(path = path, pattern = 'xls$', full.names = T)

x <- data.table::fread(trino_file, sep = '\t', na.strings = '.')
names(x)[1] <- 'gene_id'

de_file <- "RSEM.isoform.counts.matrix.Boquita_vs_Carrizales.DESeq2.DE_results.P0.01_C1.DE.subset$"

de_file <- list.files(path = path, pattern = de_file, full.names = T)


# Prepare DE data ----

# negative lFC == Carrizales
# Positive lFC == Boquita

res <- read.table(de_file, header=T, com='', row.names=1, check.names=F, sep='\t', stringsAsFactors = FALSE) %>%
  as_tibble(rownames = 'transcript')


res %>% distinct(sampleA, sampleB)

res %>% arrange(log2FoldChange)

alpha = 0.01 

lfcThreshold = 2

# Oscar: The R package DESeq2 (Love et al. 2014) was implemented to identify differentially expressed genes (DEGs) between the groups (false discovery rate < 0.01, fold change > 2, TPM > 1).

res %>% 
  mutate(cc = NA) %>% 
  mutate(cc = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, 'sigfc', cc)) %>%
  drop_na(cc) -> res.p

up <- 'Up-regulated' # transcript w/ log2FoldChange > lfcThreshold
down <- 'Down-regulated' # transcript w/ log2FoldChange < lfcThreshold

res.p %>%
  mutate(lfcT = 'Basal') %>%
  mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, down, lfcT)) %>%
  mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, up, lfcT)) -> res.p

up_df <- res.p %>% filter(grepl(up, lfcT))

down_df <- res.p %>% filter(grepl(down, lfcT))

save(up_df, down_df, 
  file = paste0(path, 'count_annot_Boquita_vs_Carrizales_up_down_genes.Rdata'))

# 2) Prepare annotation data ----

# data 

x %>% as_tibble(rownames = 'ids') -> data

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

blast <- split_blast(x, "sprot_Top_BLASTP_hit")

length(unique(blast$transcript))
# 
table(blast$domain)

x

