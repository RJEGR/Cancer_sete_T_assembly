

rm(list = ls()) # Limpiar la memoria de la sesion de R
if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE) # 

library(DESeq2)
library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

rds_f_l <- list.files(path, pattern = 'CONTRAST_',  full.names = TRUE)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

go_f <- list.files(path, pattern = 'Trinotate_report.xls.gene_ontology',  full.names = TRUE)

swiss_f <- list.files(path, pattern = 'Trinotate.xls.blastx.tsv',  full.names = TRUE)
# Database source ----

MAP <- topGO::readMappings(go_f)

swissdf <- read_tsv(swiss_f) %>%
  mutate(orf = ifelse(is.na(protein), FALSE, TRUE)) %>%
  group_by(transcript, orf) %>%
  filter(identity == max(identity)) %>%
  ungroup() %>%
  distinct(transcript, uniprot, identity, name, genus, orf) %>%
  dplyr::rename("transcript_id"="transcript", "protein_name"="name")

# DB_f <- list.files(path = path, pattern = 'Trinotate.xls$', full.names = T)
# annot <- data.table::fread(DB_f, sep = '\t', na.strings = '.')
# names(annot)[1] <- 'gene_id'


get_res <- function(dds, contrast, alpha_cutoff = 0.1) {
  
  require(DESeq2)
  require(tidyverse)
  require(IHW)
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast, alpha = alpha_cutoff)
  
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
  
 
  res <- res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("transcript_id|sample|pvalue|padj")),
      round ,digits = 2)
  
  #  take the set of p-values and then calculate the adjusted p-values using a independent Hypothesis Weighting
  
  ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = 0.1) # toma menos de 1 minuto.
  
  weighted_pvalue <- as.data.frame(ihwRes)$weighted_pvalue
  
  
  res %>% cbind(., weighted_pvalue = weighted_pvalue) %>% as_tibble(rownames = "transcript_id")
}

dds2res <- function(dds) {
  contrast <- levels(colData(dds)$Design)
  
  out <- list()
  
  for(i in 2:length(contrast)) {
    j <- i
    cat('\nContrast', contrast[1],' and ', contrast[j],'\n')
    res <- get_res(dds, contrast[c(1,j)])
    out[[j]] <- res
    
  }
  
  do.call(rbind, out) -> res
  
  return(res)

}

prep_annot_inf <- function(query.ids, annot, MAP) {
  
  n <- sum(keep <- names(MAP) %in% query.ids) / length(query.ids) 
  
  cat('\n % of ids mapped in Gen Ontology: ', n*100, '%\n')
  
  STRG2GO <- data.frame(transcript_id = rep(names(MAP[keep]),
    sapply(MAP[keep], length)),
    GO.ID = unlist(MAP[keep]), row.names = NULL) %>% as_tibble()
  
  # Sanity check: % of recovery querys with go ids
  # sum(unique(STRG2GO$ids) %in% query.ids) / length(query.ids) 
  
  llist <- function(x) {
    x <- paste(x, sep = ';', collapse = ';')
    x <- list(x)
    x <- unlist(x)
  }
  
  STRG2GO %>%
    group_by(transcript_id) %>%
    summarise(across(GO.ID, .fns = llist), .groups = "drop_last") %>% 
    ungroup() -> STRG2GO
  
  out <- data.frame(transcript_id = query.ids) %>% 
    left_join(STRG2GO) %>% 
    left_join(annot) %>%
    as_tibble()
  
  return(out)
  
}

prep_multiple_contrast_output <- function(df, select_exclusive = T) {
  
  # names(df) should contain at least cols:
  # "ids"            
  # "sampleB"       
  #  "log2FoldChange" 
  # "pvalue" or    "padj" 
  # "lfcT"
  # Ex: res or res.p output
  
  
  llist <- function(x) {
    x <- paste(x, sep = ';', collapse = ';')
    x <- list(x)
    x <- unlist(x)
  }
  
  df %>% # this is the only input needed
    dplyr::select(transcript_id, sampleB) %>%
    mutate(sampleB = gsub("[[:blank:]]","_", sampleB)) %>%
    group_by(transcript_id) %>%
    summarise(n_sam = length(sampleB), 
      across(sampleB, .fns = llist)) -> distinct_trans
  
  distinct_trans %>% group_by(n_sam) %>% tally() -> prevalence
  
  cat("\n Multiple contrast Prevalence: \n",
    "\t",paste0(prevalence$n_sam, 'sam:'), "\n", 
    "\t",prevalence$n, "\n")
  
  
  if(select_exclusive) {
    
    cat("\n Using transcripts detected by 1 sample (1sam) \n")
    
    # Select only exclusive transcripts per contrast
    
    distinct_trans %>% 
      filter(n_sam == 1) %>% 
      mutate(sampleB = unlist(sampleB)) %>%
      pull(transcript_id) -> distinct_trans_query
    
    
    
    
  } else {
    df %>% dplyr::select(transcript_id, sampleB) %>%
      group_by(sampleB) %>% distinct(transcript_id) -> distinct_trans
    
    cat("\n Prevalence: ",paste0("\n",prevalence$sampleB, ':\t', prevalence$n), "\n")
  }
  
  # 2) add annot info to the distinct transcript dataset ----
  
  cat("\n Getting annot info to the distinct transcript dataset\n")
  
  query.ids <- distinct_trans_query
  
  prep_annot_inf(query.ids, MAP = MAP, annot = swissdf) %>%
    # Add p values to the distinct trans
    left_join(df, by = 'transcript_id') -> out
  
  
  cat("\n Creating a set of", dim(out) , "\n")
  
  return(out)
  
}


# contrast2 <- levels(colData(dds)$Design)

# res <- get_res(dds, contrast2[c(1,2)])


# 8) loop the contrasts results ----


rds_f <- rds_f_l[1]

dds <- read_rds(rds_f)

res <- dds2res(dds)

alpha = 0.05; lfcThreshold = 2

# 9) Filter significant transcripts ----

res.p <- res %>% 
  mutate(filter_col = NA) %>% 
  mutate(filter_col = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, 'sigfc', filter_col)) %>%
  drop_na(filter_col) %>% select(-filter_col) %>%
  select(-lfcSE, -stat)

# res.p %>% count(sampleB) %>% view()

# str(query.ids <- unique(res.p$transcript_id))

# prep_annot_inf(query.ids, MAP = MAP, annot = swissdf)


# 10) Generate single / exclusive transcripts per condition ----

res.p.out <- prep_multiple_contrast_output(res.p)

# res.p.out %>% count(sampleB) %>% view()

# 11) Bind data count to output ----
# vst <- vst(dds) 

count <- as_tibble(assay(dds), rownames = "transcript_id")

nrow(res.p.out <- res.p.out %>% left_join(count))

res.p.out <- res.p.out %>% arrange(sampleB) 

# 12) Write output ----

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"

file_name <- strsplit(basename(rds_f), "_")[[1]][1:2]
file_name <- paste0(out_path, paste(file_name,collapse = '_'),'_Annot_count_down_up_genes.xls')

up <- 'Up-regulated' # transcript w/ log2FoldChange > lfcThreshold
down <- 'Down-regulated' # transcript w/ log2FoldChange < lfcThreshold

res.p %>%
  mutate(lfcT = 'Basal') %>%
  mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, down, lfcT)) %>%
  mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, up, lfcT)) %>%
  group_by(sampleB) %>%
  count(lfcT) %>% view()

res.p.out %>%
  mutate(lfcT = 'Basal') %>%
  mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, down, lfcT)) %>%
  mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, up, lfcT)) %>%
  group_by(sampleB) %>%
  count(lfcT) %>% view()


# xlsx::write.xlsx(res.p.out, file = file_name, sheetName = "RESULTS", row.names = FALSE)
write_excel_csv(res.p.out, file = file_name)

# CONCAT IN A SINGLE BATCH -----

dds_f <- list.files(path = "~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/",
  pattern ="CONTRAST_", full.names = T)

bind_data_sources <- function(x, alpha = 0.05, lfcThreshold = 2) {
  
  rds_f <- x
  
  dds <- read_rds(rds_f)
  
  res <- dds2res(dds)
  
  # 9) Filter significant transcripts ----
  
  res.p <- res %>% 
    mutate(filter_col = NA) %>% 
    mutate(filter_col = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, 'sigfc', filter_col)) %>%
    drop_na(filter_col) %>% select(-filter_col) %>%
    select(-lfcSE, -stat)
  
  # 10) Generate single / exclusive transcripts per condition ----
  
  res.p.out <- prep_multiple_contrast_output(res.p)
  
  
  # 11) Bind data count to output ----
  # vst <- vst(dds) 
  
  # count <- as_tibble(assay(dds), rownames = "transcript_id")
  
  # nrow(res.p.out <- res.p.out %>% left_join(count))
  
  res.p.out <- res.p.out %>% arrange(sampleB) 
  
  return(res.p.out)
  
}

out <- lapply(dds_f, bind_data_sources)

res.p <- do.call(rbind, out)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'
count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

dim(raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))
raw_count <- as_tibble(raw_count, rownames = "transcript_id")
nrow(res.p.out <- res.p %>% left_join(raw_count))

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"

file_name <- "ALL_MULTIPLE_CONTRAST"
file_name <- paste0(out_path, paste(file_name,collapse = '_'),'_Annot_count_down_up_genes.xls')


# end ----

# calculate the adjusted p-values using ihw ----
#  take the set of p-values and then calculate the adjusted p-values using a independent Hypothesis Weighting
library(IHW)

# contrast_bk <- levels(colData(dds)$Design)

# res <- get_res(dds, contrast_bk[c(1,2)])

res %>% rstatix::cor_mat(vars = c("pvalue", "padj", "weighted_pvalue"))

# ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = 0.1) # toma menos de 1 minuto.
# weighted_pvalue <- as.data.frame(ihwRes)$weighted_pvalue

# ihwResDf %>% rstatix::cor_mat(vars = c("pvalue", "adj_pvalue", "weighted_pvalue"))
# ggplot(res, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)

# Data tidying: ----
# cols: ids, Freq (optional), sampleA, sampleB, pvalue (adj), log2fc, uniprot id, identity, seq_legnth, name, henus, orf, sample_count_matrix.

# REVISAR SECCION # 1) Separate exclusive transcripts DEL SCRIPT FUNCTIONAL ANNOT..
# REVISAR SCRIPT PARSE ANNOT & functions.R PARA DETERMINAR  prep_dist_data()



# 