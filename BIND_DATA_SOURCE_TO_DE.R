

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

# view(swissdf)

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


which_gene <- res[which.min(res$padj),]$transcript_id

d1 <- plotCounts(dds, gene=which.min(res$padj), intgroup="Design", returnData=TRUE)
d2 <- plotCounts(dds, gene=which.max(res$padj), intgroup="Design", returnData=TRUE)


which_LF <- res %>% 
  filter(transcript_id %in% which_gene) %>%
  pull(log2FoldChange)

DF <- rbind(d1,d2) %>% 
  as_tibble(rownames = "LIBRARY_ID") %>%
  left_join(as_tibble(colData(dds) ))

# 
DF %>%
  # drop_na(CONTRASTE_D) %>%
  ggplot(aes(x=Design, y=count)) +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 5, alpha = 0.5) +
  # scale_y_log10() +
  stat_summary(fun.data = "mean_cl_boot", colour = "red", linewidth = 0.7, size = 1) +
  labs(y = paste0("Gene count (",which_gene ,")"), x = "Treatment") +
  theme_bw(base_family = "GillSans", base_size = 20)

# 9) Filter significant transcripts ----

res.p <- res %>% 
  mutate(filter_col = NA) %>% 
  mutate(filter_col = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, 'sigfc', filter_col)) %>%
  drop_na(filter_col) %>% select(-filter_col)

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
  # mutate(lfcT = 'Basal') %>%
  # mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, down, lfcT)) %>%
  # mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, up, lfcT)) %>%
  mutate(lfcT = sign(log2FoldChange)) %>%
  group_by(sampleB) %>%
  dplyr::count(lfcT)

res.p.out %>%
  mutate(lfcT = 'Basal') %>%
  mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, down, lfcT)) %>%
  mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, up, lfcT)) %>%
  group_by(sampleB) %>%
  dplyr::count(lfcT) %>% view()


# xlsx::write.xlsx(res.p.out, file = file_name, sheetName = "RESULTS", row.names = FALSE)
write_excel_csv(res.p.out, file = file_name)

# from file above, use:

res.p.out

# global presence of relevant proteins -----

res.p.out <- read_csv("~/Documents/DOCTORADO/human_cancer_dataset/ALL_MULTIPLE_CONTRAST_Annot_count_down_up_genes.xls")


read_xls <- function(x) { readxl::read_excel(y, sheet = x)}

y <- "~/Documents/DOCTORADO/human_cancer_dataset/ALL_MULTIPLE_CONTRAST_Annot_count_down_up_genes_LOGFC_1.xlsx"

res.p.out <- lapply(readxl::excel_sheets(y), read_xls)

res.p.out <- do.call(rbind, res.p.out)


which_histo <- res.p.out %>%  
  filter(grepl("histocompatibility", protein_name)) %>% 
  distinct(uniprot) %>% pull() %>%
  paste0("^",., collapse = "|")

which_vav <- res.p.out %>%
  filter(grepl("guanine nucleotide exchange factor ", protein_name)) %>%
  distinct(uniprot) %>% pull() %>%
  paste0("^",., collapse = "|")

which_chemok <- res.p.out %>%
  filter(grepl("chemokine", protein_name)) %>%
  distinct(uniprot) %>% pull() %>%
  paste0("^",., collapse = "|")

which_tcrep <- res.p.out %>%
  filter(grepl("T cell receptor beta ", protein_name)) %>%
  distinct(uniprot) %>% pull() %>%
  paste0("^",., collapse = "|")

which_prot <- c("VAV1", "CCL28", "ITK", "CXC19", "IL6", "JAK", "RAC2", "NCF1", "STAT1", "IL1A", "CSF1", "CDC42")

which_prot <- paste0("^",which_prot, collapse = "|")

which_prot <- paste0(which_histo, which_vav, which_chemok, which_tcrep, which_prot)

# swissdf %>% view()

# res.p.out %>% mutate(lfcT = sign(log2FoldChange)) %>% group_by(sampleA, sampleB) %>%  dplyr::count(lfcT) %>% view()

which_genes <- res.p.out %>%  
  filter(grepl(which_prot, uniprot))
  # distinct(transcript_id, .keep_all = T)

str(query.genes <- which_genes %>% distinct(transcript_id) %>% pull())

# res.p.out %>% filter(transcript_id %in% query.genes) %>% view()

out <- list()

for (i in query.genes) {
  
  transcript_id <- i
  
  # which(res$transcript_id == gene)
  
  df <- plotCounts(dds, gene=transcript_id, intgroup="Design", returnData=TRUE, normalized = T, transform = F)
  
  df <- data.frame(df, transcript_id) %>% as_tibble(rownames = "LIBRARY_ID") 
  
  out[[i]] <-  df
}

head(DF <- do.call(rbind, out))

DF <- DF %>%
  left_join(res.p.out %>% distinct(transcript_id, uniprot)) %>%
  left_join(as_tibble(colData(dds) ))

recode_to <- c(`Control` = "Control", `CON_CANCER` = "Cancer")

LOGFC_LABELLER <- which_genes  %>% 
  separate(uniprot,into = c("uniprot", "genus"), sep = "_") %>%
  arrange(log2FoldChange) %>%
  mutate(log2FoldChange = paste0(uniprot, " (",log2FoldChange,")")) %>%
  select(uniprot, log2FoldChange) %>%
  pull(log2FoldChange, uniprot)

p <- DF %>%
  dplyr::mutate(Design = dplyr::recode_factor(Design, !!!recode_to)) %>%
  separate(uniprot,into = c("uniprot", "genus"), sep = "_") %>%
  mutate(uniprot = factor(uniprot, levels = unique(names(LOGFC_LABELLER)))) %>%
  ggplot(aes(x=Design, y=log2(count+0.5), group = uniprot)) +
  facet_wrap(~ uniprot, labeller = labeller(uniprot = LOGFC_LABELLER)) +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 1, alpha = 0.5) +
  # scale_y_log10() +
  stat_summary(fun.data = "mean_cl_boot", colour = "red", linewidth = 0.7, size = 0.7, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line", colour = "red") +
  labs(y = "Gene count (Log2)", x = "Assay") +
  theme_bw(base_family = "GillSans", base_size = 14)

ggsave(p, filename = "DESEQ2LINEPLOT.png", 
  path = path, width = 8.5, height = 7.5, device = png, dpi = 300)


DF %>%
  dplyr::mutate(Design = dplyr::recode_factor(Design, !!!recode_to)) %>%
  separate(uniprot,into = c("uniprot", "genus"), sep = "_") %>%
  # mutate(uniprot = factor(uniprot, levels = unique(names(LOGFC_LABELLER)))) %>%
  dplyr::filter(uniprot == "TRBC2") %>%
  ggplot(aes(x=Design, y=count, group = uniprot)) +
  facet_wrap(~ uniprot, labeller = labeller(uniprot = LOGFC_LABELLER)) +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 3, alpha = 0.5) +
  # scale_y_log10() +
  stat_summary(fun.data = "mean_cl_boot", colour = "red", linewidth = 0.7, size = 1, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line", colour = "red", linewidth = 1) +
  labs(y = "Gene count (Normalized)", x = "Assay") +
  theme_bw(base_family = "GillSans", base_size = 16)

# AS HEATMAP

vst <- DESeq2::vst(dds) # vst if cols > 10

M <- assay(vst)
keep <- rownames(M) %in% query.genes
head(M <- M[keep,])

sample_cor = cor(M, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(sample_cor, method='euclidean')
hc_samples = hclust(sample_dist, method='complete')


hc_order <- hc_samples$labels[hc_samples$order]

# heatmap(sample_cor, col = cm.colors(12))
colData <- colData(dds) %>% as_tibble()

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(colData) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> sample_cor_long

library(ggh4x)

p <- sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  # geom_raster() + 
  # geom_text(aes(label = cor), color = 'white') +
  scale_fill_viridis_c(option = "C", name = "Pearson", direction = -1) +
  scale_x_discrete(position = 'top') +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::scale_y_dendrogram(hclust = hc_samples) +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = ggh4x::guide_axis_manual(labels = hc_order, label_size = 7)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust = -0.15, vjust = 1),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) # -> pheat


ggsave(p, filename = "DESEQ2HEATMAP.png", 
  path = path, width = 7, height = 7, device = png, dpi = 300)

genes_cor = cor(t(M), method='pearson', use='pairwise.complete.obs')
genes_dist = dist(genes_cor, method='euclidean')
hc_genes = hclust(genes_dist, method='complete')


genes_order <- hc_genes$labels[hc_genes$order]


recode_to <- res.p.out %>% filter(transcript_id %in% genes_order) %>% 
  distinct(transcript_id, uniprot, protein_name) %>%
  separate(uniprot,into = c("uniprot", "genus"), sep = "_")

# recode_to <- structure(recode_to$uniprot, names = recode_to$transcript_id)

recode_to <- structure(recode_to$protein_name, names = recode_to$transcript_id)

identical(sort(names(recode_to)),sort(genes_order))

genes_order <- recode_to[match(genes_order, names(recode_to))]

identical(names(genes_order),  hc_genes$labels[hc_genes$order])


MLONG <- M %>% 
  as_tibble(rownames = 'transcript_id') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'Reads', names_to = "LIBRARY_ID") %>%
  left_join(colData) %>% 
  left_join(res.p.out %>% distinct(transcript_id, uniprot)) %>%
  separate(uniprot,into = c("uniprot", "genus"), sep = "_") %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order))


p <- MLONG %>%
  ggplot(aes(x = LIBRARY_ID, y = transcript_id, fill = log2(Reads))) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "C", name = "Gene count (Log2)", direction = -1) +
  scale_x_discrete(position = 'top') +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  guides(y.sec = ggh4x::guide_axis_manual(labels = genes_order, label_size = 7)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = -0.15, vjust = 1, size = 5),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank())

ggsave(p, filename = "DESEQ2HEATMAP_GENES.png", 
  path = path, width = 7, height = 7, device = png, dpi = 300)


# CONCAT IN A SINGLE BATCH -----

dds_f <- list.files(path = "~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/",
  pattern ="CONTRAST_", full.names = T)

bind_data_sources <- function(x, alpha = 0.05, lfcThreshold = 1) {
  
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

names(out) <-  gsub("_DDS.rds", "", basename(dds_f))

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'
count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

dim(raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))
raw_count <- as_tibble(raw_count, rownames = "transcript_id")
nrow(res.p.out <- res.p %>% left_join(raw_count))

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"

file_name <- "ALL_MULTIPLE_CONTRAST"
file_name <- paste0(out_path, paste(file_name,collapse = '_'),'_Annot_count_down_up_genes_LOGFC_1.xls')


# xlsx::write.xlsx(res.p.out, file = file_name, sheetName = "RESULTS", row.names = FALSE)

write_excel_csv(res.p.out, file = file_name)

xlsx::write.xlsx(res.p.out, file = file_name, sheetName = "RESULTS", row.names = FALSE)

library(xlsx)

wb <- createWorkbook()
# sheetnames <- paste0("Sheet", seq_along(out)) # or names(datas) if provided
sheetnames <- names(out)

sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, out, sheets)
saveWorkbook(wb, file = file_name)

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