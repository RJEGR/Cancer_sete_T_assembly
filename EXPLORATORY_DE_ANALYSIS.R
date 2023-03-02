# RICARDO GOMEZ-REYES, FEBRAURY 28 TH, 2023
# PROPOSITO
# DEBIDO A QUE LA RELACION DE IDS OTORGADOS POR CLAUDIA ES DISTINTA A LA ENTREGADA POR FRA, REALIZAREMOS LA EVALUACION DE NOVO


rm(list = ls()) # Limpiar la memoria de la sesion de R
if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE) # 

# library(tidyverse)

# library(patchwork)


path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

mtd_f <- list.files(path, pattern = 'metadata.tsv',  full.names = TRUE)

mtd <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

mtd %>% count(CONTRASTE_A, sort = T) # %>% view()

dim(raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))

# DESEQ2 ----

# colors_fc <- c("red2", "#4169E1", "forestgreen", "grey30")


# Testing CONTRAST A
colData <- mtd

f_col <- "CONTRASTE_A"

names(colData)[names(colData) %in% f_col] <- "Design"

# 0) Remove single samples:

sam <- table(colData$Design)

sam <- names(sam[sam > 1])

colData <- colData %>% drop_na(Design) %>% filter(Design %in% sam)

x <- colData %>% pull(Design, name = LIBRARY_ID)

# 1) According to metadata, filter samples ----

ncol(count <- raw_count[names(raw_count) %in% names(x)])

# 2) Filter data by removing low-abundance genes ----

keep <- rowSums(edgeR::cpm(count) > 1) >= 2

sum(keep) # N transcripts

nrow(count <- raw_count[keep,])

count <- round(count)

# 3) Format metadata and count matrix ----

colData <- colData %>% arrange(match(LIBRARY_ID, names(count)))

sum(names(count) == colData$LIBRARY_ID) # sanity check

colData <- mutate_if(colData, is.character, as.factor)

# using relevel, just specifying the reference level:

colData <- colData %>% mutate(Design = relevel(Design, ref = "Control"))

library(DESeq2)
library(tidyverse)

table(colData$Design)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ Design )

# Note: After estimateSizeFactors(ddsFullCountTable) 
# Evaluate the size factor, coef ad dispPriorVar from different source of design as follow: 
# sizeFactors(dds), head(coef(dds)), attr(dispersionFunction(dds), "dispPriorVar")


# 4) Run DESeq in the following functions order ----

dds <- estimateSizeFactors(ddsFullCountTable) 
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds) # use Wald test when contrasting a unique pair (ex. control vs cancer)

# Note: also run as a single batch 
# dds <- estimateSizeFactors(ddsFullCountTable)
dds <- DESeq(dds)

# 5) Test transformation for downstream analyses----
# In order to test for differential expression, we operate on raw counts and use discrete distributions ... 
# However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data. Ex:

vst <- vst(dds) # wrapped version for varianceStabilizingTransformation
ntr <- DESeq2::normTransform(dds)

# rldds <- rlog(dds, intercept = mcols(vsd)$dispFit) # muy tardado, y vsd da mejor resultado

# 
# DESeq2::plotPCA(ntr, intgroup = "Design")
#
# DESeq2::plotPCA(vst, intgroup = "Design")

# BiocManager::install("vsn")

# VERIFICAR LA SD DE LOS DATOS TRANSFORMADOS, normTransform y vsd. Encontrar que transformacion disminuye la sd (linea roja) de los datos.

# ntr mantiene un sd debajo de 1 para la mayoria de los transcritos (eje x)

# vsn::meanSdPlot(assay(vst))
# vsn::meanSdPlot(assay(dds))
# vsn::meanSdPlot(assay(ntr)) # 

# 6 Generate results ----



get_res <- function(dds, contrast, alpha_cutoff = 0.1, lfc_cutoff = 0) {
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast = contrast, alpha = alpha_cutoff, lfcThreshold = lfc_cutoff)
  
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
  
  res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    as_tibble(rownames = "ids") %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("ids|sample|pvalue|padj")),
      round ,digits = 2)
}

contrast <- levels(colData(dds)$Design)

res <- get_res(dds, contrast) # # take a long


write_rds(dds, file = paste0(path, '/CONTRAST_A_DDS.rds'))

# WHEN CONTRASTING MULTIPLE GROUPS:

# The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT test can be especially helpful when performing time course analyses (not this case)

# 7) multiple contrat from designs against control -----

DESeqDataSetFromInputs <- function(colData, count, f_col = ...) {
  
  names(colData)[names(colData) %in% f_col] <- "Design"
  
  # 0) Remove single samples:
  
  sam <- table(colData$Design)
  
  sam <- names(sam[sam > 1])
  
  colData <- colData %>% drop_na(Design) %>% filter(Design %in% sam)
  
  x <- colData %>% pull(Design, name = LIBRARY_ID)
  
  # 1) According to colData, filter samples ----
  
  ncol(count <- count[names(count) %in% names(x)])
  
  # 2) Filter data by removing low-abundance genes ----
  
  keep <- rowSums(edgeR::cpm(count) > 1) >= 2
  
  nrow(count <- count[keep,])
  
  count <- round(count)
  
  cat("Data to process: ", dim(count), "\n")
  
  # 3) Format metadata and count matrix ----
  
  colData <- colData %>% arrange(match(LIBRARY_ID, names(count)))
  
  sum(names(count) == colData$LIBRARY_ID) # sanity check
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  # using relevel, just specifying the reference level:
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = "Control"))
  
  library(DESeq2)
  library(tidyverse)
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = count,
    colData = colData,
    design = ~ Design )
  
  # 4) estimateSizeF for DESeq  ----
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  return(dds)
}

dds_B <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRASTE_B") # 110533 transcripts and 56 samples
dds_B <- DESeq(dds_B)
write_rds(dds_B, file = paste0(path, '/CONTRAST_B_DDS.rds'))

dds_C <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRASTE_C") # 111,174 transcripts and 59 samples
dds_C <- DESeq(dds_C)
write_rds(dds_C, file = paste0(path, '/CONTRAST_C_DDS.rds'))

dds_D <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRASTE_D") # 111,174 transcripts and 59 samples
dds_D <- DESeq(dds_D)
write_rds(dds_D, file = paste0(path, '/CONTRAST_D_DDS.rds'))

# 8) loop the contrasts results ----

# OPEN SCRIPT BIND_DATA_SOURCE_TO_DE.R
