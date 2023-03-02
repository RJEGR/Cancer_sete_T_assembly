

# Note: After estimateSizeFactors(ddsFullCountTable) 
# Evaluate the size factor, coef ad dispPriorVar from different source of design as follow: 
# sizeFactors(dds), head(coef(dds)), attr(dispersionFunction(dds), "dispPriorVar")
rm(list = ls()) # Limpiar la memoria de la sesion de R


estimateSizeFactorsFromSamples <- function(colData, count, f_col = ...) {
  
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
  
  df <- as.data.frame(sizeFactors(dds)) %>% as_tibble(rownames  = 'LIBRARY_ID')
  
  names(df)[2] <- 'sizeFactors'
  
  df <- df %>% mutate(col = f_col)
  
  return(df)
}


path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

mtd_f <- list.files(path, pattern = 'metadata.tsv',  full.names = TRUE)


mtd <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

dim(raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))


# colData %>% count(CONTRASTE_A, sort = T) # %>% view()
# colData %>% count(CONTRASTE_B, sort = T) # %>% view()
# colData %>% count(CONTRASTE_C, sort = T) 
# colData %>% count(CONTRASTE_D, sort = T) 


A <- estimateSizeFactorsFromSamples(mtd, raw_count, f_col = "CONTRASTE_A")
B <- estimateSizeFactorsFromSamples(mtd, raw_count, f_col = "CONTRASTE_B")
C <- estimateSizeFactorsFromSamples(mtd, raw_count, f_col = "CONTRASTE_C")
D <- estimateSizeFactorsFromSamples(mtd, raw_count, f_col = "CONTRASTE_D")

rbind(A,B,C,D) %>% ggplot(aes(x = LIBRARY_ID, y = sizeFactors)) + 
  facet_grid( col ~ .) + geom_col() +
  theme_bw(base_family = "GillSans", base_size = 12)

rbind(A,B,C,D) %>% pivot_wider(names_from = col, values_from = sizeFactors) %>% view()
