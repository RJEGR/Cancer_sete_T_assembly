rm(list = ls()) # Limpiar la memoria de la sesion de R

options(stringsAsFactors = FALSE) # 

library(tidyverse)
library(patchwork)

# Define Color Palette: ----

colNames <- c("Both Surviving",	"Forward Only",	"Reverse Only",	"Dropped")

# RColorBrewer::display.brewer.all() # Revisamos todas las paletas de colores que tenemos

getPalette <- RColorBrewer::brewer.pal(length(colNames), 'Set1')

Cvalues <- structure(getPalette, names = colNames)

path <- '~/Documents/GitHub/Cancer_sete_T_assembly/'

#Theme
my_theme <-
  ggplot2::theme_bw(base_size = 10, base_family = "GillSans") +
  ggplot2::theme(
    legend.position = 'top',
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank())

# 

file <- '.csv$'

file <- list.files(path, pattern = file, full.names = TRUE)

df <- lapply(file, read.table) # leer todas las tablas

df <- do.call(rbind, df)

colnames(df) <- c('id', 'Input Read Pairs',colNames)

head(df)

# basic plot ----

barplot(df$`Input Read Pairs`)

barplot(sort(df$`Input Read Pairs`))

# Algo publicable

df %>%
  pivot_longer(cols = all_of(colNames)) %>%
  filter(!name %in% "Input Read Pairs") %>%
  ggplot(aes(x = id, y = value, fill = name)) +
  geom_col() +
  labs(y = 'M reads', x = 'Sample') +
  coord_flip() -> p1

p1

# Intermediate ----

df %>%
  pivot_longer(cols = all_of(colNames)) %>%
  filter(!name %in% "Input Read Pairs") -> df_longer

df_longer %>%
  mutate(name = factor(name, levels = colNames)) %>%
  group_by(name) %>%
  mutate(id = forcats::fct_reorder(id, `Input Read Pairs`)) %>%
  ggplot(aes(x = id, y = value, fill = name)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(y = 'M reads', x = 'Sample') +
  coord_flip() +
  scale_fill_manual('', values = Cvalues) -> p2

p2

# Proficient

caption <- 'ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36 LEADING:5 TRAILING:5'

p2 +
  labs(caption = caption) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.border = element_blank()) -> p3

p3

# Percentage

df %>%
  mutate(across(c(3:6),
      .fns = ~./`Input Read Pairs`)) %>% 
  pivot_longer(cols = all_of(colNames)) %>%
  mutate(name = factor(name, levels = colNames)) %>%
  # arrange(desc(value)) %>%
  group_by(name) %>%
  mutate(id = forcats::fct_reorder(id, `Input Read Pairs`)) %>%
  ggplot(aes(x = id, y = value*100, fill = name)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(y = '%', x = 'Sample', 
    subtitle = caption) +
  coord_flip() +
  scale_fill_manual('', values = Cvalues) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.border = element_blank()) -> p4

p4

ggsave(p4, filename = 'trimmomatic_pct.png', path = path, 
  width = 7, height = 10)

library(patchwork)

p1+p2+p3+p4

# BUSCO DATAVZ ----

rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

path <- '~/Documents/DOCTORADO/summaries_busco/summaries_snp_trans_prev/full_tables/'

pattern_f <- 'tsv'


file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

read_tsv <- function(f) {
  
  g <- gsub(".tsv", "",basename(f))
  g <- gsub("full_table_transcripts_*", "", g)
  g <- gsub("_odb9", "", g)
  g <- str_to_title(g)
  
  df <- read.delim(f, comment.char = "#", sep = "\t", header = F)
  df %>% as_tibble(df) %>% mutate(g = g)
  #return(df)
}

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

names <- c("Complete", "Duplicated", "Fragmented", "Missing")

col <- structure(col, names = rev(names))

df <- lapply(file, read_tsv)

head(df <- do.call(rbind, df))

names(df) <- c("Busco_id",	"Status",	"Query_seq_id",	"Score",	"Length", "Db")

df %>% 
  group_by(Status, Db) %>% 
  tally() %>%
  group_by(Db) %>% mutate(pct = n / sum(n)) %>%
  mutate(Status = factor(Status, levels = rev(names))) %>%
  ggplot(aes(x = Db, y = pct, fill = Status)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "% BUSCOs", x = "") +
  coord_flip() +
  scale_fill_manual("", values = col) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = 'top')


df %>% 
  drop_na(Score) %>%
  ggplot(aes(Score, Length)) +
  geom_point(alpha = 0.5) +
  facet_grid(Status ~ Db, scales = "free") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = "BUSCO Score")


df %>% 
  drop_na(Score) %>%
  ggplot(aes(Status, Score)) +
  geom_boxplot() +
  facet_grid( ~ Db, scales = "free_y") +
  theme_bw(base_size = 14, base_family = "GillSans")

theme_set(my_theme)

df %>%
  group_by(Db) %>% 
  tally() %>%
  ggplot(aes(Db, n)) +
  geom_col() +
  labs(y = "Data set Size", x = "") +
  coord_flip() +
  geom_text(aes(label = n), size = 3, hjust = -0.05, family = "GillSans") 

df %>% 
  group_by(Status, Db) %>% 
  tally() %>%
  group_by(Db) %>% mutate(pct = n / sum(n)) %>%
  mutate(Status = factor(Status, levels = rev(names))) %>%
  ggplot(aes(x = Db, y = pct, fill = Status, group = Status)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "% BUSCOs", x = "") +
  coord_flip() +
  scale_fill_manual("", values = col) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = 'top')


# add genome index

path <- '~/Documents/DOCTORADO/summaries_busco/summaries_genome/full_tables/'
file <- list.files(path, pattern = pattern_f,  full.names = TRUE)
df2 <- lapply(file, read_tsv)
head(df2 <- do.call(rbind, df2))

names(df2) <- c("Busco_id",	"Status",	"Query_seq_id",	"Score",	"Length", "Db")


df2 %>% mutate(Index = "genome") -> df2

df %>% mutate(Index = "g_snp_trans") -> df

tbl <- rbind(df, df2)

tbl %>% 
  drop_na(Score) %>%
  ggplot(aes(Score, fill = Index)) +
  geom_histogram() +
  facet_grid(Status ~ Db, scales = "free") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = "BUSCO Score")

tbl %>% 
  drop_na(Length) %>%
  ggplot(aes(Length, fill = Index)) +
  geom_histogram() +
  facet_grid(Status ~ Db, scales = "free") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = "Sequence Length")

tbl %>% 
  drop_na(Score) %>%
  group_by(Status, Db, Index) %>% 
  tally() %>%
  filter(!Db %in% 'Bacteria') %>%
  # group_by(Db, Index) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(Status, n, fill = Index)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  facet_wrap( ~ Db, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  labs(y = "", x = "")

# Quantification ----

rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

pattern_f <- 'counts.matrix'

file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

read_tsv <- function(f) {
  g <- basename(f)
  g <- sapply(strsplit(g, "\\."), `[`, 2)
  g <- str_to_title(g)
  
  df <- read.delim(f, sep = "\t", header = T, row.names = 1)
  df %>% as_tibble(df) %>% mutate(g = g)
  #return(df)
}

# mtd <- read.csv(file = paste0(path, 'Mapping_file_FINAL.csv'))

df <- lapply(file, read_tsv)

head(df <- do.call(rbind, df))

df %>% distinct(g)


count_raw <- df %>% filter(g %in% 'Isoform') %>% select(-g)

# count <- edgeR::cpm(count) %>% as.data.frame()

nrow(count_raw)

apply(count_raw, 2, function(x) sum(x > 0)) -> Total_genes

nrow(count <- count_raw[rowSums(edgeR::cpm(count_raw)) > 1,])

# How singletones are per sample?

apply(count, 2, function(x) sum(x > 0)) -> filtered_genes

cbind(as_tibble(Total_genes, rownames = 'name'), as_tibble(filtered_genes)) -> n_genes

names(n_genes) <- c('name','Raw', 'Filt')

n_genes %>% mutate(pct = Raw - Filt) -> n_genes

n_genes %>%
  arrange(desc(pct)) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  ggplot() + geom_col(aes(x = name, y = pct), fill = 'black') +
  labs(x = '', y = 'Removed genes') +
  coord_flip() +
  theme_classic() -> p1

# Is a linearity between singletones and abundances?
# Raw barplot ----

n_genes %>% arrange(desc(pct)) %>% pull(name) -> idLev

count_raw %>% 
  pivot_longer(cols = names(count_raw), 
    values_to = 'Raw') %>%
  filter(Raw > 0) %>%
  group_by(name) %>%
  summarise(Raw = sum(Raw)) -> Total_size

# vs filtered

count %>% 
  pivot_longer(cols = names(count), 
    values_to = 'Filt') %>% 
  filter(Filt > 0) %>%
  group_by(name) %>%
  summarise(Filt = sum(Filt)) -> Filt_size

Total_size %>% 
  left_join(Filt_size) %>%
  arrange(desc(Raw)) %>%
  mutate(pct = Raw - Filt) -> n_reads


n_reads %>%
  mutate(name = factor(name, levels = idLev)) %>%
  ggplot() +
  geom_col(aes(x = name, y = pct)) +
  labs(x = '', y = 'Removed reads') +
  coord_flip() +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) -> p2

# p1 + p2

n_reads %>% arrange(match(idLev,name)) %>% pull(pct) -> x
n_genes %>% arrange(match(idLev,name)) %>% pull(pct) -> y

ggdf <- data.frame(idLev, x, y)

# Linearity between genes and reads

ggplot(ggdf, aes(x, y)) + 
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = TRUE, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  geom_point() +
  # geom_text(aes(label = idLev)) +
  labs(x = 'Reads', y = 'Genes') +
  theme_bw() +
  theme(panel.border = element_blank()) -> p3

p1 + p2 + p3

#

# PCA ----

PCA <- prcomp(t(log2(count+1)), scale. = FALSE) 

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>% dist(method = "euclidean") %>% hclust() %>% cutree(., 3) %>% as_tibble(rownames = 'id') %>%
  mutate(cluster = paste0('C', value)) %>% select(-value) -> hclust_res

PCAdf %>%
  mutate(id = rownames(.)) %>%
  left_join(hclust_res) %>%
  mutate(g = substr(id, 1,1)) %>%
  # left_join(mtd) %>%
  ggplot(., aes(PC1, PC2)) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(cluster)), fill = 'grey') +
  # geom_point(size = 5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_text(aes(label = id, color = g), alpha = 1) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  coord_fixed(ratio = sd_ratio) +
  scale_color_brewer(palette = "Set1")
  # facet_grid(~ Sample_Type)


# Prevalence of features ----

prevelancedf = apply(count, 1, function(x) sum(x > 0))

mean_se = apply(count, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(count),
  mean_se) %>% 
  as_tibble(rownames = "id") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf %>% 
  arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence))) %>%
  mutate(TotalAbundance = log2(TotalAbundance+1)) %>%
  ggplot(aes(TotalAbundance)) + geom_histogram() + 
  facet_wrap(~ Prevalence, scales = 'free_y') -> p1

dat_text <- prevelancedf %>% group_by(Prevalence) %>% tally() %>% 
  mutate(cumsum = cumsum(n)) %>% arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))

p1 + geom_text(
  data    = dat_text, family = "GillSans",
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " genes")),
  hjust   = -1,
  vjust   = -2
) + theme_classic(base_size = 7, base_family = "GillSans") +
  labs(x = expression(~Log[2]~('TotalAbundance'~+1)), y = "")

# 

# boxplot


count %>% 
  pivot_longer(cols = all_of(cols), names_to = 'id', values_to = 'count') %>% 
  filter(count > 1) -> count_longer


count_longer %>% mutate(count = log2(count+1)) %>% group_by(id) %>% summarise(q = quantile(count)) -> long_summ

long_summ %>%
  ungroup() %>%
  left_join(hclust_res) %>%
  arrange(match(cluster, cmatch)) %>%
  mutate(id = factor(id, levels = unique(id))) %>%
  ggplot(aes(x = id, y = q, fill = cluster, color = cluster)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), outlier.alpha = 0.5) +
  labs(y = expression(log[2]), x = '') +
  ylim(1, 5) + theme(legend.position = 'none') -> ptop

library(patchwork)

ptop / pbottom


# Normalized vs raw data count ----
# Rarefaction will be useful in evaluations of gene discovery using next-generation sequencing technologies as an analytical approach from theoretical ecology (rarefaction) to evaluate depth of sequencing coverage relative to gene discovery ... Rarefaction suggests that normalization has little influence on the efficiency of gene discovery, at least when working with thousands of reads from a single tissue type. (Hale, M. C., McCormick, et al 2009)


rs <- rowSums(round(count))

quantile(rs)

library(vegan)

Srar <- rarefy(round(count), min(rs))

head(Srar)

# ggrare
count %>% arrange(desc(rowSums(.))) %>% 
  # slice_head(n = 1000) -> x
  slice_sample(n = 1000, replace = TRUE) -> x

x <- t(round(x))

x <- as.data.frame(x)

tot <- rowSums(x)

S <- rowSums(x > 0)

nr <- nrow(x)

# i <- 1
step = 1000 # number: increment of the sequence.

se = TRUE

rarefun <- function(i) {
  cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
  n <- seq(1, tot[i], by = step)
  if (n[length(n)] != tot[i]) {
    n <- c(n, tot[i])
  }
  y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
  if (nrow(y) != 1) {
    rownames(y) <- c(".S", ".se")
    return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
  } else {
    return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
  }
}

out <- lapply(seq_len(nr), rarefun)

out <- do.call(rbind, out)

out <- out %>% mutate(g = substr(Sample, 1, 1))

labels <- data.frame(x = tot, y = S, Sample = rownames(x)) %>% 
  mutate(g = substr(Sample, 1, 1))

ggplot(data = out, aes(x = Size, y = `.S`, 
  group = Sample, color = g)) +
  facet_wrap(~g) +
  geom_line() + 
  geom_text(data = labels, ggplot2::aes(x,  y, label = Sample, color = g), 
    size = 4, 
    hjust = 0) +
  labs(x = "Sequence Sample Size", y = "# Genes") +
  geom_ribbon(aes(ymin = .S - .se, 
    ymax = .S + .se, color = NULL, fill = g), alpha = 0.2) -> p 

p + theme_classic() +
  theme(panel.border = element_blank())


# DESEQ2 ----

hclust_res %>% mutate(g = substr(id, 1,1)) %>%
  select(-id) %>% 
  mutate(g = factor(g, levels = unique(g))) %>%
  as.data.frame() -> colData

rownames(colData) <- hclust_res$id  

library(DESeq2)

count <- round(count)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ g ) # if not rep use design = ~ 1

# dds <- estimateSizeFactors(ddsFullCountTable)

dds <- DESeq(ddsFullCountTable)

contrast <- levels(colData$g)

get_res <- function(dds, contrast) {
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast)
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,colData(dds)$conditions == sA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,colData(dds)$conditions == sB])
  
  res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    as_tibble(rownames = "ids") %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("ids|sample|pvalue|padj")),
      round ,digits = 2)
}

res <- get_res(dds, contrast)

res <- prep_DE_data(res, padj_in = , logfc_in = 2) %>%
  drop_na(cc)

# Test run_DESEQ2(count, g)


logfc_in <- 2
padj_in <- 0.05 

res %>%
  ggplot(aes(y = -log10(pvalue), x = logFC)) +
  geom_point(aes(color = cc))

