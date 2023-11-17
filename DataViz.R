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
  width = 7.2, height = 10)

library(patchwork)

# p1+p2+p3+p4

# BUSCO DATAVZ ----

rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/summaries_busco/summaries_snp_trans_prev/full_tables/'

pattern_f <- 'tsv'


file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

read_tsv <- function(f) {
  
  g <- gsub(".tsv", "",basename(f))
  g <- gsub("full_table_transcripts_*", "", g)
  g <- gsub("_odb9", "", g)
  g <- stringr::str_to_title(g)
  
  df <- read.delim(f, comment.char = "#", sep = "\t", header = F)
  df %>% as_tibble() %>% mutate(g = g)
  #return(df)
}

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

names <- c("Complete", "Duplicated", "Fragmented", "Missing")

col <- structure(col, names = rev(names))

df <- lapply(file, read_tsv)

head(df <- do.call(rbind, df))

names(df) <- c("Busco_id",	"Status",	"Query_seq_id",	"Score",	"Length", "Db")

df %>% 
  filter(Db != "Bacteria") %>%
  group_by(Status, Db) %>% 
  tally() %>%
  group_by(Db) %>% mutate(pct = n / sum(n)) %>%
  mutate(Status = factor(Status, levels = rev(names))) %>% # view()
  ggplot(aes(x = Db, y = pct, fill = Status)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "% BUSCOs", x = "", caption = 'Human cancer reference-guide transcriptome assembly') +
  coord_flip() +
  scale_fill_manual("", values = col) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = 'top') -> ps

out_path <- '~/Documents/DOCTORADO/human_cancer_dataset/'

ggsave(ps, filename = "BUSCO.png", path = out_path, height = 3, width = 5)


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

# theme_set(my_theme)

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

path <- '~/Documents/DOCTORADO/human_cancer_dataset/summaries_busco/summaries_genome/full_tables/'

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
  geom_histogram(alpha = 0.5) +
  facet_wrap(Status ~ Db, scales = "free") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = "BUSCO Score")

tbl %>% 
  drop_na(Length) %>%
  ggplot(aes(Length, fill = Index)) +
  geom_histogram() +
  facet_wrap(Status ~ Db, scales = "free") +
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
  labs(y = "Transcripts", x = "") -> ps


ggsave(ps, filename = "BUSCO_index_contrast.png", path = out_path, height = 3)


# Quantification ----


path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp'

pattern_f <- 'counts.matrix$'

file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

read_tsv <- function(f) {
  g <- basename(f)
  g <- sapply(strsplit(g, "\\."), `[`, 2)
  g <- str_to_title(g)
  
  df <- read.delim(f, sep = "\t", header = T, row.names = 1)
  # df %>% as_tibble(df, rownames = 'id') %>% mutate(g = g)
  return(df)
}

mtd <- readr::read_tsv(file = paste0(path, '/metadata.tsv'))

rownames(mtd) <- mtd$LIBRARY_ID


df <- lapply(file, read_tsv)

dim(df <- do.call(rbind, df))

# df %>% distinct(g)

# Select isoforms count

count_raw <- df # %>% filter(g %in% 'Isoform') %>% select(-g)

nrow(count_raw)

apply(count_raw, 2, function(x) sum(x > 0)) -> Total_genes

# Filter data by removing low-abundance genes

keep <- rowSums(edgeR::cpm(count_raw) > 1) >= 2


nrow(count <- count_raw[keep,])

# How singletones are per sample? ----

apply(count, 2, function(x) sum(x > 0)) -> filtered_genes

cbind(as_tibble(Total_genes, rownames = 'name'), as_tibble(filtered_genes)) -> n_genes

names(n_genes) <- c('LIBRARY_ID','Raw', 'Filt')

n_genes %>% mutate(pct_genes_retained = Filt/Raw) -> n_genes

# Total genes

n_genes %>%
  # left_join(mtd) %>%
  mutate(g = substr(LIBRARY_ID, 1,1)) %>%
  arrange(desc(Raw)) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = unique(LIBRARY_ID))) %>%
  ggplot() + geom_col(aes(x = LIBRARY_ID, y = Raw),  fill = 'black') + # 
  labs(x = '', y = 'Total Transcripts') +
  scale_y_continuous(labels = scales::comma) +
  coord_flip() +
  theme_classic(base_family = "GillSans") -> p

n_genes %>%
  mutate(g = substr(LIBRARY_ID, 1,1)) %>%
  arrange(desc(Raw)) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = unique(LIBRARY_ID))) %>%
  ggplot() + geom_col(aes(x = LIBRARY_ID, y = pct_genes_retained)) +
  labs(x = '', y = 'Retained (%)') +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  coord_flip() +
  theme_classic(base_family = "GillSans") +
  scale_fill_brewer('', palette = "Set1") +
  theme(
    legend.position = 'none',
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) -> p1

# Is a linearity between singletones and abundances?
# Raw barplot ----

n_genes %>% arrange(desc(Raw)) %>% pull(LIBRARY_ID) -> idLev

count_raw %>% 
  pivot_longer(cols = names(count_raw), 
    values_to = 'Raw', names_to = 'LIBRARY_ID') %>%
  filter(Raw > 0) %>%
  group_by(LIBRARY_ID) %>%
  summarise(Raw = sum(Raw)) -> Total_size

# vs filtered

count %>% 
  pivot_longer(cols = names(count), 
    values_to = 'Filt', names_to = 'LIBRARY_ID') %>% 
  filter(Filt > 0) %>%
  group_by(LIBRARY_ID) %>%
  summarise(Filt = sum(Filt)) -> Filt_size

Total_size %>% 
  left_join(Filt_size) %>%
  arrange(desc(Raw)) %>%
  mutate(pct_reads_retained = Filt/Raw) -> n_reads

n_reads %>%
  mutate(g = substr(LIBRARY_ID, 1,1)) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = idLev)) %>%
  ggplot() +
  geom_col(aes(x = LIBRARY_ID, y = pct_reads_retained)) +
  labs(x = '', y = 'Retained reads (%)') +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  coord_flip() +
  scale_fill_brewer('', palette = "Set1") +
  theme_classic(base_family = "GillSans") +
  theme(
    legend.position = 'top',
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) -> p2

# p+ p1 + p2

n_reads %>% left_join(n_genes, by = "LIBRARY_ID") %>% pull(pct_reads_retained) -> x
n_reads %>% left_join(n_genes, by = "LIBRARY_ID") %>% pull(pct_genes_retained) -> y

# n_reads %>% arrange(match(idLev,LIBRARY_ID)) %>% pull(pct_reads_retained) -> x
# n_genes %>% arrange(match(idLev,LIBRARY_ID)) %>% pull(pct_genes_retained) -> y

# stats?
n_reads %>% mutate(p = (pct_reads_retained * 100) / Raw)
n_genes %>% mutate(p = (pct_genes_retained * 100) / Raw)

ggdf <- data.frame(idLev, x, y)

# Linearity between genes and reads

ggplot(ggdf, aes(x, y)) + 
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = TRUE, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  geom_point() +
  # geom_text(aes(label = idLev)) +
  labs(x = 'Reads', y = 'Transcripts') +
  theme_bw(base_family = "GillSans") +
  theme(panel.border = element_blank()) -> p3

p+ p1 + p2 + p3 


# boxplot of reads ----

n_genes %>% 
  left_join(mtd) %>%
  group_by(CONTRASTE_A) %>%
  arrange(desc(Raw)) %>%
  mutate(g = substr(LIBRARY_ID, 1,1)) %>%
  # arrange(desc(pct_genes_retained)) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = unique(LIBRARY_ID))) %>%
  ggplot() + geom_col(aes(x = LIBRARY_ID, y = Raw, fill = CONTRASTE_A)) + # fill = 'black'
  labs(x = '', y = 'Total Transcripts') +
  scale_y_continuous(labels = scales::comma) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position = 'bottom',
    axis.text.x = element_text(angle = 90,
    hjust = 1, vjust = 1, size = 10)) -> pbottom

# ?gtsummary::tbl_survfit() to retrive probs

qprobs <- function(x) { 
  x <- x[x > 1]
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}

apply(log2(count+1), 2, qprobs) %>% 
  t() %>%
  as_tibble(rownames = 'LIBRARY_ID') -> probs_df




probs_df %>%
  left_join(mtd) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID,  levels = levels(pbottom$data$LIBRARY_ID))) %>%
  ggplot(., 
    aes(x = LIBRARY_ID, ymin = `5%`, lower = `25%`,
      middle = `50%`, upper = `75%`, ymax = `95%`)) +
  geom_errorbar(width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(aes(fill = CONTRASTE_A), width = 0.5, stat = 'identity', position = position_dodge(0.6)) +
  labs(y = expression(log[2]~ 'Reads'), x = '') +
  # ylim(1, 5) + 
  theme_classic(base_family = "GillSans") +
  theme(
    legend.position = 'none',
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()) -> ptop


library(patchwork)

ps <- ptop / pbottom + patchwork::plot_layout(heights = c(1,1.5))

out_path <- '~/Documents/DOCTORADO/human_cancer_dataset/'

ggsave(ps, filename = "filtered_transcripts_and_reads_plots.png", 
  path = out_path, width = 10, height = 5.5)


# CLUSTERING SAMPLES -----
# Using Bayesian Information Criterion for expectation-maximization, to select kmeans clusters instead of elbow method
# https://towardsdatascience.com/are-you-still-using-the-elbow-method-5d271b3063bd

data = log2(count+1)

library(mclust)

# d_clust <- Mclust(as.matrix(data), G=1:15, modelNames = mclust.options("emModelNames"))  # TAKE LONG TIME (~ 5 hours)

d_clust$BIC

k <- d_clust$G # 14
names(k) <- d_clust$modelName

# Top 3 models based on the BIC criterion: 
#   VEV,14    VEV,15    VEV,13 
# -20957084 -20965453 -20972947

# POR LO TANTO, EL K PREDICHO MEDIANTE EL METODO BAYESIANO FUE 14, METODO VEV
# AUNQUE POCO RARO PUES DE LOS 15 VALORES EN G RESULTO 14 EL NUMERO IDEAL. HAY QUE REVISAR!
plot(d_clust) # CHOOSE OPTION 1

# TEST NbClust method -----


# PCA ----
ncol(data <- log2(count+1))

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

d_clust <- Mclust(as.matrix(PCAdf), G=1:15, modelNames = mclust.options("emModelNames"))

PCAdf %>% 
  dist(method = "euclidean") %>% 
  hclust() %>% 
  cutree(., 14) %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>% 
  mutate(cluster = paste0('C', value)) %>% 
  dplyr::select(-value) -> hclust_res

# getPalette <- RColorBrewer::brewer.pal(n, 'Paired')

# wesanderson::wes_palette("Cavalcanti1")

# getPalette <- wesanderson::wes_palette("BottleRocket2", n, type = "continuous")

palette <- c(
  "#8e0152",
  "#c51b7d",
  "#de77ae",
  "#f1b6da",
  "#fde0ef",
  "#f7f7f7",
  "#e7d4e8",
  "#c2a5cf",
  "#9970ab",
  "#e6f5d0",
  "#a6dba0",
  "#5aae61",
  "#1b7837",
  "#00441b")

color_vector <- unique(mtd$CONTRASTE_D)
n <- length(color_vector)
getPalette <- sample(palette, n)

axis_col <- structure(getPalette, names = color_vector)

c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  
axis_col[names(axis_col) %in% "SIN_CANCER"] <- "#313695"

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  # mutate(g = substr(sample_id, 1,1)) %>%
  left_join(hclust_res) %>%
  left_join(mtd) %>%
  ggplot(., aes(PC1, PC2)) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(cluster)),
    fill = 'grey', con.colour = 'grey') +
  geom_point(size = 3, alpha = 0.9, aes(color = CONTRASTE_D)) 
  # coord_fixed(ratio = sd_ratio) 
  # scale_color_manual(values = axis_col)
  # geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  ggrepel::geom_text_repel(aes(label = LIBRARY_ID, color = CONTRASTE_B), 
    alpha = 1, size = 4, max.overlaps = 60) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_classic(base_family = "GillSans") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') -> pcaplot

ggsave(pcaplot, 
  filename = "PCA.png", path = path, 
  width = 10, height = 7)




# Correlation heatmap ----

sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

sample_cor %>% 
  as_tibble(rownames = 'sample_id') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(mtd) %>% 
  mutate(sample_id = factor(sample_id, levels = hc_order)) -> sample_cor_long

# sample_cor_long %>% distinct(sample_id, Diagnosis) %>% mutate(col = ifelse(g %in% 'C', 'red', 'blue')) -> coldf 

n <- length(unique(sample_cor_long$Diagnosis))

getPalette <- RColorBrewer::brewer.pal(n, 'Paired')

axis_col <- structure(getPalette, names = unique(sample_cor_long$Diagnosis))
axis_col <- getPalette[match(mtd$Diagnosis, names(axis_col))]
axis_col <- structure(axis_col, names = unique(sample_cor_long$sample_id))


# structure(coldf$col, names =  as.character(coldf$Sample)) -> axis_col

library(ggh4x)

sample_cor_long %>%
  ggplot(aes(x = sample_id, y = name, fill = cor)) + 
  # geom_tile(color = 'white', size = 0.2) +
  geom_raster() + 
  # geom_text(aes(label = cor), color = 'white') +
  scale_fill_viridis_c(name = "Pearson", direction = -1) +
  # scale_x_discrete(position = 'top') +
  labs(x = '', y = '') +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples) +
  theme_classic(base_size = 7, base_family = "GillSans") +
  theme(axis.text.x = element_text(angle = 90,
    hjust = 1, vjust = 1, size = 7, color = axis_col),
    axis.text.y = element_text(color = axis_col),
    axis.ticks.length = unit(5, "pt")) -> pheat

# pheat +  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') -> pheat

ggsave(pheat, 
  filename = "cor_data_matrix.png", path = path, 
  width = 6, height = 5)


# Prevalence of features ----
# TRANSCRIPT VIEW

apply(count_raw, 1, function(x) sum(x > 0)) %>% table()

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
  mutate(TotalAbundance = log2(TotalAbundance+1)) %>% # edgeR::cpm(count) 
  ggplot(aes(TotalAbundance)) + geom_histogram() + 
  facet_wrap(~ Prevalence, scales = 'free_y') -> p1

dat_text <- prevelancedf %>% group_by(Prevalence) %>% tally() %>% 
  mutate(cumsum = cumsum(n)) %>% arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))

p1 <- p1 + geom_text(
  data    = dat_text, family = "GillSans",
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " genes")),
  hjust   = -1, vjust   = -2, label.size = 0.2) + 
  theme_classic(base_size = 7, base_family = "GillSans") +
  labs(x = expression(~Log[2]~('TotalAbundance'~+1)), y = "")

#

ggplot(prevelancedf, aes(TotalAbundance, Prevalence/60)) +
  stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE) +
  scale_x_log10("Total Abundance (log10 scale)", labels = scales::comma) +  
  scale_y_continuous("Prevalence [Frac. Samples]", 
    labels = scales::percent_format(scale = 100)) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position="top") -> ps

ggsave(ps, filename = 'prevalence_and_abundant_density.png', path = path)
ggsave(p1, filename = 'prevalence_and_abundant_hist.png', path = path, width = 7,height = 5)

# summary
# apply(count_raw, 1, function(x) sum(x > 0)) %>% table() %>% data.frame()

prevelancedf %>% 
  count(Prevalence) %>% view()
  ggplot(aes(Prevalence, n)) + geom_col() +
  # ggplot(aes(Prevalence)) + 
  # geom_histogram(binwidth = 1) +
  theme_classic(base_family = "GillSans") + 
  scale_y_continuous("Number of transcripts", labels = scales::comma) -> ps

ggsave(ps, filename = 'prevalence_hist.png', path = path, width = 3, height = 3)
# 
# (omit) Normalized vs raw data count ----
# Rarefaction will be useful in evaluations of gene discovery using next-generation sequencing technologies as an analytical approach from theoretical ecology (rarefaction) to evaluate depth of sequencing coverage relative to gene discovery ... Rarefaction suggests that normalization has little influence on the efficiency of gene discovery, at least when working with thousands of reads from a single tissue type. (Hale, M. C., McCormick, et al 2009)

library(vegan)


# ggrare
count %>% arrange(desc(rowSums(.))) %>% 
  # slice_head(n = 1000) -> x
  slice_sample(n = 1000, replace = TRUE) -> m

x <- t(round(m))

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
  # facet_wrap(~g) +
  geom_line() + 
  geom_text(data = labels, ggplot2::aes(x,  y, label = Sample, color = g), 
    size = 4, 
    hjust = 0) +
  labs(x = "Sequence Sample Size", y = "# Genes") +
  geom_ribbon(aes(ymin = .S - .se, 
    ymax = .S + .se, color = NULL, fill = g), alpha = 0.2) -> p 

p + theme_classic() +
  theme(panel.border = element_blank())

# Normalized (hard to run, therefore use n samples)

x <-  edgeR::cpm(m[, 1:4]) %>% as.data.frame()

x <- t(round(x))

x <- as.data.frame(x)

tot <- rowSums(x)

S <- rowSums(x > 0)

nr <- nrow(x)

outN <- lapply(seq_len(nr), rarefun)

outN <- do.call(rbind, outN)

out %>% filter(Sample %in% unique(outN$Sample)) -> outR

outN %>% 
  as_tibble() %>%
  mutate(g = substr(Sample, 1, 1)) %>%
  mutate(g = paste0(g, '-N')) %>%
  rbind(., outR) -> outNp

unique(outNp$Sample)

ggplot(data = outNp, aes(x = Size, y = `.S`, 
  group = Sample, color = g)) +
  facet_wrap(~g, scales = 'free_x') +
  geom_line() + 
  geom_ribbon(aes(ymin = .S - .se, 
    ymax = .S + .se, color = NULL, fill = g), alpha = 0.2) +
  theme_classic() +
  labs(x = "Sequence Sample Size", y = "# Genes")

# DESEQ2 ----


colors_fc <- c("red2", 
  "#4169E1",
  "forestgreen", "grey30")

colData <- mtd %>% arrange(match(sample_id, names(count)))

colData <- mutate_if(colData, is.character, as.factor)

# using relevel, just specifying the reference level:

colData$g2 <- relevel(colData$g2, ref = "Control")

levels(colData$g2)

library(DESeq2)

count <- round(count)

table(colData$g2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ g1 ) # if not rep use design = ~ 1

# dds <- estimateSizeFactors(ddsFullCountTable)

dds <- DESeq(ddsFullCountTable)

# write_rds(dds, file = paste0(path, '/Cancer_vs_Control_dds.rds'))

# dds <- read_rds("~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/Cancer_vs_Control_dds.rds")

contrast <- levels(colData(dds)$g1)

get_res <- function(dds, contrast, alpha_cutoff = 0.1) {
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast, alpha = alpha_cutoff)
  
  
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

prep_DE_data <- function(res, alpha, lfcThreshold) {
  
  FC_cols <- c('logFC','log2FoldChange')
  pv_cols <- c('P.Value','pvalue', 'PValue')
  pvad_cols <- c('padj', 'FDR', 'adj.P.Val')
  
  sam <- c('sampleA',  'sampleB')
  samv <- c('baseMeanA', 'baseMeanB')
  
  rename_to <- c('logFC', 'pvalue', 'padj')
  # 
  # sigfc <- "<b>P</b>-value & Log<sub>2</sub> FC"
  # pv <- "<b>P</b>-value"
  # fc <- "Log<sub>2</sub> FC"
  # 
  # sigfc <- expression(p - value ~ and ~ log[2] ~ FC)
  # pv <- "p-value"
  # fc <- expression(Log[2] ~ FC)
  # 
  sigfc <- "p - value ~ and ~ log[2] ~ FC"
  pv <- "p-value"
  fc <- "Log[2] ~ FC"

  
  sample_NS <- function(res, n) {
    
    x <- filter(res, !(cc %in% c('NS', fc)))
    
    y <- sample_n(filter(res, cc == fc), n)
    
    z <- sample_n(filter(res, cc == 'NS'), n) 
    
    dat_sampled <- rbind(x,y,z)
    
    return(dat_sampled)
  }
  
  res %>%
    # drop_na() %>%
    select_at(vars(ids, contains(sam), contains(samv), contains(FC_cols), 
      contains(pv_cols),
      contains(pvad_cols))) %>%
    rename_at(vars(contains(FC_cols),
      contains(pv_cols),
      contains(pvad_cols)), ~ rename_to) %>%
    arrange(pvalue) -> res
  
  
  res$cc <- 'NS'
  res[which(abs(res$logFC) > lfcThreshold), 'cc'] <- fc
  res[which(abs(res$padj) <= alpha), 'cc'] <- pv
  res[which(res$padj <= alpha & abs(res$logFC) > lfcThreshold), 'cc'] <- sigfc
  
  up <- 'Up-regulated'
  down <- 'Down-regulated'
  
  res$lfcT <- 'Basal'
  res[which(res$padj < alpha & res$logFC > lfcThreshold), 'lfcT'] <- up
  res[which(res$padj < alpha & res$logFC < -lfcThreshold), 'lfcT'] <- down
  
  res %>%
    # sample_NS(.,1000) %>%
    mutate(cc = factor(cc, levels = c(sigfc, pv, fc, "NS"))) %>%
    mutate(lfcT = factor(lfcT, levels = c(up, down, 'Basal')))

}

res <- get_res(dds, contrast)

res.p <- prep_DE_data(res, alpha = 0.05, lfcThreshold = 2) %>%
  drop_na(cc)

# Test run_DESEQ2(count, g)


logfc_in <- 2

padj_in <- 0.05 

# volcano ----

res.p %>%
  mutate(pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = logFC, y = pvalue)) +
  geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = colors_fc) + 
  labs(x= expression(Log[2] ~ "Fold Change"), 
    y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top") -> pvol
  # geom_abline(slope = 0, intercept = -log10(padj_in), linetype="dashed", alpha=0.5) 
  # geom_vline(xintercept = 0, linetype="dashed", alpha=0.5)

ggsave(pvol, 
  filename = "volcano.png", path = path, 
  width = 8, height = 5)

# res %>% arrange(log2FoldChange)
# - logfc == Cancer samples
# + logfc == Control samples

res.p %>%
  mutate(g = ifelse(logFC > 0, 'Control', 'Cancer')) %>%
  group_by(cc,g ) %>% tally() %>%
  filter(cc != 'NS') %>%
  ggplot() +
  geom_bar(aes(x = g, y = n, fill = cc), 
    stat = 'identity', position = position_dodge2()) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = '# Transcripts', x = '') +
  theme(legend.position = 'top') -> pbar

# ggsave(pbar, 
#   filename = "diffExp_bar.png", path = path, 
#   width = 5, height = 5)

res.p %>%
  filter(cc != 'NS') %>%
  filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  ggplot() +
  geom_histogram(aes(padj, fill = cc)) +
  # facet_wrap( cc ~ ., scales = 'free_x') +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = 'none') +
  labs(y = '# Transcripts') -> padjp

pbar / padjp -> savep

ggsave(savep, 
  filename = "diffExp_bar.png", path = path, 
  width = 5, height = 5)

# and bar plot

# res.p %>%
#   filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
#   pivot_longer(cols = c('baseMeanA', 'baseMeanB')) %>%
#   group_by(name) %>%
#   # arrange(desc(logFC)) %>%
#   slice_head(n = 20) %>%
#   ggplot(aes(x = ids, y = value, fill = name)) +
#   geom_col(position = position_dodge2()) +
#   labs(y = 'Base Mean', x = 'Genes') +
#   theme_classic() +
#   theme(
#     legend.position = 'top',
#     panel.border = element_blank(),
#     axis.text.x = element_blank(),
#     # axis.ticks.x = element_blank(),
#     axis.line.x = element_blank()) +
#   facet_grid(~name)

# heatmap ----
# returno to heatmap to get positions


path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp'
file_name <- paste0(path, '/multiple_contrast_vs_control_res.rds')
res <- read_rds(file_name) # Para leer archivos tipo rds (Rdata) desde R. El formato rds es una version de archivos compridos de R

res %>% distinct(sampleA, sampleB)

res %>% arrange(log2FoldChange)

# res %>% sample_n(1000) %>%
#   ggplot() +
#   geom_point(aes(y = log2FoldChange, x = baseMean, color = -log10(padj)))


alpha = 0.05 

lfcThreshold = 2

# Queremos genes que fueron superiores al umbral lfcThreshold & < alpha == FC+sig

res %>% 
  mutate(cc = NA) %>% 
  mutate(cc = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, 'sigfc', cc)) %>%
  drop_na(cc) -> res.p

up <- 'Up-regulated'
down <- 'Down-regulated'

res.p %>%
  mutate(lfcT = 'Basal') %>%
  mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, up, lfcT)) %>%
  mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, down, lfcT)) -> res.p


res.p %>%
  # filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  pull(ids) -> diffexpList

keep <- rownames(count) %in% diffexpList

dim(count[keep,] -> DEcount)

sample_cor = cor(DEcount, method='pearson', use='pairwise.complete.obs')

# superheat::superheat()

dist.method <- 'euclidean'
linkage.method <- 'complete'

m <- log2(DEcount+1)

hc_samples <- hclust(dist(t(m), method = dist.method), 
  method = linkage.method)

# sample_dist = dist(t(DEcount), method='euclidean')
# sample_dist = dist(sample_cor, method='euclidean')

hc_sam_order <- hc_samples$labels[hc_samples$order]

hc_genes <- hclust(dist(m, method = dist.method), 
  method = linkage.method)

hc_genes_ord <- hc_genes$labels[hc_genes$order]


DEcount %>% 
  as_tibble(rownames = 'TRANSCRIPT') %>%
  pivot_longer(cols = colnames(DEcount), values_to = 'count', names_to = 'LIBRARY_ID') %>%
  # mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_sam_order)) %>%
  # mutate(id = factor(id, levels = hc_genes_ord)) %>%
  left_join(mtd) %>%
  filter(count > 0) -> countLong

countLong %>% 
  distinct(LIBRARY_ID, CONTRASTE_A) %>% 
  # mutate(col = ifelse(g1 %in% 'C', 'red', 'blue')) %>%
  arrange(match(LIBRARY_ID, hc_sam_order)) -> coldf 


n <- length(unique(coldf$CONTRASTE_A))

getPalette <- RColorBrewer::brewer.pal(n+1, 'Paired')

axis_col <- structure(getPalette, names = unique(coldf$CONTRASTE_A))
axis_col <- getPalette[match(mtd$CONTRASTE_A, names(axis_col))]
axis_col <- structure(axis_col, names = unique(coldf$LIBRARY_ID))



library(ggh4x)

countLong %>%
  # sample_n(1000) %>%
  ggplot(aes(x = LIBRARY_ID, y = TRANSCRIPT, fill = log2(count+1))) + 
  geom_raster() + 
  theme_classic(base_size = 12, base_family = "GillSans") +
  scale_fill_viridis_c(name = expression(Log[2]~'(count+1)')) +
  labs(x = '', y = '') +
  ggh4x::scale_y_dendrogram(hclust = hc_genes) +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  theme(axis.text.x = element_text(angle = 90, 
    hjust = 1, vjust = 1, size = 7, color = axis_col), # 
    axis.ticks.length = unit(10, "pt"),
    # legend.position = 'top',
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) -> pheat

# pheat +facet_grid(~ Diagnosis, scales = 'free_x', space = 'free_x')

pheat + guides(fill = guide_colorbar(barheight = unit(4, "in"),
  ticks.colour = "black",
  frame.colour = "black",
  label.theme = element_text(size = 12))) -> pheat

# pheat + facet_grid(~ g1, scales = 'free_x', space = 'free_x')

ggsave(pheat, 
  filename = "diffExp_hetmap.png", path = path, 
  width = 8, height = 8)

# Functional annotation ----
# test https://yulab-smu.top/biomedical-knowledge-mining-book/index.html

rm(list = ls())

library(tidyverse)

library(topGO)

runtopGO <- function(topGOdata, topNodes = 20, conservative = TRUE) {
  
  RFisher <- runTest(topGOdata, 
    algorithm = "classic", 
    statistic = "fisher")
  
  # To make this test conservative. Next we will test the enrichment using the Kolmogorov-Smirnov test. We will use the both the classic and the elim method.
  
  if(conservative) 
  {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks")
    
    RKS.elim <- runTest(topGOdata, algorithm = "elim", 
      statistic = "ks")
    
    
    
    
    allRes <- GenTable(topGOdata, 
      classicFisher = RFisher,
      classicKS = RKS, 
      elimKS = RKS.elim,
      orderBy = "elimKS", 
      ranksOf = "classicFisher", 
      topNodes = topNodes) 
  } else {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks")
    
    test.stat <- new("weightCount",
      testStatistic = GOFisherTest,
      name = "Fisher test", sigRatio = "ratio")
    
    weights <- getSigGroups(topGOdata, test.stat)
    
    allRes <- GenTable(topGOdata,
      classic = RFisher,
      KS = RKS,
      weight = weights,
      orderBy = "weight",
      ranksOf = "classic",
      topNodes = topNodes)
    
    # allRes <- GenTable(object = topGOdata, 
    #                    elimFisher = RFisher,
    #                    topNodes = topNodes)
  }
  
  return(allRes)
}

description <- "complete topGO enrichment using split_annot"

path <- "~/Documents/DOCTORADO/human_cancer_dataset/annot/"

# write_rds(res.p, file = paste0(path, 'Cancer_vs_Control_res_p.rds'))
res.p <- read_rds(paste0(path, 'Cancer_vs_Control_res_p.rds'))

go_file <- paste0(path, 'Trinotate_report.xls.gene_ontology')

MAP <- topGO::readMappings(go_file)

# geneNames <- sort(unique(as.character(names(MAP))))

# Load semantic similarity data

hsGO <- read_rds(paste0(path, '/hsGO_BP.rds'))

# Prepare semantic similarity data

# hsGO <- GOSemSim::godata('org.Hs.eg.db', ont="BP")
# hsGOMF <- GOSemSim::godata('org.Hs.eg.db', ont="MF")

# saveRDS(hsGO, file = paste0(path, '/hsGO_BP.rds'))
# saveRDS(hsGOMF, file = paste0(path, '/hsGO_MF.rds'))


# countLong %>% distinct(id) %>% pull(.) %>% as.character() -> query.genes


res.p %>% filter(cc ==  'p - value ~ and ~ log[2] ~ FC') -> res.p

res.p %>% filter(logFC > 0) -> query.res.g1
res.p %>% filter(logFC < 0) -> query.res.g2

query.p.1 <- query.res.g1 %>% pull(padj)
query.names.1 <- query.res.g1 %>% pull(ids)

query.p.2 <- query.res.g2 %>% pull(padj)
query.names.2 <- query.res.g2 %>% pull(ids)


str(query.names.1) # 3412 diffExp genes (1157 vs 2255 genes) and 
str(query.names.2) # 3412 diffExp genes (1157 vs 2255 genes) and 


topEnrichment <- function(query.p, query.names, MAP, cons = T, onto = "BP") {
  
  names(query.p) <- query.names
  
  # keep MAP of query genes
  keep <- names(MAP) %in% names(query.p) 
  
  # MAP <- MAP[keep]
  
  topGOdata <- new("topGOdata", 
    ontology = onto, 
    description = description,
    allGenes = query.p,
    # geneSel = function(x) { x == 1 },
    geneSel = function(x) x,
    annot = annFUN.gene2GO,
    # mapping = hsGO, # omit this flag
    gene2GO = MAP[keep])
  
  # run TopGO results 
  
  allGO <- usedGO(topGOdata)
  
  allRes <- runtopGO(topGOdata, topNodes = length(allGO), conservative = cons)
  
  # make p adjustable
  
  p.adj.ks <- p.adjust(allRes$classicKS , method="BH")
  
  allRes <- cbind(allRes, p.adj.ks)
  
  allRes$Term <- gsub(" [a-z]*\\.\\.\\.$", "", allRes$Term)
  allRes$Term <- gsub("\\.\\.\\.$", "", allRes$Term)
  
  return(allRes)

  
}

# it perform better than other ways

goEnrich_g1 <- topEnrichment(query.p.1, query.names.1, MAP)
goEnrich_g2 <- topEnrichment(query.p.2, query.names.2, MAP)

nrow(goEnrich_g1)
nrow(goEnrich_g2)

# sampleA == Control (C), sampleB == Cancer (P)
# sampleA == Control (C), sampleB == Cancer (P)
# baseMeanA == +FC, baseMeanB == -FC


write.table(goEnrich_g1, file = paste0(path, 'go_enrichment_control_vs_cancer_control.csv'))
write.table(goEnrich_g2, file = paste0(path, 'go_enrichment_control_vs_cancer_cancer.csv'))

quantile(as.numeric(goEnrich_g1$classicKS))
quantile(goEnrich_g2$p.adj.ks)

nrow(topGO1.p <- goEnrich_g1[which(goEnrich_g1$classicKS <= 0.05),])
nrow(topGO2.p <- goEnrich_g2[which(goEnrich_g2$p.adj.ks <= 0.05),])

#

# length(MAP1 <- MAP[names(MAP) %in% query.names.1])
# length(MAP2 <- MAP[names(MAP) %in% query.names.2])
# 
# STRG2GO <- data.frame(ID = rep(names(MAP1),
#   sapply(MAP1, length)),
#   GO.ID = unlist(MAP1))


# After enrichment analysis, introduce reduction of terms using 
# Take a while 

# test pipeline from termAnalysis using the moduleTraitRelationshipHeatmap.R code

library(GOSemSim)
library(rrvgo)

str(goid1 <- sort(unique(topGO1.p$GO.ID))) # using only significant genes
str(goid2 <- sort(unique(topGO2.p$GO.ID))) # using only significant genes


scores1 <- -log(as.numeric(topGO1.p$classicKS))
scores2 <- -log(as.numeric(topGO2.p$classicKS))

names(scores1) <- topGO1.p$GO.ID
names(scores2) <- topGO2.p$GO.ID

sort(names(scores2))

# SimMatrix <- calculateSimMatrix(goid, orgdb = 'org.Hs.eg.db', ont="BP", method = 'Wang')

SimMatrix1 <- GOSemSim::termSim(goid1, goid1, semData = hsGO,  method = "Wang")
SimMatrix2 <- GOSemSim::termSim(goid2, goid2, semData = hsGO,  method = "Wang")

reducedTerms1 <- reduceSimMatrix(SimMatrix, scores1, threshold = 0.9, orgdb = 'org.Hs.eg.db')
reducedTerms2 <- reduceSimMatrix(SimMatrix2, scores2, threshold = 0.9, orgdb = 'org.Hs.eg.db')

scatterPlot <- function (simMatrix, reducedTerms, size = "score", addLabel = TRUE, 
  labelSize = 3) 
{
  if (!all(sapply(c("ggplot2", "ggrepel"), requireNamespace, 
    quietly = TRUE))) {
    stop("Packages ggplot2, ggrepel and/or its dependencies not available. ", 
      "Consider installing them before using this function.", 
      call. = FALSE)
  }
  x <- cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE, 
    k = 2)
  df <- cbind(as.data.frame(x$points), reducedTerms[match(rownames(x$points), 
    reducedTerms$go), c("term", "parent", "parentTerm", 
      "size")])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = V1, y = V2, color = parentTerm)) + 
    ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) + 
    ggplot2::scale_color_discrete(guide = 'none') + 
    ggplot2::scale_size_continuous(guide = 'none', range = c(0, 25)) + 
    ggplot2::scale_x_continuous(name = "") + 
    ggplot2::scale_y_continuous(name = "") + 
    ggplot2::theme_minimal(base_size = 12, base_family = "GillSans") + 
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
      axis.text.y = ggplot2::element_blank()) +
    labs(x = 'Dimension 1', y = 'Dimension 2')
  if (addLabel) {
    p + ggrepel::geom_label_repel(aes(label = parentTerm), 
      data = subset(df, parent == rownames(df)), box.padding = grid::unit(1, 
        "lines"), size = labelSize, family = "GillSans")
  }
  else {
    p
  }
}


scatterPlot(SimMatrix1, reducedTerms1) -> sim1plot
scatterPlot(SimMatrix2, reducedTerms2) -> sim2plot

library(patchwork)

sim1plot + sim2plot -> psave

ggsave(psave, filename = 'reducedTerms_semantic.png', 
  path = path, 
  width = 9.5, height = 7.5, dpi = 250)


# as barplots
# g1 == Control (C) and g2 == Cancer (P)

rbind(data.frame(reducedTerms1, g = 'Control'),
  data.frame(reducedTerms2, g = 'Cancer')) %>%
  as_tibble() -> reducedTerms

reducedTerms %>%
  group_by(g, parentTerm) %>%
  tally() %>%
  group_by(g) %>%
  arrange(desc(n)) %>%
  mutate(parentTerm = factor(parentTerm, 
    levels = unique(parentTerm))) %>%
  ggplot() +
  geom_col(aes(x = parentTerm, y = n, fill = g)) +
  labs(x = '', y = '# GO terms') +
  coord_flip() +
  facet_grid(g ~ ., scales = 'free_y') +
  theme_bw(base_size = 12, base_family = "GillSans") -> psave

ggsave(psave, filename = 'topGO_reduced_term.png', 
  path = path, 
  width = 12, height = 5, dpi = 250)

#
rbind(data.frame(topGO1.p, g = 'Control'),
  data.frame(topGO1.p, g = 'Cancer')) %>%
  as_tibble() -> topGO.p

# Multiple contrast ----
path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp'

mtd <- read.csv(file = paste0(path, '/metadata.csv'))

rownames(mtd) <- mtd$sample_id

colData <- mtd %>% arrange(match(sample_id, names(count)))

colData <- colData %>% mutate(g2 = ifelse(g1 == 'C', 'Control', g2))

colData <- mutate_if(colData, is.character, as.factor)

# using relevel, just specifying the reference level:

colData$g2 <- relevel(colData$g2, ref = "Control")

levels(colData$g2)

library(DESeq2)

count <- round(count)

table(colData$g2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ g2 ) # if not rep use design = ~ 1

# write_rds(dds, file = paste0(path, '/multiple_contrast_vs_control_dds.rds'))
# dds <- read_rds("/Cancer_vs_Control_dds.rds")
dds <- read_rds('~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/multiple_contrast_vs_control_dds.rds')

# dds <- DESeq(ddsFullCountTable)

contrast <- levels(colData(dds)$g2)

out <- list()

for(i in 2:length(contrast)) {
  j <- i
  cat('\nContrast', contrast[1],' and ', contrast[j],'\n')
  res <- get_res(dds, contrast[c(1,j)])
  out[[j]] <- res
  
}

do.call(rbind, out) -> res

# write_rds(res, file = '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/multiple_contrast_vs_control_res.rds')

res %>% distinct(sampleA, sampleB)

logfc_in <- 2

padj_in <- 0.05 

# volcano ----
mtd %>% group_by(Diagnosis, g2) %>% tally()

res.p %>% distinct(sampleA,sampleB)

res.p %>%
  mutate(pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = logFC, y = pvalue)) +
  geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = colors_fc) + 
  labs(x= expression(Log[2] ~ "Fold Change"), 
    y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top") +
  facet_grid(~ sampleB) -> pvol

ggsave(pvol, filename = 'multiple_contrast_control_vs_cancer.png',path = path,  
  width = 9.5, height = 3.5, dpi = 250)

#how up and down genes ----
alpha  = 0.05
lfcThreshold = 2

res.p <- prep_DE_data(res, alpha  = 0.05, lfcThreshold = 2)

write_rds(res.p, path = path, file = '')

title <- 'Differentially expressed genes'
subtitle <- expression(p[adj] < 0.05 ~ and ~ abs ~ (log[2] ~ FC) > 2)

res.p %>%
  filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  group_by(sampleB, lfcT) %>%
  tally() %>%
  ggplot(aes(x = sampleB, y = n, fill = lfcT)) +
  geom_col(color = 'black') +
  scale_fill_manual(name = "", values = c('#1f78b4', '#a6cee3')) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = '# Transcripts', x = '', subtitle = subtitle, title = title) +
  theme(legend.position = 'top') -> psave

ggsave(psave, filename = "multiple_contrast_control_vs_cancer_up-dpwn.png", path = path,
  width = 5, height = 5)

# AND VISUALIZE
  
# Pbar
res.p %>% arrange(logFC)

res.p %>%
  mutate(g = ifelse(logFC > 0, 'Control', 'Cancer')) %>%
  filter(g != 'Control') %>%
  group_by(cc, sampleB) %>% 
  tally() %>%
  filter(cc != 'NS') %>%
  ggplot(aes(x = sampleB, y = n)) +
  geom_bar(aes(fill = cc), 
    stat = 'identity', position = position_dodge2(width = 1)) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = '# DE Transcripts', x = '') +
  theme(legend.position = 'top') -> pbar

pbar + geom_text(aes(label = n), hjust = 0.5, family = "GillSans", 
  position = position_dodge2(width = 1),
  vjust = -0.5, size = 3) -> pbar

ggsave(pbar,
  filename = "multiple_contrast_control_vs_cancer_bar.png", path = path,
  width = 6.5, height = 5)

res.p %>%
  mutate(g = ifelse(logFC > 0, 'Control', 'Cancer')) %>%
  filter(g != 'Control') %>%
  filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  ggplot() +
  facet_wrap(~sampleB, scales = 'free_y', nrow = 1) +
  geom_histogram(aes(padj, fill = cc)) +
  # facet_wrap( cc ~ ., scales = 'free_x') +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = 'none') +
  labs(y = '# Transcripts') -> padjp

padjp + theme(panel.border = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.x = element_blank(),
  strip.text.x = element_blank()) -> padjp

library(patchwork)

pbar / padjp + patchwork::plot_layout(heights = c(3,1))-> savep

ggsave(savep, 
  filename = "multiple_contrast_control_vs_cancer_bar_hist.png", path = path, 
  width = 5, height = 5)

# heatmap AND pca
res.p %>% filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>% distinct(lfcT)

res.p %>%
  filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  # filter(sampleB %in% 'Cancer_IDC') %>%
  distinct(ids) %>%
  pull(ids) -> diffexpList

# mtd %>% filter(g2 %in% c('Control', 'Cancer_IDC')) %>% pull(sample_id) -> which_sam

str(diffexpList)

m <- DESeq2::counts(dds, normalized=T)
dim(DEcount <- m[rownames(m) %in% diffexpList,])

# PCA -----

# Heatmap -----

data = log2(DEcount+1)

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

gcolor <- structure(c('#d53e4f', '#3288bd'), names = c('P', 'C'))

PCAdf %>%
  mutate(sample_id = rownames(.)) %>%
  left_join(mtd) %>%
  ggplot(., aes(PC1, PC2)) +
  ggforce::geom_mark_hull(aes(group = as.factor(g1), label = as.factor(g1)), fill = 'grey') +
  geom_point(aes(color = g2), size = 5, alpha = 0.9) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_classic(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # coord_fixed(ratio = sd_ratio) +
  scale_color_brewer(palette = 'Set1') -> plotPca
  # scale_color_manual('', values = gcolor) -> plotPca


ggsave(plotPca, 
  filename = "multiple_contrast_control_vs_cancer_PCA.png", path = path, 
  width = 7, height = 5)


# DEcount <- DEcount[, colnames(DEcount) %in% which_sam]

# keep <- rownames(count) %in% diffexpList

# dim(count[keep,] -> DEcount)

sample_cor = cor(DEcount, method='pearson', use='pairwise.complete.obs')

# superheat::superheat()

dist.method <- 'euclidean'
linkage.method <- 'complete'

m <- log2(DEcount+1)

mtd %>% filter(g1 %in% c('C')) %>% pull(sample_id) -> control_sam
mtd %>% filter(g1 %in% c('P')) %>% pull(sample_id) -> cancer_sam


dim(m1 <- DEcount[, colnames(DEcount) %in% control_sam])
dim(m2 <- DEcount[, colnames(DEcount) %in% cancer_sam])

hc_control_samples <- hclust(dist(t(m1), method = dist.method),method = linkage.method)
hc_cancer_samples <- hclust(dist(t(m2), method = dist.method),method = linkage.method)

hc_sam_order <- hc_samples$labels[hc_samples$order]

hc_genes <- hclust(dist(m, method = dist.method), 
  method = linkage.method)

hc_genes_ord <- hc_genes$labels[hc_genes$order]


DEcount %>% 
  as_tibble(rownames = 'id') %>%
  pivot_longer(cols = colnames(DEcount), values_to = 'count', names_to = 'sample_id') %>%
  # mutate(sample_id = factor(sample_id, levels = hc_sam_order)) %>%
  # mutate(id = factor(id, levels = hc_genes_ord)) %>%
  left_join(mtd) %>%
  filter(count > 0) -> countLong

# preliminary bp
countLong %>% 
  # mutate(count = log2(count+1)) %>%
  ggplot(aes(x = sample_id, y = count)) +
  geom_boxplot()
# gcolor <- structure(c('#d53e4f', '#3288bd'), names = c('P', 'C'))


countLong %>% 
  distinct(sample_id, g1) %>% 
  mutate(col = ifelse(g1 %in% 'C', '#d53e4f', '#3288bd')) %>%
  arrange(match(sample_id, hc_sam_order)) -> coldf 

structure(coldf$col, names =  as.character(coldf$sample_id)) -> axis_col

library(ggh4x)



countLong %>%
  mutate(count = log2(count+1)) %>%
  ggplot(aes(x = sample_id, y = id, fill = count)) + 
  geom_raster() + 
  theme_classic(base_size = 12, base_family = "GillSans") +
  # ggsci::scale_fill_gsea(name = expression(Log[2]~'(cpm+1)'), reverse = T) +
  scale_fill_gradient2(midpoint = 0, low = '#ff7f00', mid = 'white', high = '#984ea3') +
  labs(x = '', y = '') +
  # ggh4x::scale_y_dendrogram(hclust = hc_genes) +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  theme(axis.text.x = element_text(angle = 90, 
    hjust = 1, vjust = 1, size = 7, color = axis_col), # 
    axis.ticks.length = unit(10, "pt"),
    # legend.position = 'top',
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) -> pheat

pheat


pheat + guides(fill = guide_colorbar(barheight = unit(4, "in"),
  ticks.colour = "black",
  frame.colour = "black",
  label.theme = element_text(size = 12))) -> pheat

# pheat + facet_grid(~ g1, scales = 'free_x', space = 'free_x')

ggsave(pheat, 
  filename = "multiple_contrast_control_vs_cancer_hetmap.png", path = path, 
  width = 8, height = 8)

# ven diagram

res.p

res.p %>% 
  distinct(ids, sampleB) %>% 
  group_by(GO.ID, Module) %>%
  count() %>%
  pivot_wider(names_from = 'Module', values_from = 'n', 
    values_fill = 0) -> euler

length(unique(euler$GO.ID)) == nrow(euler)

# rownames(euler) <- euler$GO.ID


eulerdf <- as.data.frame(euler)
eulerdf$GO.ID <- NULL

library(UpSetR)

labels <- names(eulerdf)
n <- length(labels)
grid.col <- labels # ggsci::pal_rickandmorty()(n)
names(grid.col) <- labels


# png(filename = paste0(dir, "intersect_families.png"), width = 780, height = 780, res = 150)
upset(eulerdf, number.angles = 0, point.size = 3.5, line.size = 2, 
  nintersects = 40, nsets = n,
  mainbar.y.label = "Intersections", sets.x.label = "G.O per Module", 
  text.scale = c(1.3, 2, 1, 2, 2, 2),
  order.by = c("freq", "degree"))

# dev.off()

