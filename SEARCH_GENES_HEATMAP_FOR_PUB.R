
# LOAD MTD
# LOAD COUNT_MATRIX
# LOAD SEARCHED PROTEINS

rm(list = ls()) 

if(!is.null(dev.list())) dev.off()

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

ANNOT <- read_rds(paste0(path, "/WHICH_PROTEINS.rds"))

query.ids <- ANNOT %>% distinct(transcript_id) %>% pull()

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

mtd_f <- list.files(path, pattern = 'metadata.tsv',  full.names = TRUE)

dim(.raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))

keep <- rownames(.raw_count) %in% query.ids

dim(COUNT <- as(.raw_count[keep,], "matrix") )


# THEN

MTD <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

# clustering


COUNT <- DESeq2::varianceStabilizingTransformation(round(COUNT))

# CLUSTERING SAMPLES

sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

# CLUSTERING GENES

genes_dist = dist(COUNT, method='euclidean')

hc_genes = hclust(genes_dist, method='complete')

genes_order <- hc_genes$labels[hc_genes$order]

recode_to <- ANNOT %>% filter(transcript_id %in% genes_order) %>% 
  distinct(transcript_id, protein_name)

recode_to <- structure(recode_to$protein_name, names = recode_to$transcript_id)

identical(sort(names(recode_to)),sort(genes_order))

genes_order <- recode_to[match(genes_order, names(recode_to))]

identical(names(genes_order),  hc_genes$labels[hc_genes$order])


COUNT %>% 
  as_tibble(rownames = 'transcript_id') %>%
  pivot_longer(-transcript_id, names_to = "LIBRARY_ID") %>%
  left_join(MTD) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> COUNT_LONG


COUNT_LONG <- COUNT_LONG %>%
  left_join(distinct(ANNOT, transcript_id, protein_name), by = "transcript_id") 
  # mutate(Family = ifelse(!grepl("Cluster", Family), toupper(Family), Family)) %>%
  # group_by(LIBRARY_ID, transcript_id) %>% 
  # summarise(value = sum(value)) 

library(ggh4x)

COUNT_LONG %>% distinct(CONTRASTE_C)
COUNT_LONG %>% distinct(CONTRASTE_D)

recode_to <- c(  `Control` = "(A) Control",
  `METASTASIS` = "(B) Metastasis", `NO METASTASIS`= "(C) No metastasis",
  `Grado I` = "(B) Stage I", `Grado II` = "(C) Stage II", `Grado III` = "(D) Stage III", 
  `Indiferenciado` = NA)

P <- COUNT_LONG %>%
  mutate(WHICH_CONTRAST = CONTRASTE_C) %>% # USE CONTRASTE_C OR CONTRASTE_D
  dplyr::mutate(WHICH_CONTRAST = dplyr::recode_factor(WHICH_CONTRAST, !!!recode_to)) %>%
  drop_na(WHICH_CONTRAST) %>%
  ggplot(aes(x = LIBRARY_ID, y = transcript_id, fill = log2(value))) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", name = "Log2(Count)", direction = -1, na.value = "white") +
  ggh4x::facet_nested( ~ WHICH_CONTRAST, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = genes_order, label_size = 4.5)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
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
    panel.grid.major.x = element_blank()) 

ggsave(P, filename = 'HEATMAP_FOR_PUB_CONTRAST_C.png', path = path, width = 10, height = 7, device = png, dpi = 300)

# IF COLLAPSED GENES

# (OPTIONAL) COLLAPSE BY 

raw_count <- .raw_count[keep,] %>% as_tibble(rownames = 'transcript_id') %>%
  left_join(distinct(ANNOT, transcript_id, protein_name)) %>%
  group_by(protein_name) %>%
  summarise_at(vars(names(.raw_count)), sum)


COUNT <- raw_count %>% 
  select(any_of(names(.raw_count))) %>%
  round() %>%
  as("matrix") 

rownames(COUNT) <- raw_count$protein_name

# THEN

MTD <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

# clustering


COUNT <- DESeq2::varianceStabilizingTransformation(round(COUNT))

# CLUSTERING SAMPLES

sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

# CLUSTERING GENES

genes_dist = dist(COUNT, method='euclidean')

hc_genes = hclust(genes_dist, method='complete')

genes_order <- hc_genes$labels[hc_genes$order]

COUNT %>% 
  as_tibble(rownames = 'protein_name') %>%
  pivot_longer(-protein_name, names_to = "LIBRARY_ID") %>%
  left_join(MTD) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> COUNT_LONG

library(ggh4x)



recode_to_c <- c(  `Control` = "(A) Control",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Grado I` = "(B) Stage I", `Grado II` = "(C) Stage II", `Grado III` = "(D) Stage III",
  `Indiferenciado` = NA)

recode_to_d <- c(  `Control` = "",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Indiferenciado` = NA)


P <- COUNT_LONG %>%
  mutate(WHICH_CONTRAST = CONTRASTE_D) %>% # USE CONTRASTE_C OR CONTRASTE_D
  dplyr::mutate(CONTRASTE_C = dplyr::recode_factor(CONTRASTE_C, !!!recode_to_c)) %>%
  dplyr::mutate(CONTRASTE_D = dplyr::recode_factor(CONTRASTE_D, !!!recode_to_d)) %>%
  drop_na(WHICH_CONTRAST) %>%
  ggplot(aes(x = LIBRARY_ID, y = protein_name, fill = log2(value))) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", 
    name = "Log2(Count)", direction = -1, na.value = "white") +
  ggh4x::facet_nested( ~ CONTRASTE_D+CONTRASTE_C, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = genes_order, label_size = 7)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 

# P
ggsave(P, filename = 'HEATMAP_FOR_PUB_CONTRAST_C_D.png', path = path, width = 12, height = 4.5, device = png, dpi = 300)

# COUNT_LONG %>% select_if(!grepl("CONTRASTE", names(COUNT_LONG)))

# 
# LINE-PLOT BY CANCER RELATED GENES ====

OUT1 <- read_tsv(paste0(path, "/CONTRAST_C_GRADOS_HISTOLOGICOS.xls"), )
OUT2 <- read_csv(paste0(path, "/CONTRAST_D_METASTASIS_NO_METASTASIS.xls"), )


query.ids <- rbind(OUT1, OUT2) %>%
  filter(padj < 0.05) %>%
  distinct(transcript_id) %>%
  pull()


# barplot(cnts <- DESeq2::counts(dds, normalized = T, replaced = F)[query.ids,])
# barplot(cnts <- DESeq2::varianceStabilizingTransformation(round(cnts)))

keep <- rownames(.raw_count) %in% query.ids

dim(cnts <- as(.raw_count[keep,], "matrix") )

DF <- cnts %>% as_tibble(rownames = "transcript_id") %>% 
  pivot_longer(-transcript_id, names_to = "LIBRARY_ID", values_to = "count")


DF <- ANNOT %>% 
  separate(uniprot, into = c("uniprot", "sp"), sep = "_") %>%
  distinct(transcript_id, uniprot) %>% 
  right_join(DF, by = "transcript_id") %>% 
  group_by(LIBRARY_ID, uniprot) %>% 
  summarise(count = sum(count)) %>%
  left_join(MTD, by = "LIBRARY_ID")

fun.data.trend <- "mean_sdl" # "mean_cl_boot", "mean_se"

recode_to <- c(`CON_CANCER` = "Cancer",`Control` = "(A) Control",
  `METASTASIS` = "(B) Metastasis", `NO METASTASIS`= "(C) No metastasis",
  `Grado I` = "(B) Stage I", `Grado II` = "(C) Stage II", `Grado III` = "(D) Stage III", 
  `Indiferenciado` = NA)


DF %>%
  ungroup() %>%
  mutate(EDAD = as.numeric(EDAD)) %>%
  mutate(Age = "") %>%
  mutate(Age = ifelse(between(EDAD, 0, 19), "0-19", Age)) %>%
  mutate(Age = ifelse(between(EDAD, 21, 29), "20-29", Age)) %>%
  mutate(Age = ifelse(between(EDAD, 30, 39), "20-39", Age)) %>%
  mutate(Age = ifelse(between(EDAD, 40, 69), "40-69", Age)) %>%
  mutate(Age = ifelse(between(EDAD, 70, 90), "70-90", Age)) %>%
  # ungroup() %>% count(Age) %>%
  filter(count > 0) %>%
  dplyr::mutate(CONTRASTE_A = dplyr::recode_factor(CONTRASTE_A, !!!recode_to, .ordered = T)) %>%
  dplyr::mutate(CONTRASTE_C = dplyr::recode_factor(CONTRASTE_C, !!!recode_to, .ordered = T)) %>%
  drop_na(CONTRASTE_C) %>%
  # mutate(y = count) %>%
  mutate(y  = log2(round(count)+1)) %>%
  ggplot(aes(x = CONTRASTE_A, y = y, 
    fill = CONTRASTE_C,
    color = CONTRASTE_C,
    group = CONTRASTE_C)) +
  facet_wrap(~ uniprot, scales = "free_y") +
  geom_point(aes(shape = Age), position=position_jitterdodge(dodge.width=0.9)) +
  # stat_summary(fun = mean, geom = "point") +
  stat_summary(position = position_dodge(width=0.9), 
    fun.data = fun.data.trend, linewidth = 0.7, size = 0.7, alpha = 0.7) +
  labs(y = "Read Count (log2)", x = "") +
  # scale_color_manual("", values = scale_col_pH) +
  # scale_fill_manual("", values =  scale_col_pH) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = "top",
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 0, size = 10))

# sample-heatmap

COUNT <- DESeq2::vst(round(COUNT))

# CLUSTERING SAMPLES ====

dim(COUNT <- as(.raw_count, "matrix") )

sample_cor = cor(COUNT, method='pearson', use='pairwise.complete.obs')

sample_dist = dist(sample_cor, method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

# heatmap(sample_cor, col = cm.colors(12))

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(MTD) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> sample_cor_long

library(ggh4x)

recode_to_c <- c(  `Control` = "(A) Control",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Grado I` = "(B) Stage I", `Grado II` = "(C) Stage II", `Grado III` = "(D) Stage III",
  `Indiferenciado` = NA)

recode_to_d <- c(  `Control` = "",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Indiferenciado` = NA)

sample_cor_long %>%
  # dplyr::mutate(CONTRASTE_C = dplyr::recode_factor(CONTRASTE_C, !!!recode_to_c)) %>%
  # dplyr::mutate(CONTRASTE_D = dplyr::recode_factor(CONTRASTE_D, !!!recode_to_d)) %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) + 
  # geom_tile(color = 'white', size = 0.2) +
  geom_tile(color = 'white') +
  # ggh4x::facet_nested( ~ CONTRASTE_D+CONTRASTE_C, nest_line = F, scales = "free", space = "free") +
  # geom_text(aes(label = cor), color = 'white') +
  scale_fill_viridis_c(option = "B", name = "", direction = -1) +
  # scale_x_discrete(position = 'top') +
  # scale_y_discrete(position = "right") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left") +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 4.5),
         x.sec = guide_axis_manual(labels = hc_order, label_size = 4.5)) +
  guides(fill = guide_legend(title = "", nrow = 1)) +
  labs(x = '', y = '') +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(legend.position = 'top',
    axis.text.x = element_text(angle = 90,
      hjust = -0.15, vjust = 1)) -> p

p <- p + theme(strip.background = element_rect(fill = 'white', color = 'white'),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks.x = element_blank(),
  panel.grid.major = element_blank())



ggsave(p, filename = 'SAMPLE_HEATMAP_FOR_PUB.png', path = path, width = 5, height = 5, device = png, dpi = 300)
# OR BY GROUP

COUNT %>% 
  as_tibble(rownames = 'protein_name') %>%
  pivot_longer(-protein_name, names_to = "LIBRARY_ID") %>%
  left_join(MTD) -> COUNT_LONG

COUNT_LONG %>%
  group_by(CONTRASTE_C, CONTRASTE_D, protein_name) %>%
  summarise(value = sum(value)) -> COUNT_LONG

COUNT_LONG %>%
  pivot_wider(names_from = CONTRASTE_C, values_from = value, values_fill = 0) %>%
  group_by(protein_name) %>%
  # rstatix::cor_test(`Grado I`, `Grado II`)

# library(rstatix)
# 
# sample_cor %>%
#   cor_reorder() %>%
#   # pull_lower_triangle() %>%
#   cor_plot(label = F)

