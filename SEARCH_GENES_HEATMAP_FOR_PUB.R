
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

dim(raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))

keep <- rownames(raw_count) %in% query.ids
 
dim(COUNT <- as(raw_count[keep,], "matrix") )

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

recode_to <- ANNOT %>% filter(transcript_id %in% genes_order) %>% distinct(transcript_id, protein_name)

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



