# SEARCH Cytokine ----

# write_excel_csv(res.p.out, file = file_name)
library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

swiss_f <- list.files(path, pattern = 'Trinotate.xls.blastx.tsv',  full.names = TRUE)

swissdf <- read_tsv(swiss_f) %>%
  mutate(orf = ifelse(is.na(protein), FALSE, TRUE)) %>%
  group_by(transcript, orf) %>%
  filter(identity == max(identity)) %>%
  ungroup() %>%
  distinct(transcript, uniprot, identity, name, genus, orf) %>%
  dplyr::rename("transcript_id"="transcript", "protein_name"="name")


out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"

file_name <- "ALL_MULTIPLE_CONTRAST"

file_name <- paste0(out_path, paste(file_name,collapse = '_'),
  '_Annot_count_down_up_genes.xls')

res.p.out <- read_csv(file_name)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

mtd_f <- list.files(path, pattern = 'metadata.tsv',  full.names = TRUE)

mtd <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

dim(raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))


which_col <- names(raw_count)

# global presence of Cytokine -----

global_cyk_swiss <- swissdf %>%  filter(grepl("Cytokine", protein_name)) %>% distinct(transcript_id, .keep_all = T)

query.genes <- global_cyk_swiss$transcript_id

sum(keep <- rownames(raw_count) %in% query.genes)

data <- raw_count[keep,]

dist.method <- 'euclidean'; linkage.method <- 'complete'

m <- log2(data+1)

hc_genes <- hclust(dist(m, method = dist.method), 
  method = linkage.method)

hc_genes_ord <- hc_genes$labels[hc_genes$order]

#

sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method= dist.method)

hc_samples = hclust(sample_dist, method= linkage.method)


hc_sam_order <- hc_samples$labels[hc_samples$order]

# hc_order <- hc_samples$labels[hc_samples$order]

data %>% 
  as_tibble(rownames = 'transcript_id') %>%
  left_join(global_cyk_swiss) %>%
  pivot_longer(cols = colnames(data), values_to = 'count', names_to = 'LIBRARY_ID') %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_sam_order)) %>%
  mutate(transcript_id = factor(transcript_id, levels = hc_genes_ord)) %>%
  left_join(mtd, by = "LIBRARY_ID") %>%
  filter(count > 0) -> countLong

library(ggh4x)

countLong %>%
  mutate(count = log2(count+1)) %>%
  ggplot(aes(x = LIBRARY_ID, y = transcript_id, fill = count)) + # protein_name
  geom_raster() + 
  theme_classic(base_size = 12, base_family = "GillSans") +
  # ggsci::scale_fill_gsea(name = expression(Log[2]~'(cpm+1)'), reverse = T) +
  scale_fill_gradient2(name = expression(Log[2]~'(x+1)'),
    midpoint = 0, low = '#ff7f00', mid = 'white', high = '#984ea3') +
  labs(x = '', y = '') +
  # ggh4x::scale_y_dendrogram(hclust = hc_genes) +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  theme(axis.text.x = element_text(angle = 90, 
    hjust = 1, vjust = 1, size = 7), # 
    axis.ticks.length = unit(10, "pt"),
    # legend.position = 'top',
    panel.border = element_blank(),
    # axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) -> p

p <- p +facet_grid(protein_name ~ ., scales = "free")

p

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"
ggsave(p, path = out_path, filename = "Cytokine_heatmap.png", width = 14, height = 10) 

# # pheatmap ----

annotation_col <- mtd %>% 
  data.frame(row.names = .$LIBRARY_ID) %>%
  select(CONTRASTE_A, CONTRASTE_B, CONTRASTE_C,CONTRASTE_D)


# IF LABEL ROW BY PROTEIN NAME

prot_data <- data %>% 
  as_tibble(rownames = 'transcript_id') %>%
  left_join(global_cyk_swiss) %>% 
  group_by(protein_name) %>% summarise_at(vars(names(data)), sum) %>%
  data.frame(row.names = .$protein_name) %>% select(-protein_name)
  
pheatmap::pheatmap(log2(prot_data+1),
  annotation_col = annotation_col)

# ELSE
annotation_row <- data %>% 
  as_tibble(rownames = 'transcript_id') %>%
  left_join(global_cyk_swiss) %>% 
  select(-any_of(names(data))) %>% 
  data.frame(row.names = .$transcript_id) %>%
  select( protein_name) %>%
  mutate_if(is.character, as.factor)

ann_colors = list(
  CONTRASTE_A = c(Control = "blue", CON_CANCER = "red")
#   CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
#   GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)

library(pheatmap)

log2(data+1) %>% pheatmap::pheatmap(
  annotation_row = annotation_row, 
  row_split = annotation_row$protein_name,
  annotation_col = annotation_col,
  annotation_colors = ann_colors
  # cluster_rows = T,
  # cluster_cols = T,
  # cutree_cols = 3,
  # cellwidth = 7,
  # cex = 1,
  # border_color = T
)

# m %>%
#   # data.frame(row.names = rownames(my_gene_col), .) %>% 
#   pheatmap(., 
#     annotation_row = my_gene_col, 
#     annotation_colors = annotation_colors,
#     cluster_rows = T,
#     cluster_cols = T,
#     # cutree_rows = 3,
#     # cutree_cols = 3,
#     na_col = "white",
#     cellwidth = 30,
#     cex = 1,
#     border_color = T,
#     angle_col = c("45"),
#     # legend_breaks = -1:4
#     # annotation_name_row = T,
#     color = pal,
#     fontsize = 14, 
#     show_rownames = T


# after DE test
res.p.out %>% filter(grepl("Cytokine", protein_name)) %>% 
  distinct(transcript_id, .keep_all = T) %>% view()
  select(transcript_id, any_of(which_col)) %>% as.data.frame() -> data

# data %>% view()
rownames(data) <- data$transcript_id

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(mtd) -> PCAdf

PCAdf %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = CONTRASTE_A), size = 5, alpha = 0.7) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # coord_fixed(ratio = sd_ratio) +
  guides(color=guide_legend("",nrow=2)) +
  see::scale_color_metro() #-> ps


# Heatmap -----


sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]


# SAMPLE COR ----
sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(mtd) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> sample_cor_long


library(ggh4x)

sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) + 
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
    hjust = 1, vjust = 1, size = 7),
    # axis.text.y = element_text(color = axis_col),
    axis.ticks.length = unit(5, "pt"))

# GENES MAP

dist.method <- 'euclidean'
linkage.method <- 'complete'

m <- log2(data+1)

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



# out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"
# ggsave(ps, path = out_path, filename = file_name, width = 12) 

