
# LOAD MTD
# LOAD COUNT_MATRIX
# LOAD SEARCHED PROTEINS

rm(list = ls()) 

if(!is.null(dev.list())) dev.off()

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

ANNOT <- read_rds(paste0(path, "/WHICH_PROTEINS.rds"))

query.ids <- ANNOT %>% distinct(transcript_id) %>% pull()

swiss_f <- list.files(path, pattern = 'Trinotate.xls.blastx.tsv',  full.names = TRUE)

swissdf <- read_tsv(swiss_f) %>%
  mutate(orf = ifelse(is.na(protein), FALSE, TRUE)) %>%
  group_by(transcript, orf) %>%
  filter(identity == max(identity)) %>%
  ungroup() %>%
  distinct(transcript, uniprot, identity, name, genus, orf) %>%
  dplyr::rename("transcript_id"="transcript", "protein_name"="name")

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

mtd_f <- list.files(path, pattern = 'metadata.tsv',  full.names = TRUE)

dim(.raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))

# LOAD ONLY UP-expresed in CANCER

.RES <- read_rds(paste0(path, "CONTRAST_C_AND_D_FOR_PUB.rds")) %>%
  do.call(rbind, .) %>% 
  filter(padj < 0.05 & log2FoldChange < -20)

# ACCORDING TO 

# INCLUIDE ONLY EXCLUSIVELY W/ANNOT

nrow(P.RES.ANNOT <- .RES %>% count(transcript_id) %>% 
  filter(n == 1) %>% left_join(swissdf) %>%
  select(-n))

SAVEDF <- .RES%>%
  right_join(P.RES.ANNOT) %>%
  separate(uniprot, into = c("uniprot", "sp"), sep = "_") %>%
  select(-genus)
  

SAVEDF %>%  
  ggplot(aes(x = sampleB, y = uniprot, fill = log2FoldChange, color = log2FoldChange)) + 
  geom_tile( linewidth = 0.2)

write_tsv(SAVEDF, paste0(path, "EXCLUSIVE_CONTRAST_C_AND_D_FOR_PUB.tsv"))


# MUST MATCH 1392 exlcusive transcripts id

nrow(P.RES.ANNOT <- P.RES.ANNOT %>% drop_na(protein_name)) # w/ 1224 annot

# BIND w/ CANCER TRANSCRIPTS

P.RES.ANNOT <- rbind(P.RES.ANNOT, ANNOT) %>% 
  separate(uniprot, into = c("uniprot", "sp"), sep = "_") %>%
  distinct(transcript_id, protein_name, uniprot)


str(query.ids <- P.RES.ANNOT %>%  distinct(transcript_id) %>% pull()) # 1336

sum(keep <- rownames(.raw_count) %in% query.ids)

dim(COUNT <- as(.raw_count[keep,], "matrix") )

# THEN

MTD <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

keep <- MTD$LIBRARY_ID[!MTD$CONTRASTE_C %in% "Indiferenciado"]

dim(COUNT <- as(COUNT[,keep], "matrix") )

# GO TO COLLAPSED GENES HEATMAP

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

# IF COLLAPSED GENES ====

# (OPTIONAL) COLLAPSE BY 

P.RES.ANNOT %>% group_by(protein_name)

raw_count <- COUNT %>% as_tibble(rownames = 'transcript_id') %>%
  left_join(distinct(P.RES.ANNOT, transcript_id, protein_name)) %>%
  group_by(protein_name) %>%
  summarise_at(vars(colnames(COUNT)), sum)


COUNT <- raw_count %>% 
  select(any_of(colnames(COUNT))) %>%
  round() %>%
  as("matrix") 

rownames(COUNT) <- raw_count$protein_name

# THEN

MTD <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

# clustering

COUNT <- DESeq2::varianceStabilizingTransformation(round(COUNT))

# CLUSTERING SAMPLES

sample_cor = cor(COUNT, method='pearson', use='pairwise.complete.obs')

sample_dist = dist(sample_cor, method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

# CLUSTERING GENES

genes_dist = dist(COUNT, method='euclidean')

hc_genes = hclust(genes_dist, method='complete')

genes_order <- hc_genes$labels[hc_genes$order]

# sample_cor %>% 
#   as_tibble(rownames = 'LIBRARY_ID') %>%
#   pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
#   left_join(MTD) %>% 
#   mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) %>%
#   mutate(name = factor(name, levels = hc_order))-> sample_cor_long
# 

# heatmap(COUNT)
dim(COUNT)

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

DFANNOT <- P.RES.ANNOT %>% 
  separate(uniprot, into = c("uniprot", "sp"), sep = "_") %>%
  distinct(protein_name, uniprot)

labels <- DFANNOT %>%
  arrange(match(protein_name, genes_order)) %>%
  mutate(pullids = paste0(uniprot, " (", protein_name, ")")) %>%
  pull(pullids, name = protein_name)


P <- COUNT_LONG %>%
  filter(value > 1) %>%
  mutate(WHICH_CONTRAST = CONTRASTE_D) %>% # USE CONTRASTE_C OR CONTRASTE_D
  dplyr::mutate(CONTRASTE_C = dplyr::recode_factor(CONTRASTE_C, !!!recode_to_c)) %>%
  dplyr::mutate(CONTRASTE_D = dplyr::recode_factor(CONTRASTE_D, !!!recode_to_d)) %>%
  mutate(value = log2(value + 1))%>%
  ggplot(aes(x = LIBRARY_ID, y = protein_name, fill = value, color = value)) + 
  geom_tile( linewidth = 0.2) +
  ggsci::scale_fill_material(name = "Log2(Count)", "blue-grey") +
  ggsci::scale_color_material(name = "Log2(Count)", "blue-grey") +
  # scale_fill_viridis_c(option = "B", 
  #   name = "Log2(Count)", direction = -1, na.value = "white") +
  # scale_color_viridis_c(option = "B", 
  #   name = "Log2(Count)", direction = -1, na.value = "white") +
  # # ggh4x::facet_nested( ~ CONTRASTE_D+CONTRASTE_C, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = labels, label_size = 7)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "bottom",
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

P <- P + guides(
        fill = guide_colorbar(barwidth = unit(2, "in"),
          barheight = unit(0.05, "in"), label.position = "bottom",
          alignd = 0.5,
          ticks.colour = "black", ticks.linewidth = 0.5,
          frame.colour = "black", frame.linewidth = 0.5,
          label.theme = element_text(family = "GillSans", size = 7)))
    

# TOP PLOT ====

recode_to <- c(  `Control` = "Control",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Grado I` = "Stage I", `Grado II` = "Stage II", `Grado III` = "Stage III", 
  `Indiferenciado` = NA)

TOPDF <- COUNT_LONG %>%
  mutate(CONTRAST = CONTRASTE_C) %>%
  distinct(LIBRARY_ID, CONTRAST, CONTRASTE_D, EDAD) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  mutate(label = ifelse(CONTRASTE_D %in% "METASTASIS", "*", "")) %>%
  mutate(y = 1)

topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID, color = CONTRAST)) +
  geom_point(shape = 15, size = 2) +
  geom_text(aes(label = label),  vjust = -0.5, hjust = 0.5,
    color = "black", size = 2, family =  "GillSans") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::guide_dendro()
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  see::scale_color_pizza(name = "", reverse = T) +
  theme(legend.position = 'top',
    panel.border = element_blank(),
    plot.background = element_rect(fill='transparent', color = 'transparent'),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank()) 

library(patchwork)


psave <- topplot/ plot_spacer()/P + plot_layout(heights = c(0.6, -0.5, 5))


ggsave(psave, filename = 'HEATMAP_FOR_PUB_CONTRAST_C_D3.png', path = path, width = 10, height = 5.5, device = png, dpi = 300)


# ggsave(P, filename = 'HEATMAP_FOR_PUB_CONTRAST_C_D.png', path = path, width = 12, height = 4.5, device = png, dpi = 300)

# COUNT_LONG %>% select_if(!grepl("CONTRASTE", names(COUNT_LONG)))

# 
# LINE-PLOT BY CANCER RELATED GENES ====

dim(OUT1 <- read_csv(paste0(path, "/CONTRAST_C_GRADOS_HISTOLOGICOS.xls")))
dim(OUT2 <- read_csv(paste0(path, "/CONTRAST_D_METASTASIS_NO_METASTASIS.xls")))


query.ids <- rbind(OUT1, OUT2) %>%
  filter(padj < 0.05) %>%
  dplyr::rename("transcript_id" = "Name") %>%
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

COUNT <- DESeq2::vst(as(round(.raw_count), "matrix"))
 
# CLUSTERING SAMPLES ====

keep <- MTD$LIBRARY_ID[!MTD$CONTRASTE_C %in% "Indiferenciado"]

dim(COUNT <- as(COUNT[,keep], "matrix") )

sample_cor = cor(COUNT, method='pearson', use='pairwise.complete.obs')

sample_dist = dist(sample_cor, method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]


# heatmap(sample_cor, col = cm.colors(12))

# sample_cor %>% 
#   as_tibble(rownames = 'LIBRARY_ID') %>%
#   pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
#   left_join(MTD) %>% 
#   mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) %>%
#   mutate(name = factor(name, levels = hc_order))->     frame.colour = "black", frame.linewidth = 0.5,
#     label.theme = element_text(family = "GillSans", size = 7)))) 

library(ggh4x)

recode_to_c <- c(  `Control` = "(A) Control",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Grado I` = "(B) Stage I", `Grado II` = "(C) Stage II", `Grado III` = "(D) Stage III",
  `Indiferenciado` = NA)

recode_to_d <- c(  `Control` = "",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Indiferenciado` = NA)

sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, color = cor, fill = cor)) + 
  geom_tile() + 
  ggsci::scale_color_material(name = "", "blue-grey") +
  ggsci::scale_fill_material(name = "", "blue-grey") +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(
    y.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  labs(x = '', y = '') +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(legend.position = 'bottom',
    axis.text.x = element_text(angle = 90,
      hjust = 1, vjust = 0.5, size = 4)
    ) -> p

p <- p +  guides(
  colour = guide_colorbar(barwidth = unit(2, "in"),
    barheight = unit(0.05, "in"), label.position = "bottom",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(family = "GillSans", size = 7))) 

p <- p + theme(
  # strip.background = element_rect(fill = 'white', color = 'white'),
  panel.border = element_blank(),
  plot.background = element_rect(fill='transparent', color = 'transparent'),
  plot.margin = unit(c(0,0,0,0), "pt"),
  panel.grid.minor = element_blank(),
  # axis.ticks.x = element_blank(),
  panel.grid.major = element_blank())

p

# ADD TOP dendograp  ====
# scales::show_col(see::palette_pizza("default")(5))
#
recode_to <- c(  `Control` = "Control",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Grado I` = "Stage I", `Grado II` = "Stage II", `Grado III` = "Stage III", 
  `Indiferenciado` = NA)

topplot <- sample_cor_long %>%
  mutate(CONTRAST = CONTRASTE_C) %>%
  distinct(LIBRARY_ID, CONTRAST, CONTRASTE_D) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  mutate(label = ifelse(CONTRASTE_D %in% "METASTASIS", "*", "")) %>%
  ggplot(aes(y = 1, x = LIBRARY_ID, color = CONTRAST)) +
  geom_point(shape = 15) +
  geom_text(aes(label = label),  vjust = -0.5, hjust = 0.5,
    color = "black", size = 2, family =  "GillSans") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::guide_dendro()
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  see::scale_color_pizza(name = "", reverse = T) +
  theme(legend.position = 'top',
    panel.border = element_blank(),
    plot.background = element_rect(fill='transparent', color = 'transparent'),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank()) 

library(patchwork)


psave <- topplot/ plot_spacer()/p + plot_layout(heights = c(0.6, -0.5, 5))

ggsave(psave, filename = 'SAMPLE_HEATMAP_FOR_PUB_STAR.png', path = path, width = 3.7, height = 4.5, device = png, dpi = 300)

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
  group_by(protein_name) 
  # rstatix::cor_test(`Grado I`, `Grado II`)

# library(rstatix)
# 
# sample_cor %>%
#   cor_reorder() %>%
#   # pull_lower_triangle() %>%
#   cor_plot(label = F)

# PCA ====

COUNT <- DESeq2::vst(as(round(.raw_count), "matrix"))

keep <- MTD$LIBRARY_ID[!MTD$CONTRASTE_C %in% "Indiferenciado"]

dim(COUNT <- as(COUNT[,keep], "matrix") )

# COUNT <- DESeq2::varianceStabilizingTransformation(round(as(.raw_count[keep,], "matrix")))


# ncol(data <- log2(COUNTS+1))

PCA = prcomp(t(COUNT), center = T, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

# Visualize eigenvalues 
# Show the percentage of variances explained by each principal component.

barplot(PCA$sdev)


k <- 4

PCAdf %>% 
  dist(method = "euclidean") %>% 
  hclust() %>% 
  cutree(., k) %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>% 
  mutate(cluster = paste0('C', value)) %>% 
  dplyr::select(-value) -> hclust_res


# scales::show_col(see::pizza_colors())

col_values <- c("#768947", "#D3BEAA", "#CE3722", "#642118")

recode_to <- c(  `Control` = "Control",
  `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
  `Grado I` = "Stage I", `Grado II` = "Stage II", `Grado III` = "Stage III", 
  `Indiferenciado` = NA)
  
PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(MTD) %>%
  mutate(sample_group = CONTRASTE_C) %>%
  dplyr::mutate(sample_group = dplyr::recode_factor(sample_group, !!!recode_to)) %>%
  ggplot(., aes(PC1, PC2)) +
  # coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(CONTRASTE_A)),
    fill = 'grey89', color = NA) +
  geom_point(size = 5, alpha = 0.7, aes(color = sample_group)) +
  # geom_text( family = "GillSans",
  #   mapping = aes(label = paste0(hpf, " hpf")), size = 2.5) +
  labs(caption = '') +
  # ylim(-250, 250) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  see::scale_color_pizza(name = "", reverse = T) +
  # scale_color_manual("", values = col_values) +
  theme_classic(base_family = "GillSans", base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') -> p


ggsave(p, filename = 'SAMPLE_PCA_FOR_PUB.png', path = path, width = 4, height = 4, device = png, dpi = 300)

library(latticeExtra)

sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

colData <- MTD %>% arrange(match(LIBRARY_ID, hc_order))

identical(hc_order, colData$LIBRARY_ID)

groups <- factor(colData$CONTRASTE_D)

hc1 <- as.dendrogram(hc_samples)
ord.hc1 <- order.dendrogram(hc1)
hc2 <- reorder(hc1, groups[ord.hc1])
ord.hc2 <- order.dendrogram(hc2)
colors <- trellis.par.get("superpose.polygon")$col


levelplot(t(scale(sample_dist))[, ord.hc2], scales = list(x = list(rot = 90)),
  colorkey = F, legend = list(right = list(fun = dendrogramGrob,
    args = list(x = hc2, ord = ord.hc2, side = "right",
      size = 10, size.add = 0.5, add = list(rect = list(col = "transparent",
        fill = colors[groups])),
      type = "rectangle"))))



