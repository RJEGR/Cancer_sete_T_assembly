# VOLCANO PLOTS FROM DESEQ2

rm(list = ls()) 

if(!is.null(dev.list())) dev.off()

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

.RES <- read_rds(paste0(path, "CONTRAST_C_AND_D_FOR_PUB.rds"))

.RES <- do.call(rbind, .RES)

RES.P <- .RES %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)

RES.P %>% group_by(sampleB) %>% dplyr::count()

# RES.P %>%
#   ggplot(aes(log10(baseMeanB), log2FoldChange)) +
#   facet_grid(~ sampleB) +
#   geom_point()

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

# IDS

dim(OUT1 <- read_csv(paste0(path, "../CONTRAST_C_GRADOS_HISTOLOGICOS.xls")))
dim(OUT2 <- read_csv(paste0(path, "../CONTRAST_D_METASTASIS_NO_METASTASIS.xls")))


query.ids <- rbind(OUT1, OUT2) %>%
  filter(padj < 0.05) %>%
  dplyr::rename("transcript_id" = "Name") %>%
  distinct(transcript_id) %>%
  pull()


# 

sigfc <- "p - value ~ and ~ log[2] ~ FC"
pv <- "p-value"
fc <- "Log[2] ~ FC"

.RES$cc <- 'NS'
.RES[which(abs(.RES$log2FoldChange) > 2), 'cc'] <- fc
.RES[which(abs(.RES$padj) <= 0.05), 'cc'] <- pv
.RES[which(.RES$padj <= 0.05 & abs(.RES$log2FoldChange) > 2), 'cc'] <- sigfc

colors_fc <- c("forestgreen",  "#4169E1", "red2")

c("#DADADA", "#D4DBC2")

recode_to <- c(`METASTASIS` = "(A) Metastasis", `NO METASTASIS`= "(B) No metastasis",
  `Grado I` = "(C) Stage I", `Grado II` = "(D) Stage II", `Grado III` = "(E) Stage III")


p <- .RES %>% 
  filter(cc != "NS") %>%
  dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  facet_grid(~ sampleB) +
  geom_rect(
    aes(xmin=-30, xmax = -1, ymin = 1, ymax = Inf), fill = '#DADADA') +
  geom_rect(
    aes(xmin=1, xmax = 30, ymin = 1, ymax = Inf), fill = '#D4DBC2') +
  labs(x= expression(Log[2] ~ "Fold Change"), 
    y = expression(-Log[10] ~ "padj")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  geom_abline(slope = 0, intercept = -log10(0.05), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 1, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = -1, linetype="dashed", alpha=0.5) +
  annotate("text", x = -15, y = 30, label = "Control", family = "GillSans") +
  annotate("text", x = 15, y = 30, label = "Cancer", family = "GillSans") + 
  geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = colors_fc) +
  xlim(-30, 30)

p <- p + theme(legend.position = "top",
  panel.border = element_blank(),
  strip.background = element_rect(fill = 'grey89', color = 'white'),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  # panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank(),
  # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
  # axis.text.y = element_text(angle = 0, size = 5),
  axis.text.x = element_text(angle = 0))

ggsave(p, filename = 'DESEQ2VOLCANO.png', 
  path = path, width = 10, height = 3.5, device = png, dpi = 300)

# AS DENSITY

p2 <- RES.P %>%
  # filter(log2FoldChange < 0) %>%
  dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  ggplot(aes(x = log2(baseMean), y = log2FoldChange)) +
  facet_grid(~ sampleB, scales = "free") +
  geom_point(aes(alpha = -log10(padj)), color = 'grey') +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # scale_fill_viridis_c() +
  geom_density_2d(aes(color = ..level..), linewidth = 0.5) +
  scale_color_viridis_c() +
  # stat_density_2d(
  #   geom = "raster",
  #   aes(fill = after_stat(density)),
  #   contour = FALSE
  # ) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(angle = 0))

# p2

ggsave(p2, filename = 'DESEQ2DENSITY.png', 
  path = path, width = 7, height = 3.5, device = png, dpi = 300)

# ONLY CANCER ASSOCIATED ====

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

ANNOT <- read_rds(paste0(path, "/WHICH_PROTEINS.rds"))


colors_fc <- c("forestgreen",  "#4169E1", "red2", "grey89")

colors_fc <- structure(colors_fc, names = c(fc, sigfc, pv, "NS"))

ANNOT <- ANNOT %>% 
  separate(uniprot, into = c("uniprot", "sp"), sep = "_") %>%
  distinct(transcript_id, uniprot) 

ANNOT %>% distinct(uniprot)

.RES %>%
  mutate(cc = ifelse(padj < 0.05, "p-value","")) %>%
  left_join(ANNOT, by = "transcript_id") %>%
  filter(transcript_id %in% query.ids) %>%
  dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  ggplot(aes(x = log2(baseMean), y = log2FoldChange)) + 
  # geom_rect(
    # aes(ymin=-0, ymax = 4, xmin = 0, xmax = Inf), fill = 'grey89') +
  geom_rect(
    aes(ymin=-0, ymax = -4, xmin = 0, xmax = Inf), fill = 'grey89') + # D4DBC2
  annotate("text", x = 6.5, y = 3, label = "Control", family = "GillSans", color = "grey40") +
  annotate("text", x = 6.5, y = -3, label = "Cancer", family = "GillSans",  color = "grey40") +
  ggrepel::geom_text_repel(aes(color = cc, label = uniprot), size = 2, family = "GillSans", max.overlaps = 50) +
  # geom_text(aes(color = cc, label = uniprot), family = "GillSans", size = 3) +
  facet_grid(~ sampleB) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  # geom_abline(slope = 0, intercept = -1, linetype="dashed", alpha=0.5) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = 1, linetype="dashed", alpha=0.5) +
  # geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = c("grey20","red")) +
  ylim(-4, 4) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "none",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0)) -> p

ggsave(p, filename = 'SEARCH_GENES_FOR_PUB.png', path = path, width = 10, height = 3.5, device = png, dpi = 300)


# also

SUBSET_RES <- .RES %>%
  mutate(cc = ifelse(padj < 0.05, "p-value","")) %>%
  left_join(ANNOT, by = "transcript_id") %>%
  filter(transcript_id %in% query.ids) %>%
  dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to))

p2 +
  # ylim(c(-10,10))
  geom_rect(
    aes(ymin=-30, ymax = -1, xmin = 1, xmax = Inf), fill = '#DADADA', alpha = 0.02) +
  geom_rect(
    aes(ymin=1, ymax = 30, xmin = 1, xmax = Inf), fill = '#DADADA', alpha = 0.02) +
  ggrepel::geom_text_repel(data = SUBSET_RES, aes(label = uniprot), 
    size = 2, family = "GillSans", max.overlaps = 100) 

# HEATMAP OF 
# CONTINUE w/ SEARCH_GENES_HEATMAP_FOR_PUB 

RES.P %>%
  filter(log2FoldChange < -2) %>%
  dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  ggplot(aes(log2FoldChange, fill = sampleB)) +
  geom_histogram(color = "grey20") +
  facet_grid(~ sampleB, scales = "free_x") +
  # stat_ecdf() +
  # xlim(c(-30, -2)) +
  labs(y = "Number of transcripts (Up-expressed") +
  see::scale_color_pizza(name = "", reverse = T) +
  see::scale_fill_pizza(name = "", reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  see::scale_color_pizza(name = "", reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(angle = 0)) -> p


ggsave(p, filename = 'HISTOGRAM_UPEXPRESSED_GENES_FOR_PUB.png', path = path, width = 10, height = 3.5, device = png, dpi = 300)
