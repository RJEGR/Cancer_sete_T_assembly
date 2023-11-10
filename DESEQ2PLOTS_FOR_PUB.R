# VOLCANO PLOTS FROM DESEQ2

rm(list = ls()) 

if(!is.null(dev.list())) dev.off()

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

.RES <- read_rds(paste0(path, "CONTRAST_C_AND_D_FOR_PUB.rds"))

.RES <- do.call(rbind, .RES)

RES.P <- .RES %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)

RES.P %>% group_by(sampleB) %>% dplyr::count()

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

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
  annotate("text", x = -10, y = 30, label = "Control", family = "GillSans") +
  annotate("text", x = 10, y = 30, label = "Cancer", family = "GillSans") + 
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
  path = path, width = 10, height = 4, device = png, dpi = 300)

# AS DENSITY

RES.P %>%
  # filter(log2FoldChange < 0) %>%
  dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  ggplot(aes(x = log2(baseMean), y = log2FoldChange)) +
  facet_wrap(~ sampleB, scales = "free") +
  # geom_point() +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) +
  scale_fill_viridis_c()
  # geom_density_2d_filled(alpha = 0.5) +
  # geom_density_2d(linewidth = 0.25, colour = "black")
