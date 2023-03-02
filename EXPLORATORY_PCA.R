# PCA using vsd
library(tidyverse)
library(DESeq2)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

rds_f_l <- list.files(path, pattern = 'CONTRAST_',  full.names = TRUE)

# file_name <- strsplit(basename(dds_file), "_")[[1]][1:2]
# file_name <- paste0(out_path, paste(file_name,collapse = '_'),'_Annot_count_down_up_genes.xls')

prep_pca_df <- function(dds_file) {
  
  group <- strsplit(basename(dds_file), "_")[[1]][1:2]
  group <- paste(group,collapse = '_')
  
  dds <- read_rds(dds_file)
  
  vst <- vst(dds)
  
  count <- assay(vst)
  colData <- colData(dds) %>% as_tibble()
  
  # data = log2(DEcount+1)
  data <- count
  
  PCA = prcomp(t(data), center = FALSE, scale. = FALSE)
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], group)
  
  # gcolor <- structure(c('#d53e4f', '#3288bd'), names = c('P', 'C'))
  
  PCAdf %>%
    mutate(LIBRARY_ID = rownames(.)) %>%
    left_join(colData) -> PCAdf
  
  return(PCAdf)
  
  PCAdf %>%
    ggplot(., aes(PC1, PC2)) +
    geom_point(aes(color = Design), size = 5, alpha = 0.7) +
    labs(caption = '') +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme_classic(base_family = "GillSans", base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
    # coord_fixed(ratio = sd_ratio) +
    scale_color_brewer(palette = 'Set1') -> p
  
  print(p)
  
}

# CONTRAST A ----


dds_file <- rds_f_l[1]

file_name <- strsplit(basename(dds_file), "_")[[1]][1:2]
file_name <- paste0(paste(file_name,collapse = '_'),'_PCA.png')

dds <- read_rds(dds_file)

vst <- vst(dds)

count <- assay(vst)
colData <- colData(dds) %>% as_tibble()

# data = log2(DEcount+1)
data <- count

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(colData) -> PCAdf

# color_vector <- as.character(unique(PCAdf$Design))
# n <- length(color_vector)
# getPalette <- RColorBrewer::brewer.pal(n, 'Set2')
# axis_col <- structure(getPalette, names = color_vector)

# axis_col[names(axis_col) %in% "Control"] <- "#313695"

PCAdf %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = Design), size = 5, alpha = 0.7) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # coord_fixed(ratio = sd_ratio) +
  guides(color=guide_legend("",nrow=1)) +
  see::scale_color_metro() -> ps

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"
ggsave(ps, path = out_path, filename = file_name, width = 12) 


# CONTRAST B ----

dds_file <- rds_f_l[2]

file_name <- strsplit(basename(dds_file), "_")[[1]][1:2]
file_name <- paste0(paste(file_name,collapse = '_'),'_PCA.png')

dds <- read_rds(dds_file)

vst <- vst(dds)

count <- assay(vst)
colData <- colData(dds) %>% as_tibble()

# data = log2(DEcount+1)
data <- count

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

# gcolor <- structure(c('#d53e4f', '#3288bd'), names = c('P', 'C'))

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(colData) -> PCAdf

color_vector <- as.character(unique(PCAdf$Design))
n <- length(color_vector)
getPalette <- RColorBrewer::brewer.pal(n, 'Set2')
axis_col <- structure(getPalette, names = color_vector)

axis_col[names(axis_col) %in% "Control"] <- "#313695"
  
PCAdf %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = Design), size = 5, alpha = 0.7) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # coord_fixed(ratio = sd_ratio) +
  guides(color=guide_legend("",nrow=2)) +
  see::scale_color_metro() -> ps

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"
ggsave(ps, path = out_path, filename = file_name, width = 12) 

# CONTRAST C -----

dds_file <- rds_f_l[3]

file_name <- strsplit(basename(dds_file), "_")[[1]][1:2]
file_name <- paste0(paste(file_name,collapse = '_'),'_PCA.png')

dds <- read_rds(dds_file)

vst <- vst(dds)

count <- assay(vst)
colData <- colData(dds) %>% as_tibble()

# data = log2(DEcount+1)
data <- count

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(colData) -> PCAdf

# color_vector <- as.character(unique(PCAdf$Design))
# n <- length(color_vector)
# getPalette <- RColorBrewer::brewer.pal(n, 'Set2')
# axis_col <- structure(getPalette, names = color_vector)

# axis_col[names(axis_col) %in% "Control"] <- "#313695"
  
PCAdf %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = Design), size = 5, alpha = 0.7) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # coord_fixed(ratio = sd_ratio) +
  guides(color=guide_legend("",nrow=1)) +
  see::scale_color_metro() -> ps

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"
ggsave(ps, path = out_path, filename = file_name, width = 12) 


# CONTRAST D -----

dds_file <- rds_f_l[4]

file_name <- strsplit(basename(dds_file), "_")[[1]][1:2]
file_name <- paste0(paste(file_name,collapse = '_'),'_PCA.png')

dds <- read_rds(dds_file)

vst <- vst(dds)

count <- assay(vst)
colData <- colData(dds) %>% as_tibble()

# data = log2(DEcount+1)
data <- count

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(colData) -> PCAdf

# color_vector <- as.character(unique(PCAdf$Design))
# n <- length(color_vector)
# getPalette <- RColorBrewer::brewer.pal(n, 'Set2')
# axis_col <- structure(getPalette, names = color_vector)

# axis_col[names(axis_col) %in% "Control"] <- "#313695"

PCAdf %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = Design), size = 5, alpha = 0.7) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # coord_fixed(ratio = sd_ratio) +
  guides(color=guide_legend("",nrow=1)) +
  see::scale_color_metro() -> ps

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"
ggsave(ps, path = out_path, filename = file_name, width = 12) 

# INTERACTIVE VERSION ----
# https://wilkelab.org/SDS375/slides/color-scales.html#41
# https://easystats.github.io/see/articles/seecolorscales.html#overview-of-palette-colors

out <- lapply(rds_f_l, prep_pca_df)

head(PCAdf <- do.call(rbind, out))

library(rbokeh)

PCAdf <- out[[4]]
hover <- PCAdf %>% select(LIBRARY_ID, SAMPLE_ID,DESCRIPTION, EDAD, Design)
# 
# xlabel <- paste0("PC1, VarExp: ", percentVar[1], "%")
# ylabel <- paste0("PC2, VarExp: ", percentVar[2], "%")



p <- figure(width = 1000, height = 600, legend_location = "bottom_right") %>%
  ly_points(x = PC1, y = PC2, data = PCAdf, color = Design, 
    hover = list(LIBRARY_ID, SAMPLE_ID,DESCRIPTION, EDAD, Design)) 
  # x_axis(label = xlabel) %>% y_axis(label = ylabel) %>%
  # set_palette(discrete_color = pal_color(c("red", "#313695")))

p


plot_bokeh <- function(x) {
  figure(width = 300, height = 300) %>%
    ly_points(x = PC1, y = PC2, color = Design, data = x,
      hover = list(LIBRARY_ID, SAMPLE_ID,DESCRIPTION, EDAD, Design))
}

# plot_bokeh(out[[1]])

# ?rbokeh::grid_plot()

figs <- lapply(out, plot_bokeh)


grid_plot(figs, nrow = 2, same_axes = TRUE, width = 1500, height = 700) %>%
  rbokeh::theme_legend(border_line_alpha = 0,background_fill_alpha = 0)
