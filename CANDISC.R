# CANONICAL DISCRIMINATORY ANALYSIS TO "DRAW" CORRELATION BETWEEN LIBRARY SOURCE AND GENE EXPRESSION:
# ENCONTRAR UN PERFIL DE EXPRESION SIMILAR ENTRE CONJUNTO DE MUESTRAS:

set.seed(202309)

library(tidyverse)
library(candisc)


path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

rds_f_l <- list.files(path, pattern = 'CONTRAST_',  full.names = TRUE)

rds_f <- rds_f_l[1]

dds <- read_rds(rds_f)

# 1)
M <- counts(dds, normalized = T)[query.genes, ]
M <- DESeq2::varianceStabilizingTransformation(round(M))

# # 2)
# 
# vst <- DESeq2::vst(dds) # vst if cols > 10
# M <- assay(vst)
# keep <- rownames(M) %in% query.genes
# head(M <- M[keep,])

colData <- colData(dds) %>% as_tibble()

# A) ====

DESIGN <- "Design"

data <- M %>% as_tibble(rownames = "transcript_id") %>% 
  pivot_longer(-transcript_id, values_to = 'Reads', names_to = "LIBRARY_ID") %>%
  left_join(colData) %>%
  mutate_if(is.character, as.factor) 
# group_by_at(c(DESIGN, "transcript_id")) %>%
# summarise(Reads = sum(Reads))

table(colData$Design)

dat <- data %>%
  pivot_wider(names_from = all_of(DESIGN), values_from = Reads, values_fill = 0)

mod <- lm(cbind(`Control`,`CON_CANCER`) ~ LIBRARY_ID , data=dat)

term <- "LIBRARY_ID"

can <- candisc(mod, data = dat, term = term, ndim = "2")

heplot(candisc(mod, data = dat, term = term,  ndim = "1"))

# plot(can)

ellipse_df <- can$means %>%
  data.frame() %>%
  rownames_to_column(var= term) %>%
  left_join(colData, by = term)

data.frame(can$structure+can$coeffs.std ) %>%
  rownames_to_column(var= 'Metric') -> segment_df

can$scores %>% as_tibble() %>%
  mutate(g = DESIGN) %>%
  left_join(colData, by = term) -> scores_df

x_label <- paste0("Can1 (", round(can$pct[1], 1), "%)")
y_label <- paste0("Can2 (", round(can$pct[2], 1), "%)")

ggplot() + 
  geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  geom_point(data = scores_df,
    aes_string(x = "Can1",y = "Can2", color = DESIGN),  alpha = 0.5) +
  ggforce::geom_mark_circle(
    data = ellipse_df, aes_string(x = "Can1", y = "Can2", color = DESIGN, fill = DESIGN)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  geom_segment(data = segment_df, 
    aes(x = 0, y = 0, xend = Can1, yend = Can2),
    arrow = arrow(length = unit(0.07, "cm"))) +
  geom_text(data = segment_df, 
    aes(x = Can1+sign(Can1)*0.2, 
      y = Can2+sign(Can1)*0.2, label = Metric), 
    family = "GillSans") +
  coord_fixed() +
  labs(x = x_label, y = y_label) +
  theme(panel.border = element_blank(), legend.position = 'none') -> p

p + 
  facet_grid(~ g) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey50', color = 'white')) +
  xlim(c(-5, 5)) +
  ylim(c(-2.5, 2.5)) +
  see::scale_fill_see(reverse = F) +
  see::scale_color_see(reverse = F)  -> p1

# B) ====

DESIGN <- "CONTRASTE_B"

table(colData$CONTRASTE_B)

dat <- data %>%
  pivot_wider(names_from = all_of(DESIGN), values_from = Reads, values_fill = 0)

mod <- lm(cbind(`Control`,`CARCINOMA`, `CARCINOMA DUCTAL`, `CARCINOMA DUCTAL INFILTRANTE`, `CARCINOMA INFILTRANTE`, `CARCINOMA INTRADUCTAL`, `CARCINOMA LOBULILLAR INFILTRANTE`, `CARCINOMA MEDULAR`, `CARCINOMA PAPILAR`) ~ LIBRARY_ID , data=dat)

term <- "LIBRARY_ID"

can <- candisc(mod, data = dat, term = term, ndim = "2")

heplot(candisc(mod, data = dat, term = term,  ndim = "1"))

# plot(can)

ellipse_df <- can$means %>%
  data.frame() %>%
  rownames_to_column(var= term) %>%
  left_join(colData, by = term)

data.frame(can$structure+can$coeffs.std ) %>%
  rownames_to_column(var= 'Metric') -> segment_df

can$scores %>% as_tibble() %>%
  mutate(g = DESIGN) %>%
  left_join(colData, by = term) -> scores_df

x_label <- paste0("Can1 (", round(can$pct[1], 1), "%)")
y_label <- paste0("Can2 (", round(can$pct[2], 1), "%)")

ggplot() + 
  geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  geom_point(data = scores_df,
    aes_string(x = "Can1",y = "Can2", color = DESIGN),  alpha = 0.5) +
  ggforce::geom_mark_circle(
    data = ellipse_df, aes_string(x = "Can1", y = "Can2", color = DESIGN, fill = DESIGN)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  geom_segment(data = segment_df, 
    aes(x = 0, y = 0, xend = Can1, yend = Can2),
    arrow = arrow(length = unit(0.07, "cm"))) +
  ggrepel::geom_text_repel(data = segment_df, 
    aes(x = Can1+sign(Can1)*0.2, 
      y = Can2+sign(Can1)*0.2, label = Metric), 
    family = "GillSans", size = 3, max.overlaps = 20) +
  coord_fixed() +
  labs(x = x_label, y = y_label) +
  theme(panel.border = element_blank(), legend.position = 'none') -> p

p + 
  facet_grid(~ g) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey50', color = 'white')) +
  xlim(c(-5, 5)) +
  ylim(c(-2.5, 2.5)) +
  see::scale_fill_see(reverse = F) +
  see::scale_color_see(reverse = F) -> p2

p2

# C) ====

DESIGN <- "CONTRASTE_C"

table(colData$CONTRASTE_C)

dat <- data %>%
  pivot_wider(names_from = all_of(DESIGN), values_from = Reads, values_fill = 0)


mod <- lm(cbind(`Control`,`Grado I`, `Grado II`, `Grado III`, `Indiferenciado`) ~ LIBRARY_ID , data=dat)

term <- "LIBRARY_ID"
# term <- "protein_name"

can <- candisc(mod, data = dat, term = term, ndim = "2")

heplot(candisc(mod, data = dat, term = term,  ndim = "1"))

# plot(can)

ellipse_df <- can$means %>%
  data.frame() %>%
  rownames_to_column(var= term) %>%
  left_join(colData, by = term)

data.frame(can$structure+can$coeffs.std ) %>%
  rownames_to_column(var= 'Metric') -> segment_df

can$scores %>% as_tibble() %>%
  mutate(g = DESIGN) %>%
  left_join(colData, by = term) -> scores_df

x_label <- paste0("Can1 (", round(can$pct[1], 1), "%)")
y_label <- paste0("Can2 (", round(can$pct[2], 1), "%)")

ggplot() + 
  geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  geom_point(data = scores_df,
    aes_string(x = "Can1",y = "Can2", color = DESIGN),  alpha = 0.5) +
  ggforce::geom_mark_circle(
    data = ellipse_df, aes_string(x = "Can1", y = "Can2", color = DESIGN, fill = DESIGN)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  geom_segment(data = segment_df, 
    aes(x = 0, y = 0, xend = Can1, yend = Can2),
    arrow = arrow(length = unit(0.07, "cm"))) +
  ggrepel::geom_text_repel(data = segment_df, 
    aes(x = Can1+sign(Can1)*0.2, 
      y = Can2+sign(Can1)*0.2, label = Metric), 
    family = "GillSans") +
  coord_fixed() +
  labs(x = x_label, y = y_label) +
  theme(panel.border = element_blank(), legend.position = 'none') -> p

p + 
  facet_grid(~ g) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey50', color = 'white')) +
  # xlim(c(-5, 5)) +
  # ylim(c(-2.5, 2.5)) +
  see::scale_fill_see(reverse = F) +
  see::scale_color_see(reverse = F) -> p3


# D) ====

DESIGN <- "CONTRASTE_D"

table(colData$CONTRASTE_D)


dat <- data %>%
  pivot_wider(names_from = all_of(DESIGN), values_from = Reads, values_fill = 0)


mod <- lm(cbind(`Control`,`METASTASIS`, `NO METASTASIS`) ~ LIBRARY_ID , data=dat)

can <- candisc(mod, data = dat, term = term, ndim = "2")

heplot(candisc(mod, data = dat, term = term,  ndim = "1"))

# plot(can)

ellipse_df <- can$means %>%
  data.frame() %>%
  rownames_to_column(var= term) %>%
  left_join(colData, by = term)

data.frame(can$structure+can$coeffs.std ) %>%
  rownames_to_column(var= 'Metric') -> segment_df

can$scores %>% as_tibble() %>%
  mutate(g = DESIGN) %>%
  left_join(colData, by = term) -> scores_df

x_label <- paste0("Can1 (", round(can$pct[1], 1), "%)")
y_label <- paste0("Can2 (", round(can$pct[2], 1), "%)")

ggplot() + 
  geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  geom_point(data = scores_df,
    aes_string(x = "Can1",y = "Can2", color = DESIGN),  alpha = 0.5) +
  ggforce::geom_mark_circle(
    data = ellipse_df, aes_string(x = "Can1", y = "Can2", color = DESIGN, fill = DESIGN)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  geom_segment(data = segment_df, 
    aes(x = 0, y = 0, xend = Can1, yend = Can2),
    arrow = arrow(length = unit(0.07, "cm"))) +
  geom_text(data = segment_df, 
    aes(x = Can1+sign(Can1)*0.2, 
      y = Can2+sign(Can1)*0.2, label = Metric), 
    family = "GillSans") +
  coord_fixed() +
  labs(x = x_label, y = y_label) +
  theme(panel.border = element_blank(), legend.position = 'none') -> p

p + 
  facet_grid(~ g) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey50', color = 'white')) +
  # xlim(c(-5, 5)) +
  # ylim(c(-2.5, 2.5)) +
  see::scale_fill_see(reverse = F) +
  see::scale_color_see(reverse = F) -> p4




library(patchwork)

p1 + p2 + p3 + p4 + plot_layout(nrow = 2)
