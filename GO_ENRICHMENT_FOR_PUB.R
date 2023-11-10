# BOOSTRAP TOP i GENES BASED ON DIFFERENTIAL GENE EXPRESSION ANALYSIS (CONTRAST C AND D)
# Generar figura con el top 20 de genes DEGS para cada contraste:
# ¿existe un top 20 de DEG entre: a)control y metástasis, b) grados histológicos?; y b) ¿puede relacionarse esos top 20 a alguna(s) vía(s) metabólica(s) en particular?

rm(list = ls()) 

if(!is.null(dev.list())) dev.off()

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

.RES <- read_rds(paste0(path, "CONTRAST_C_AND_D_FOR_PUB.rds"))

.RES <- do.call(rbind, .RES)

RES.P <- .RES %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)

RES.P %>% group_by(sampleB) %>% dplyr::count()

# UPSET 1ST

recode_to <- c(`METASTASIS` = "(A) Metastasis", `NO METASTASIS`= "(B) No metastasis",
  `Grado I` = "(C) Stage I", `Grado II` = "(D) Stage II", `Grado III` = "(E) Stage III")

library(ggupset)

UPSETDF <- RES.P %>% 
  filter(log2FoldChange < 0 ) %>% # ONLY UP-EXPRESSED IN CANCER SAMPLES (i.e. DOWN-EXP. IN Ctrl)
  dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  group_by(transcript_id) %>%
  summarise(across(sampleB, .fns = list), n = n())

# Levels <- RES.P %>% distinct(sampleB) %>% pull()

UPSETDF %>%
  mutate(col = ifelse(n == 1, "A", "B")) %>%
  ggplot(aes(x = sampleB, fill = col)) + # , fill = SIGN
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 3.5) +
  scale_x_upset(order_by = "degree", reverse = F) +
  theme_bw(base_family = "GillSans") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans", base_size = 16) +
  axis_combmatrix(levels = recode_to) +
  labs(x = '', y = 'Number of transcripts (up-expressed)') +
  # scale_color_manual("", values = col) +
  scale_fill_manual("", values =  c("red", "black")) +
  guides(fill = guide_legend(title = "", nrow = 1)) -> p1

p1 <- p1 + theme(legend.position = "none",
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank())

ggsave(p1, filename = 'GO_ENRICHMENT_FOR_PUB_UPSET.png', 
  path = path, width = 10, height = 6, device = png, dpi = 300)


# PULL DISTINCT TRASCRIPTS

nrow(QUERY <- UPSETDF %>% filter(n == 1) %>% unnest(sampleB) ) # 1392

# sanity check, 

nrow(QUERY %>% distinct(transcript_id)) # 1392

# ANNOT

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

path <- "~/Documents/DOCTORADO/human_cancer_dataset/annot/"

orgdb <- "org.Hs.eg.db"

# semdata <- GOSemSim::godata(orgdb, ont="BP")

# write_rds(semdata, file = paste0(path, "hsGO_PB.rds"))

semdata <- read_rds(paste0(path, "hsGO_PB.rds"))

go_file <- paste0(path, 'Trinotate_report.xls.gene_ontology')

MAP <- topGO::readMappings(go_file)

STRG2GO <- data.frame(transcript_id = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  right_join(QUERY) %>%
  group_by(sampleB, transcript_id) %>%
  summarise(across(GO.ID, .fns = paste_go), n = n())

gene2GO <- split(strsplit(STRG2GO$GO.ID, ";") , STRG2GO$transcript_id)

gene2GO <- lapply(gene2GO, unlist)


# BOOSTRAPING

boostrap_enrichment <- function(STRG2GO, which_sam = "(A) Metastasis") {
  
  
  query.names <- STRG2GO %>% filter(sampleB %in% which_sam) %>% distinct(transcript_id) %>% pull()
  
  query.p <- RES.P %>% filter(transcript_id %in% query.names) %>% 
    group_by(transcript_id) %>% sample_n(1) %>%
    pull(padj, name = transcript_id)
  
  query.p <- query.p[match(query.names, names(query.p))]
  
  identical(names(query.p), query.names)
  
  allRes <- list()
  
  for (i in c(20, 50, 75, 100, length(query.p))) {
    
    
    df <- GOenrichment(query.p, query.names, gene2GO, Nodes = i, onto = "BP")
    
    allRes[[i]] <- data.frame(df, Top = i)
  }
  
  data <- do.call(rbind, allRes) %>% as_tibble() %>% mutate(sampleB = which_sam)
  
  return(data)
}

# boostrap_enrichment(STRG2GO, which_sam = "(A) Metastasis")

OUT <- lapply(recode_to, function(x) boostrap_enrichment(STRG2GO, which_sam = x))

OUT <- do.call(rbind, OUT) %>% as_tibble() 

# BIND FIVE CANCER UP-EXPRESSED

nrow(QUERY <- UPSETDF %>% filter(n == 5) %>% mutate(sampleB = "FiveCancer"))

STRG2GO <- data.frame(transcript_id = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  right_join(QUERY) %>%
  group_by(sampleB, transcript_id) %>%
  summarise(across(GO.ID, .fns = paste_go), n = n())

gene2GO <- split(strsplit(STRG2GO$GO.ID, ";") , STRG2GO$transcript_id)

gene2GO <- lapply(gene2GO, unlist)

OUTFIVE <- boostrap_enrichment(STRG2GO, which_sam = "FiveCancer")

OUT <- rbind(select_at(OUT, names(OUTFIVE)), OUTFIVE)

str(GO.IDS <- OUT %>% distinct(GO.ID) %>% pull() %>% sort())

SEM <- SEMANTIC_SEARCH(GO.IDS, orgdb = "org.Hs.eg.db",semdata = semdata)

OUT <- OUT %>% left_join(SEM, by = c("GO.ID" = "go"))

write_rds(OUT, file = paste0(path, "/GO_ENRICHMENT_FOR_PUB.rds"))

OUT <- read_rds(paste0(path, "/GO_ENRICHMENT_FOR_PUB.rds"))

# PLOT

OUT %>%
  filter(p.adj.ks < 0.05) %>%
  mutate(sampleB = factor(sampleB, levels = recode_to)) %>%
  mutate(Top = ifelse(!Top %in% c(20, 50, 75, 100), "All", Top)) %>%
  mutate(Top = factor(as.character(Top), levels = c(20, 50, 75, 100, "All"))) %>%
  ggplot(aes(x = Top, y = Term, fill = -log10(p.adj.ks))) +
  geom_tile(color = 'white', linewidth = 0.2) +
  facet_grid(parentTerm+sampleB ~., switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  # guides(fill = "none") +
  labs(y = "Biological process (Up-expressed)", x = "Top Enrichment") +
  scale_fill_viridis_c("-log10(padj)", option = "inferno", direction = -1) +
  theme(legend.position = "top",
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    strip.background.x = element_rect(fill = 'grey89', color = 'white'),
    strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
    strip.background.y = element_rect(fill = 'white', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.text.x = element_text(angle = 90, size = 10)) -> p


ggsave(p, filename = 'GO_ENRICHMENT_FOR_PUB2.png', path = path, width = 10, height = 12, device = png, dpi = 300)

# BY BARS

# https://trinkerrstuff.wordpress.com/2016/12/23/ordering-categories-within-ggplot2-facets/

# https://juliasilge.com/blog/reorder-within/

OUT %>%
  filter(Top == 50) %>%
  group_by(sampleB) %>%
  mutate(size = size / max(size)) %>%
  filter(size > 0) %>% 
  arrange(desc(size), .by_group = T) %>%
  mutate(Label = Term, row_number = row_number(Label)) %>% 
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  ggplot(aes(y = Label, x = size, color = -log10(p.adj.ks))) + # 
  facet_wrap(~ sampleB, nrow = 3, ncol = 2, scales = "free_y") +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.5) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  labs(y = "Biological process (Up-expressed)", x = "Enrichment ratio") +
  scale_color_viridis_c("-log10(padj)", option = "inferno") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 7)
    # axis.text.x = element_text(angle = 90, size = 7)
    ) -> p2

# p2
  

ggsave(p2, filename = 'GO_ENRICHMENT_FOR_PUB_BAR2.png', path = path, width = 10, height = 12, device = png, dpi = 300)

