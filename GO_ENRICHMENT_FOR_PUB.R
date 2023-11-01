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
  ggplot(aes(x = sampleB)) + # , fill = SIGN
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 3.5) +
  scale_x_upset(order_by = "freq", reverse = F) +
  theme_bw(base_family = "GillSans") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans", base_size = 16) +
  axis_combmatrix(levels = recode_to) +
  labs(x = '', y = 'Number of transcripts (up-expressed') +
  # scale_color_manual("", values = col) +
  # scale_fill_manual("", values =  col) +
  guides(fill = guide_legend(title = "", nrow = 1)) -> p1

p1 <- p1 + theme(legend.position = "top",
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank())

# PULL DISTINCT TRASCRIPTS

nrow(QUERY <- UPSETDF %>% filter(n == 1) %>% unnest(sampleB) )

# sanity check, 

QUERY %>% distinct(transcript_id)

# ANNOT

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

path <- "~/Documents/DOCTORADO/human_cancer_dataset/annot/"

go_file <- paste0(path, 'Trinotate_report.xls.gene_ontology')

MAP <- topGO::readMappings(go_file)

STRG2GO <- data.frame(transcript_id = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  right_join(QUERY) %>%
  group_by(sampleB, transcript_id) %>%
  summarise(across(GO.ID, .fns = paste_go), n = n())


STRG2GO <- split(strsplit(STRG2GO$GO.ID, ";") , STRG2GO$sampleB)

STRG2GO <- lapply(STRG2GO, unlist)

GOenrichment(p, q, STRG2GO, Nodes = 10, onto = "BP", mapping = NULL)


for (i in 1:length(STRG2GO)) {
  
  q <- STRG2GO[i]
  
  cat("\nUsing ",names(q), " Contrast...\n")
  
  
  gene2GO <- STRG2GO[names(STRG2GO) %in% names(q)][[1]]
  
  n <- length(query.go)
  
  cat("\nUsing ",n, " GO terms\n")
  
  # p <- query.p[i]
  
  query.p <- rep(0.05, n)
  
  description <- "complete topGO enrichment using split_annot"
  
  
  topGOdata <- new("topGOdata", 
    ontology = "BP", 
    description = description,
    allGenes = query.p,
    # geneSel = function(x) { x == 1 },
    geneSel = function(x) x,
    annot = annFUN.gene2GO,
    # mapping = hsGO, # omit this flag
    gene2GO = gene2GO)
  
  
  df <- GOenrichment(p, q, STRG2GO, Nodes = 10, onto = "BP", mapping = NULL)
  
  allRes[[i]] <- data.frame(df, Name = q)
}

DF1 <- do.call(rbind, allRes) %>% as_tibble() %>% mutate(facet = CONTRAST)

# STRG2GO <- data.frame(ids = rep(names(MAP[keep]),
#   sapply(MAP[keep], length)),
#   GO.ID = unlist(MAP[keep]), row.names = NULL) %>% as_tibble()
# 
