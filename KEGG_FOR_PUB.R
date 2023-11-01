

# RICARDO GOMEZ-REYES
# SEARCHING TOP 20 DEGS W/ EGGNOG OR GO
# Generar figura con el top 20 de genes DEGS para cada contraste:
# ¿existe un top 20 de DEG entre: a)control y metástasis, b) grados histológicos?; y b) ¿puede relacionarse esos top 20 a alguna(s) vía(s) metabólica(s) en particular?

rm(list = ls()) 

if(!is.null(dev.list())) dev.off()

library(tidyverse)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

.DF <- read_rds(paste0(path, "CONTRAST_C_AND_D_FOR_PUB.rds"))

DF <- do.call(rbind, .DF)

# LOAD ANNOT:

path <- "~/Documents/DOCTORADO/human_cancer_dataset/annot/"

print(x <- read_rds(paste0(path, '/trinotate_eggnog.rds')))

print(cogs <- read_rds(paste0(path, '/cogs.rds')))

#download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")

egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")

names(egg) <- c("db", "nog", "proteins", "species", "class", "description")

get_eggnog <- function (x, ids, by = "transcript_id") 
{
  # trinotateR::plot_NOGs()
  
  if (by == "transcript_id") {
    
    x1 <- x %>% filter(transcript_id %in% ids)
    
    # y <- unique(x1[!is.na(eggnog), .(transcript_id, eggnog)])
    eggnog <- x1 %>% drop_na(eggnog) %>% distinct(eggnog) %>% pull(eggnog)
    
  }
  else {
    x1 <- x %>% filter(gene_id %in% ids)
    
    y <- unique(x1[!is.na(eggnog), .(gene_id, eggnog)])
    
  }
  # nogs <- gsub("(.*)\\^.*", "\\1", eggnog)
  
  nogs <- eggnog
  
  
  # y %>% separate(col = eggnog, sep = "(.*)\\^.*", into = c('nogs', 'eggnog'))
  
  n <- match(nogs, egg$nog)
  
  y <- table(unlist(strsplit(egg$class[n], "")))
  
  y <- data.frame(y)
  
  names(y) <- c('code', 'Freq')
  
  
  
  return(y)
}


# MAKE EGG/NOGG PLOTS

DF.P <- DF %>% filter(padj < 0.05 & abs(log2FoldChange) >= 2) 

DF.P %>% count(sampleB)

# UP GENES IN CANCER

UP_DF <- DF.P %>% filter(log2FoldChange <= 2)

UP_DF %>% count(sampleB)

# Separete by contrast

sam_group <- unique(UP_DF$sampleB)

out <- list()

for(i in sam_group) {
  
  j <- i
  
  UP_DF %>% filter(sampleB %in% j) %>%
    distinct(transcript_id) %>% pull(transcript_id) -> which_ids
  
  kegg_df <- get_eggnog(x, which_ids)
  out[[j]] <- kegg_df
  
}

kegg_df <- do.call(rbind, out)

# Up genes ----
kegg_df %>% as_tibble(rownames = 'group') %>%
  mutate(group = sapply(strsplit(group, "[.]"), "[", 1)) -> kegg_df_up

# Down genes ----
# DOWN IN CANCER
DOWN_DF <- DF.P %>% filter(log2FoldChange >= 2)

DOWN_DF %>% count(sampleB)

out <- list()

for(i in sam_group) {
  j <- i
  
  DOWN_DF %>% filter(sampleB %in% j) %>%
    distinct(transcript_id) %>%
    pull(transcript_id) -> which_ids
  
  kegg_df <- get_eggnog(x, which_ids)
  out[[j]] <- kegg_df
  
}

kegg_df <- do.call(rbind, out)

kegg_df %>% as_tibble(rownames = 'group') %>%
  mutate(group = sapply(strsplit(group, "[.]"), "[", 1)) -> kegg_df_down


kegg_df <- rbind(
  data.frame(kegg_df_up, lfcT = "Up-regulated"),
  data.frame(kegg_df_down, lfcT = "Down-regulated")) %>%
  mutate(lfcT = factor(lfcT, levels = c("Up-regulated", "Down-regulated")))


col_palette <- cogs %>% distinct(code, clrs)
col_palette <- structure(col_palette$clrs, names = col_palette$code)


kegg_df %>% 
  group_by(lfcT, group) %>%
  summarise(n = sum(Freq))
# ungroup() %>% summarise(N = sum(n))

DF %>% filter(padj < 0.05 & abs(log2FoldChange) >= 2) %>% 
  mutate(lfcT = sign(log2FoldChange)) %>%
  group_by(sampleB, lfcT) %>% tally()

recode_to <- c(`METASTASIS` = "(A) Metastasis", `NO METASTASIS`= "(B) No metastasis",
  `Grado I` = "(C) Stage I", `Grado II` = "(D) Stage II", `Grado III` = "(E) Stage III")

cogs %>% left_join(kegg_df) %>%
  drop_na(Freq) %>%
  # filter(Freq > 2) %>%
  filter(!grepl('unknown', name)) %>%
  group_by(group, lfcT) %>%
  arrange(desc(Freq)) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  dplyr::mutate(group = dplyr::recode_factor(group, !!!recode_to)) %>%
  mutate(group = factor(group, levels = recode_to)) %>%
  ggplot(aes(y = Freq, x = name, fill = code)) + 
  geom_col() +
  labs(x = '' , y = 'Number of transcripts') +
  # facet_wrap(lfcT ~ group, scales = 'free_y') +
  facet_grid(lfcT ~ group) +
  coord_flip() +
  geom_text(aes(label = Freq), size = 2, hjust = -0.05, family = "GillSans") +
  theme_bw(base_family = "GillSans") -> psave

psave + scale_fill_manual(values = col_palette) +
  theme(legend.position = 'none') -> psave

psave <- psave + 
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.minor = element_blank())

ggsave(psave, path = path, filename = 'EGGNOG_FOR_PUB.png', 
  width = 10, height = 6, device = png) 
#
