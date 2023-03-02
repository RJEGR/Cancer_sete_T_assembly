

rm(list = ls()) # Limpiar la memoria de la sesion de R
if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE) # 
path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/'

count_f <- list.files(path, pattern = 'counts.matrix$',  full.names = TRUE)

mtd_f <- list.files(path, pattern = 'metadata.tsv',  full.names = TRUE)

mtd <- readr::read_tsv(mtd_f) %>% mutate_all(list(~ str_replace(., "SIN_CANCER", "Control")))

dim(raw_count <- read.delim(count_f, sep = "\t", header = T, row.names = 1))

keep <- rowSums(edgeR::cpm(raw_count) > 1) >= 2

sum(keep) # N transcripts

nrow(count <- raw_count[keep,])

count <- round(count)

# Prevalence of features ----

apply(raw_count, 1, function(x) sum(x > 0)) %>% table()

prevelancedf = apply(raw_count, 1, function(x) sum(x > 0))

mean_se = apply(raw_count, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(raw_count),
  mean_se) %>% 
  as_tibble(rownames = "id") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf


query.ids <- prevelancedf %>% filter(Prevalence == 60) %>% pull(id)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

x <- read_rds(paste0(path, '/trinotate_eggnog.rds'))
# names(x) <- c("gene_id", "transcript_id", "protein", "eggnog")

cogs <- read_rds(paste0(path, '/cogs.rds'))

download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")

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

kegg_df <- get_eggnog(x = x, ids = query.ids, by = "transcript_id")


col_palette <- cogs %>% distinct(code, clrs)
col_palette <- structure(col_palette$clrs, names = col_palette$code)



cogs %>% left_join(kegg_df) %>%
  drop_na(Freq) %>%
  filter(Freq > 2) %>%
  # filter(!grepl('unknown', name)) %>%
  arrange(desc(Freq)) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  ggplot(aes(y = Freq, x = name, fill = code)) + 
  geom_col() +
  labs(x = '' , y = 'Number of transcripts') +
  coord_flip() +
  geom_text(aes(label = Freq), size = 2, hjust = -0.05, family = "GillSans") +
  theme_bw(base_family = "GillSans") -> psave

psave + scale_fill_manual(values = col_palette) +
  theme(legend.position = 'none') -> psave

par(family = "GillSans")

out_path <- "~/Documents/DOCTORADO/human_cancer_dataset/"
ggsave(psave, path = out_path, filename = 'eggnog_prevalence_transcripts.png', 
  width = 10, height = 6) 

# Add aditional features as contrast ----

df <- xlsx::read.xlsx("~/Documents/DOCTORADO/human_cancer_dataset/CONTRAST_A_Annot_count_down_up_genes.xlsx",sheetIndex = "RESULTS") # tarda muchisimo

sam_group <- unique(df$sampleB)


df %>% as_tibble()

out <- list()

for(i in sam_group) {
  j <- i
  
  up_df %>% filter(sampleB %in% j) %>%
    distinct(ids) %>%
    pull(ids) -> which_ids
  
  kegg_df <- get_eggnog(x, which_ids)
  out[[j]] <- kegg_df
  
}

str(out)

kegg_df <- do.call(rbind, out)

# Up genes ----
kegg_df %>% as_tibble(rownames = 'group') %>%
  mutate(group = sapply(strsplit(group, "[.]"), "[", 1)) -> kegg_df_up

# Down genes ----

out <- list()

for(i in sam_group) {
  j <- i
  
  down_df %>% filter(sampleB %in% j) %>%
    distinct(ids) %>%
    pull(ids) -> which_ids
  
  kegg_df <- get_eggnog(x, which_ids)
  out[[j]] <- kegg_df
  
}

kegg_df <- do.call(rbind, out)

kegg_df %>% as_tibble(rownames = 'group') %>%
  mutate(group = sapply(strsplit(group, "[.]"), "[", 1)) -> kegg_df_down


kegg_df <- rbind(
  data.frame(kegg_df_up, lfcT = up),
  data.frame(kegg_df_down, lfcT = down))


col_palette <- cogs %>% distinct(code, clrs)
col_palette <- structure(col_palette$clrs, names = col_palette$code)


kegg_df %>% group_by(lfcT, group) %>%
  summarise(n = sum(Freq))
# ungroup() %>% summarise(N = sum(n))

res.p %>% group_by(sampleB, lfcT) %>% tally()

cogs %>% left_join(kegg_df) %>%
  drop_na(Freq) %>%
  filter(Freq > 2) %>%
  filter(!grepl('unknown', name)) %>%
  group_by(group, lfcT) %>%
  arrange(desc(Freq)) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
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

ggsave(psave, path = path, filename = 'eggnog_up_down_sign.png', 
  width = 10, height = 6) 