
# Relacionar resultados: abundancia + annotation + DE
# EggNOG: A database of orthology relationships, functional annotation, and gene evolutionary histories.


# Load libraries and functions -----

library(tidyverse)

# Cargar resultados:
# Resultados de DE ----

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp'
file_name <- paste0(path, '/multiple_contrast_vs_control_res.rds')
res <- read_rds(file_name) # Para leer archivos tipo rds (Rdata) desde R. El formato rds es una version de archivos compridos de R

res %>% distinct(sampleA, sampleB)

res %>% arrange(log2FoldChange)

# res %>% sample_n(1000) %>%
#   ggplot() +
#   geom_point(aes(y = log2FoldChange, x = baseMean, color = -log10(padj)))


alpha = 0.05 

lfcThreshold = 2

# Queremos genes que fueron superiores al umbral lfcThreshold & < alpha == FC+sig

res %>% 
  mutate(cc = NA) %>% 
  mutate(cc = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, 'sigfc', cc)) %>%
  drop_na(cc) -> res.p

up <- 'Up-regulated' # transcript w/ log2FoldChange > lfcThreshold
down <- 'Down-regulated' # transcript w/ log2FoldChange < lfcThreshold

res.p %>%
  mutate(lfcT = 'Basal') %>%
  mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, down, lfcT)) %>%
  mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, up, lfcT)) -> res.p

# res.p %>% group_by(sampleB, lfcT) %>% tally()

up_df <- res.p %>% filter(grepl(up, lfcT))
down_df <- res.p %>% filter(grepl(down, lfcT))

# Matriz de abundacia ----

path <- '~/Documents/DOCTORADO/human_cancer_dataset/DiffExp'

pattern_f <- 'counts.matrix$'

file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

read_tsv <- function(f) {
  g <- basename(f)
  g <- sapply(strsplit(g, "\\."), `[`, 2)
  g <- str_to_title(g)
  
  df <- read.delim(f, sep = "\t", header = T, row.names = 1)
  # df %>% as_tibble(df, rownames = 'id') %>% mutate(g = g)
  return(df)
}

mtd <- read.csv(file = paste0(path, '/metadata.csv'))

rownames(mtd) <- mtd$sample_id


df <- lapply(file, read_tsv)

head(df <- do.call(rbind, df))

count_raw <- df # %>% filter(g %in% 'Isoform') %>% select(-g)

nrow(count_raw)

# Filter data by removing low-abundance genes

keep <- rowSums(edgeR::cpm(count_raw) > 1) >= 2


nrow(data <- count_raw[keep,])

# Anotacion ----



path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

# trino_file <- list.files(path = path, pattern = 'xls$', full.names = T)
# x <- data.table::fread(trino_file, sep = '\t', na.strings = '.')
# names(x)[1] <- 'gene_id'
# 

x <- read_rds(paste0(path, '/trinotate_eggnog.rds'))

cogs <- read_rds(paste0(path, '/cogs.rds'))


download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")

egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")

names(egg) <- c("db", "nog", "proteins", "species", "class", "description")

str(egg)


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

ids <- res.p %>% distinct(ids) %>% pull(ids)

kegg_df <- get_eggnog(x = x, ids = ids)

cogs %>% left_join(kegg_df) %>% 
  arrange(desc(Freq)) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  ggplot(aes(x = name, y = Freq)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = Freq), size = 3, hjust = -0.05, family = "GillSans") 

# Separete by contrast

sam_group <- unique(res.p$sampleB)

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
#


# (omit)
# x %>% mutate(eggnog = gsub("(.*)\\^.*", "\\1", eggnog)) %>% select(transcript_id, eggnog) -> x_clean


# write_rds(x_clean, file = paste0(path, '/trinotate_eggnog.rds'))
# write_rds(x_clean, file = paste0(path, '/trinotate_eggnog.rds'))
# data(cogs, package = 'trinotateR')
# write_rds(cogs, file = paste0(path, '/cogs.rds'))


x %>% drop_na(eggnog) %>%
  left_join(egg, by = c('eggnog', 'nog'))
  
# # count freq of genes in count data

# Prevalence of features ----

prevelancedf = apply(data, 1, function(x) sum(x > 0))

PrevalenceL <- sort(unique(prevelancedf))

mean_se = apply(data, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

identical(names(prevelancedf), rownames(mean_se))

TotalAbundance <- rowSums(data)

identical(names(TotalAbundance), rownames(mean_se))

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = TotalAbundance,
  mean_se) %>% 
  as_tibble(rownames = "id") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf %>% 
  group_by(Prevalence = as.character(Prevalence)) %>%
  tally() %>% 
  mutate(Prevalence = factor(Prevalence, levels = PrevalenceL)) %>%
  ggplot() +
  geom_col(aes(x = Prevalence, y = n)) +
  labs(y = 'Number of transcripts')
  

# prevelancedf %>% 
#   arrange(Prevalence) %>%
#   mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
#   mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence))) %>%
#   mutate(TotalAbundance = log2(TotalAbundance+1)) %>%
#   ggplot(aes(TotalAbundance)) + geom_histogram() + 
#   facet_wrap(~ Prevalence, scales = 'free_y') -> p1
# 
# dat_text <- prevelancedf %>% group_by(Prevalence) %>% tally() %>% 
#   mutate(cumsum = cumsum(n)) %>% arrange(Prevalence) %>%
#   mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
#   mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))
# 
# p1 + geom_text(
#   data    = dat_text, family = "GillSans",
#   mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " genes")),
#   hjust   = -1, vjust   = -2, label.size = 0.2) + 
#   theme_classic(base_size = 7, base_family = "GillSans") +
#   labs(x = expression(~Log[2]~('TotalAbundance'~+1)), y = "")


# Prepare data for Functional annot ----

save(x, up_df, down_df, data, prevelancedf, mtd, file = paste0(path, '/count_annot_multiple_contrast_vs_control_up_down_genes.Rdata'))


# #Experiment w  complex heatmap 

# devtools::install_github("karthik/wesanderson")
# https://wesandersonpalettes.tumblr.com/

# Gradient color

pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

# image(volcano, col = pal)

dds <- read_rds("~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/Cancer_vs_Control_dds.rds")

m <- DESeq2::counts(dds, normalized=T)

rbind(up_df, down_df) %>% distinct(ids) %>% pull() -> query.ids

dim(m <- m[rownames(m) %in% query.ids,])

dist.method <- 'euclidean'

linkage.method <- 'complete'

m <- log2(m+1)

hc_samples <- hclust(dist(t(m), method = dist.method), 
  method = linkage.method)

hc_sam_ord <- hc_samples$labels[hc_samples$order]

hc_genes <- hclust(dist(m, method = dist.method), 
  method = linkage.method)

hc_genes_ord <- hc_genes$labels[hc_genes$order]


# data[rownames(data) %in% query.ids ,] %>% 

m %>% 
  as_tibble(rownames = 'ids') %>%
  pivot_longer(cols = names(data), names_to = 'sample_id') %>%
  left_join(mtd) %>%
  mutate(ids = factor(ids ,levels = hc_genes_ord)) %>%
  mutate(sample_id = factor(sample_id, levels = hc_sam_ord)) %>%
  mutate(value = log2(value+1)) %>%
  ggplot(aes(x=sample_id, y = ids, fill = value)) +
  geom_raster() +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  theme_classic() +
  theme(axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 90, 
      hjust = 1, vjust = 1, size = 7)) +
  scale_fill_gradientn('Log2(x+1)', colours = pal) +
  ggh4x::facet_grid2(~g2, scales = 'free_x', space = 'free_x')

# upset ----
# library(ggupset)
# 

# devtools::install_github("const-ae/ggupset")

# test
down_df %>% distinct(sampleB) %>% pull() -> sets

down_df %>% dplyr::select(ids, sampleB, log2FoldChange) %>% 
  mutate(log2FoldChange = 1) %>%
  pivot_wider(names_from = sampleB, values_from = log2FoldChange, values_fill = 0) %>% 
  dplyr::select(-ids) %>%
  as.data.frame() %>%
  UpSetR::upset(sets = sets, keep.order = TRUE, order.by = 'degree')


rbind(up_df, down_df) %>% 
  dplyr::select(ids, sampleB, lfcT) %>%
  group_by(ids, lfcT) %>%
  # distinct() %>% # Esto no resulta tan necesario ya que necesitamos trabajar con todo los ids
  summarise(across(sampleB, .fns = list)) %>% # Super!!!!
  ggplot(aes(x=sampleB, fill = lfcT, color = lfcT)) +
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 1.5) +
  scale_x_upset(order_by = "degree") +
  theme_classic(base_family = "GillSans") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans") +
  labs(x = '', y = 'Number of transcripts') +
  scale_fill_manual(name = "", values = c('#1f78b4', '#a6cee3')) +
  scale_color_manual(name = "", values = c('#1f78b4', '#a6cee3')) +
  theme(legend.position = 'top') -> psave

ggsave(psave, filename = 'multiple_contrast_control_vs_cancer_up-down_upset.png', path = path,
  width = 8,height = 3.7)
