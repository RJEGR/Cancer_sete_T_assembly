rm(list = ls()) # Limpiar la memoria de la sesion de R

options(stringsAsFactors = FALSE) # 

library(tidyverse)
library(patchwork)

# Define Color Palette: ----

colNames <- c("Both Surviving",	"Forward Only",	"Reverse Only",	"Dropped")

RColorBrewer::display.brewer.all() # Revisamos todas las paletas de colores que tenemos

getPalette <- RColorBrewer::brewer.pal(length(colNames), 'Set1')

Cvalues <- structure(getPalette, names = colNames)

path <- '~/Documents/GitHub/Cancer_sete_T_assembly/'

#Theme
my_theme <-
  ggplot2::theme_bw(base_size = 10, base_family = "GillSans") +
  ggplot2::theme(
    legend.position = 'top',
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank())

# 

file <- '.csv$'

file <- list.files(path, pattern = file, full.names = TRUE)

df <- lapply(file, read.table) # leer todas las tablas

df <- do.call(rbind, df)

colnames(df) <- c('id', 'Input Read Pairs',colNames)

head(df)

# basic plot ----

barplot(df$`Input Read Pairs`)

barplot(sort(df$`Input Read Pairs`))

# Algo publicable

df %>%
  pivot_longer(cols = all_of(colNames)) %>%
  filter(!name %in% "Input Read Pairs") %>%
  ggplot(aes(x = id, y = value, fill = name)) +
  geom_col() +
  labs(y = 'M reads', x = 'Sample') +
  coord_flip() -> p1

p1

# Intermediate ----

df %>%
  pivot_longer(cols = all_of(colNames)) %>%
  filter(!name %in% "Input Read Pairs") -> df_longer

df_longer %>%
  mutate(name = factor(name, levels = colNames)) %>%
  group_by(name) %>%
  mutate(id = forcats::fct_reorder(id, `Input Read Pairs`)) %>%
  ggplot(aes(x = id, y = value, fill = name)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(y = 'M reads', x = 'Sample') +
  coord_flip() +
  scale_fill_manual('', values = Cvalues) -> p2

p2

# Proficient

caption <- 'ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36 LEADING:5 TRAILING:5'

p2 +
  labs(caption = caption) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.border = element_blank()) -> p3

p3

# Percentage

df %>%
  mutate(across(c(3:6),
      .fns = ~./`Input Read Pairs`)) %>% 
  pivot_longer(cols = all_of(colNames)) %>%
  mutate(name = factor(name, levels = colNames)) %>%
  # arrange(desc(value)) %>%
  group_by(name) %>%
  mutate(id = forcats::fct_reorder(id, `Input Read Pairs`)) %>%
  ggplot(aes(x = id, y = value*100, fill = name)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(y = '%', x = 'Sample', 
    subtitle = caption) +
  coord_flip() +
  scale_fill_manual('', values = Cvalues) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.border = element_blank()) -> p4

p4

ggsave(p4, filename = 'trimmomatic_pct.png', path = path, 
  width = 7, height = 10)

library(patchwork)

p1+p2+p3+p4

# BUSCO DATAVZ ----

rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

path <- '~/Documents/DOCTORADO/summaries_busco/summaries_snp_trans_prev/full_tables/'

pattern_f <- 'tsv'


file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

read_tsv <- function(f) {
  
  g <- gsub(".tsv", "",basename(f))
  g <- gsub("full_table_transcripts_*", "", g)
  g <- gsub("_odb9", "", g)
  g <- str_to_title(g)
  
  df <- read.delim(f, comment.char = "#", sep = "\t", header = F)
  df %>% as_tibble(df) %>% mutate(g = g)
  #return(df)
}

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

names <- c("Complete", "Duplicated", "Fragmented", "Missing")

col <- structure(col, names = rev(names))

df <- lapply(file, read_tsv)

head(df <- do.call(rbind, df))

names(df) <- c("Busco_id",	"Status",	"Query_seq_id",	"Score",	"Length", "Db")

df %>% 
  group_by(Status, Db) %>% 
  tally() %>%
  group_by(Db) %>% mutate(pct = n / sum(n)) %>%
  mutate(Status = factor(Status, levels = rev(names))) %>%
  ggplot(aes(x = Db, y = pct, fill = Status)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "% BUSCOs", x = "") +
  coord_flip() +
  scale_fill_manual("", values = col) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = 'top')


df %>% 
  drop_na(Score) %>%
  ggplot(aes(Score, Length)) +
  geom_point(alpha = 0.5) +
  facet_grid(Status ~ Db, scales = "free") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = "BUSCO Score")


df %>% 
  drop_na(Score) %>%
  ggplot(aes(Status, Score)) +
  geom_boxplot() +
  facet_grid( ~ Db, scales = "free_y") +
  theme_bw(base_size = 14, base_family = "GillSans")

theme_set(my_theme)

df %>%
  group_by(Db) %>% 
  tally() %>%
  ggplot(aes(Db, n)) +
  geom_col() +
  labs(y = "Data set Size", x = "") +
  coord_flip() +
  geom_text(aes(label = n), size = 3, hjust = -0.05, family = "GillSans") 

df %>% 
  group_by(Status, Db) %>% 
  tally() %>%
  group_by(Db) %>% mutate(pct = n / sum(n)) %>%
  mutate(Status = factor(Status, levels = rev(names))) %>%
  ggplot(aes(x = Db, y = pct, fill = Status, group = Status)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "% BUSCOs", x = "") +
  coord_flip() +
  scale_fill_manual("", values = col) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = 'top')


# add genome index

path <- '~/Documents/DOCTORADO/summaries_busco/summaries_genome/full_tables/'
file <- list.files(path, pattern = pattern_f,  full.names = TRUE)
df2 <- lapply(file, read_tsv)
head(df2 <- do.call(rbind, df2))

names(df2) <- c("Busco_id",	"Status",	"Query_seq_id",	"Score",	"Length", "Db")


df2 %>% mutate(Index = "genome") -> df2

df %>% mutate(Index = "g_snp_trans") -> df

tbl <- rbind(df, df2)

tbl %>% 
  drop_na(Score) %>%
  ggplot(aes(Score, fill = Index)) +
  geom_histogram() +
  facet_grid(Status ~ Db, scales = "free") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = "BUSCO Score")

tbl %>% 
  drop_na(Length) %>%
  ggplot(aes(Length, fill = Index)) +
  geom_histogram() +
  facet_grid(Status ~ Db, scales = "free") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = "Sequence Length")

tbl %>% 
  drop_na(Score) %>%
  group_by(Status, Db, Index) %>% 
  tally() %>%
  filter(!Db %in% 'Bacteria') %>%
  # group_by(Db, Index) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(Status, n, fill = Index)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  facet_wrap( ~ Db, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  labs(y = "", x = "")
