rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

# Define Color Palette: ----

colNames <- c("Both Surviving",	"Forward Only",	"Reverse Only",	"Dropped")

getPalette <- RColorBrewer::brewer.pal(length(colNames), 'Paired')

Cvalues <- structure(getPalette, names = colNames) 

path <- '~/Documents/GitHub/Cancer_sete_T_assembly/'

# 

file <- '.csv$'

file <- list.files(path, pattern = file, full.names = TRUE)

df <- lapply(file, read.table)

df <- do.call(rbind, df)

colnames(df) <- c('id', 'Input Read Pairs',colNames)

# basic plot ----

barplot(df$`Input Read Pairs`)

barplot(sort(df$`Input Read Pairs`))

#

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
caption <- ' ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36 LEADING:5 TRAILING:5'

p2 +
  labs(caption = caption) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.border = element_blank()) -> p2

# Percentage

df %>%
  mutate(across(c(3:6),
      .fns = ~./`Input Read Pairs`)) %>% 
  pivot_longer(cols = all_of(colNames)) %>%
  mutate(name = factor(name, levels = colNames)) %>%
  arrange(desc(value)) %>%
  group_by(name) %>%
  mutate(id = forcats::fct_reorder(id, `Input Read Pairs`)) %>%
  ggplot(aes(x = id, y = value*100, fill = name)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(y = '%', x = 'Sample') +
  coord_flip() +
  scale_fill_manual('', values = Cvalues) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.border = element_blank()) -> p3

library(patchwork)

p1+p2+p3
