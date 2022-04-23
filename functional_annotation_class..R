
# Ricardo Gomez-Reyes April 2022
# Transcriptomic Cancer analysis

# Functional annotation enrichment for identify a common modular structure between G.O terms hierarchy

# In many cases the list of differentially expressed genes is not sufficient for accurate inference of the underlying biology. Additional biological knowledge needs to be included to enhance the interpretation of such a list of genes. With the development of biological knowledge databases, information for augmenting gene expression data is available. Biologically interesting sets of genes, for example genes that belong to a pathway or genes known to have the same biological function, can now be compiled. A popular choice for gene sets are genes collected under Gene Ontology (GO) terms (Alexa et al, 2006)

# read full methods and models here: https://academic.oup.com/bioinformatics/article/22/13/1600/193669#1877071

# separate exclusive (distinct) differential expressed transcripts found in the multiple cancer comparison vs control

# up-regulated means transcripts well expressed in cancer samples
# down-regulated means transcripts well expressed in control or down-expressed in cancer

# 1) Evaluate a descriptive  functional annotation of Diff Exp transcripts
# 2) run the gene enrichment ontology

rm(list = ls()) # Limpiar la memoria de la sesion de R

options(stringsAsFactors = FALSE) # 

library(tidyverse)
library(topGO)
library(ggupset)
library(patchwork)

# source("/Users/cigom/Documents/GitHub/Cancer_sete_T_assembly/functions.R")

url <- 'https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R'

source(url)

path <- '~/Documents/DOCTORADO/human_cancer_dataset/annot/'

load(paste0(path, '/count_annot_multiple_contrast_vs_control_up_down_genes.Rdata'))

go_file <- paste0(path, 'Trinotate_report.xls.gene_ontology')

MAP <- topGO::readMappings(go_file)

swissdf <- read_tsv(paste0(path, '/Trinotate.xls.blastx.tsv'))

# 1) Separate exclusive transcripts ----

prep_dist_data <- function(df, select_exclusive = T) {
  
  # names(df) should contain at least cols:
  # "ids"            
  # "sampleB"       
  #  "log2FoldChange" 
  # "pvalue" or    "padj" 
  # "lfcT"
  

  if(select_exclusive) {
    
    # Select only n_sam == 1 ----
    
    df %>% # this is the only input needed
      dplyr::select(ids, sampleB) %>% 
      group_by(ids) %>%
      summarise(n_sam = length(sampleB), 
        across(sampleB, .fns = list)) -> distinct_trans
    
    # Sanity check of correct grouping exclusive ids
    # Also you can compare against the upset png value bars
    
    distinct_trans %>% group_by(n_sam) %>% tally() -> prevalence
    
    distinct_trans %>% 
      filter(n_sam == 1) %>% 
      mutate(sampleB = unlist(sampleB)) -> distinct_trans
    
    cat("\n Prevalence: ",
      paste0(prevalence$n_sam, 's'), "\n", 
      "\t",prevalence$n, "\n")
    
  } else {
    df %>% dplyr::select(ids, sampleB) %>% 
      group_by(sampleB) %>% distinct(ids) -> distinct_trans
    
    distinct_trans %>% group_by(sampleB) %>% tally() -> prevalence
    
    cat("\n Prevalence: ",
      paste0(prevalence$sampleB, 's'), "\n", 
      "\t",prevalence$n, "\n")
    
  }
 

  
  
  
  # 2) Subset of GO ids based on the distinct transcript dataset ----
  
  distinct_trans %>% pull(ids) %>% unique(.) -> query.ids
  
  n <- sum(keep <- names(MAP) %in% query.ids) / length(query.ids) 
  
  cat('\n % of ids mapped in Gen Ontology: ', n*100, '%\n')
  
  STRG2GO <- data.frame(ids = rep(names(MAP[keep]),
    sapply(MAP[keep], length)),
    GO.ID = unlist(MAP[keep]), row.names = NULL) %>% as_tibble()
  
  # Sanity check: % of recovery querys with go ids
  # sum(unique(STRG2GO$ids) %in% query.ids) / length(query.ids) 
  
  distinct_trans %>% 
    left_join(STRG2GO, by = "ids") %>%
    group_by(ids, sampleB) %>%
    summarise(across(GO.ID, .fns = list), .groups = "drop_last") %>%
    ungroup() -> distinct_trans
  
  # Add p values to the distinct trans
  distinct_trans %>% 
    left_join(df %>% dplyr::select(ids, padj, log2FoldChange), 
      by = "ids") -> distinct_trans
  
  
  return(distinct_trans)
  
}

distinct_trans_up <- prep_dist_data(up_df)

distinct_trans_down <- prep_dist_data(down_df)

# Sanity check: how much exclusive transcripts w/ go id ?

distinct_trans_up %>% group_by(sampleB) %>% tally()

distinct_trans_down %>% group_by(sampleB) %>% tally()

# 3) run go enrichment ----

sam_group <- unique(distinct_trans_up$sampleB)

# df <- distinct_trans %>% filter(sampleB %in% sam_group[1])
# allRes <- GOenrichment(df, Nodes = 20)

outup <- list()
# up-regulated means transcripts well expressed in cancer samples

for(i in sam_group) {
  
  j <- i
  
  df <- distinct_trans_up %>% filter(sampleB %in% j)
  
  allRes <- GOenrichment(df, Nodes = 20)
  
  allRes <- allRes %>% mutate(sampleB = j) %>% as_tibble()
  outup[[j]] <- allRes
  
}

# down-regulated means transcripts well expressed in control or down-expressed in cancer

outdown <- list()

# up

for(i in sam_group) {
  
  j <- i
  
  df <- distinct_trans_up %>% filter(sampleB %in% j)
  
  allRes <- GOenrichment(df, Nodes = 20)
  
  allRes <- allRes %>% mutate(sampleB = j) %>% as_tibble()
  outdown[[j]] <- allRes
  
}


do.call(rbind, outup) -> allRes_up

do.call(rbind, outdown) -> allRes_down

write.table(allRes_up, file = paste0(path, 'multiple_contrast_control_vs_cancer_up_topgo.csv'))

write.table(allRes_down, file = paste0(path, 'multiple_contrast_control_vs_cancer_down_topgo.csv'))

distinct_trans_up %>% group_by(sampleB) %>% tally()
distinct_trans_down %>% group_by(sampleB) %>% tally()


save(allRes_up, allRes_down, distinct_trans_up, distinct_trans_down,
  file = paste0(path, '/multiple_contrast_vs_control_up_down_topgo.Rdata'))


save(x, up_df, down_df, data, prevelancedf, mtd, file = paste0(path, '/count_annot_multiple_contrast_vs_control_up_down_genes.Rdata'))

# nrow(topGOdown.p <- allRes_down[which(as.numeric(allRes_down$classicKS) < 0.05),])


pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

# image(volcano, col = pal)

allRes_up %>% 
  filter(!sampleB %in% 'CON_CANCER') %>%
  group_by(sampleB) %>%
  arrange(desc(Annotated)) %>%
  mutate(Term = factor(Term, levels = unique(Term))) %>%
  # mutate(Annotated = Annotated / max(Annotated)) %>%
  ggplot() +
  geom_col(aes(x = Term, y = Annotated)) + # fill = -log10(p.adj.ks)
  coord_flip() +
  facet_wrap(~ sampleB, scales = 'free') +
  theme_bw(base_family = "GillSans") +
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1")) -> psave
  # scale_fill_gradientn(colours = pal)
  
# due to down regulates are up in cancer, lets reduce terms here!

# Term reduction by

library(GOSemSim)
library(rrvgo)

str(goid1 <- sort(unique(allRes_down$GO.ID))) # using only significant genes


scores1 <- -log(as.numeric(allRes_down$classicKS))

names(scores1) <- allRes_down$GO.ID

hsGO <- read_rds(paste0(path, '/hsGO_BP.rds'))


# SimMatrix <- calculateSimMatrix(goid, orgdb = 'org.Hs.eg.db', ont="BP", method = 'Wang')

SimMatrix1 <- GOSemSim::termSim(goid1, goid1, semData = hsGO,  method = "Wang")

# SimMatrix2 <- GOSemSim::termSim(goid2, goid2, semData = hsGO,  method = "Wang")

reducedTerms1 <- reduceSimMatrix(SimMatrix1, scores1, threshold = 0.9, orgdb = 'org.Hs.eg.db')

# and plot

reducedTerms1 %>%
  arrange(desc(size)) %>%
  mutate(parentTerm = factor(parentTerm, levels = unique(parentTerm))) %>% 
  ggplot() +
  geom_col(aes(x = parentTerm, y = size)) +
  coord_flip() +
  theme_bw(base_family = "GillSans") -> psave2

psave2
psave

# Network of interactions ----
# Test

library(tidygraph)
library(igraph)
library(ggraph)

allRes_up

distinct_trans_up %>% 
  filter(sampleB %in% 'Cancer_IDC') %>%
  arrange(log2FoldChange) %>%
  head(10) -> net_input



net <- data.frame(ids = rep(net_input$ids,
  sapply(net_input$GO.ID, length)),
  GO.ID = unlist(net_input$GO.ID), row.names = NULL) %>% as_tibble()

net %>%
  left_join(allRes_up %>% dplyr::select(GO.ID,Term, Annotated))  -> net

g <- as_tbl_graph(net, directed = T)

g %>%
  activate(nodes) %>%
  left_join(net_input, by = c('name'='ids')) %>%
  mutate(
    betweenness = betweenness(.), 
    degree = centrality_degree(),
    centrality = components(.)$membership,
    transitivity = transitivity(.)) -> g

# La centralidad de grado («degree centrality»)
# La cercanía («closeness»)
# La intermediación («betweenness»)

hist(V(g)$centrality)
hist(V(g)$degree)
hist(V(g)$betweenness)

# group <- igraph::cluster_edge_betweenness(g)$membership

# plot(net)

g %>% 

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')

ggraph(layout) +
  ggforce::geom_mark_hull(
    aes(x, y, group = centrality),
    color = NA,
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25) -> gplot


gplot +
  geom_edge_link(aes(color = Term))+
  geom_node_point() + #aes(color = abs(log2FoldChange), size = degree)) + # , alpha = GS
  geom_node_text(aes(label = name), repel = TRUE, size = 2) +
  scale_color_viridis_c(name = 'log2FC', option = 'magma', direction = -1) +
  # scale_color_continuous(name = 'Gene Correlation') + # values = nodeCol
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top")+
  coord_fixed()
