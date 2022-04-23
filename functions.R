runtopGO <- function(topGOdata, topNodes = 20, conservative = TRUE) {
  
  RFisher <- runTest(topGOdata, 
    algorithm = "classic", 
    statistic = "fisher")
  
  # To make this test conservative. Next we will test the enrichment using the Kolmogorov-Smirnov test. We will use the both the classic and the elim method.
  
  if(conservative) 
  {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks")
    
    RKS.elim <- runTest(topGOdata, algorithm = "elim", 
      statistic = "ks")
    
    
    
    
    allRes <- GenTable(topGOdata, 
      classicFisher = RFisher,
      classicKS = RKS, 
      elimKS = RKS.elim,
      orderBy = "elimKS", 
      ranksOf = "classicFisher", 
      topNodes = topNodes) 
  } else {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks")
    
    test.stat <- new("weightCount",
      testStatistic = GOFisherTest,
      name = "Fisher test", sigRatio = "ratio")
    
    weights <- getSigGroups(topGOdata, test.stat)
    
    allRes <- GenTable(topGOdata,
      classic = RFisher,
      KS = RKS,
      weight = weights,
      orderBy = "weight",
      ranksOf = "classic",
      topNodes = topNodes)
    
    # allRes <- GenTable(object = topGOdata, 
    #                    elimFisher = RFisher,
    #                    topNodes = topNodes)
  }
  
  return(allRes)
}

GOenrichment <- function(df, cons = T, onto = "BP", Nodes = Inf) {
  
  # df needs to include follow cols: "ids"            "sampleB"        "GO.ID" (list)      and    "padj" 
  
  df %>% pull(padj) -> query.p
  df %>% pull(ids) -> query.names
  
  names(query.p) <- query.names
  
  # keep MAP of query genes
  gene2GO <- df$GO.ID
  names(gene2GO) <- query.names
  
  description <- "complete topGO enrichment using split_annot"
  
  
  topGOdata <- new("topGOdata", 
    ontology = onto, 
    description = description,
    allGenes = query.p,
    # geneSel = function(x) { x == 1 },
    geneSel = function(x) x,
    annot = annFUN.gene2GO,
    # mapping = hsGO, # omit this flag
    gene2GO = gene2GO)
  
  # run TopGO results 
  
  allGO <- usedGO(topGOdata)
  
  if(is.infinite(Nodes)) {
    topNodes <- length(allGO)
  } else {
    topNodes <- Nodes
  }

  
  allRes <- runtopGO(topGOdata, topNodes = Nodes, conservative = cons)
  
  # make p adjustable
  
  p.adj.ks <- p.adjust(allRes$classicKS , method="BH")
  
  allRes <- cbind(allRes, p.adj.ks)
  
  allRes$Term <- gsub(" [a-z]*\\.\\.\\.$", "", allRes$Term)
  allRes$Term <- gsub("\\.\\.\\.$", "", allRes$Term)
  
  return(allRes)
  
  
}

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
