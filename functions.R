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
