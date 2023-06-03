require(org.Mm.eg.db)
#require(gage)
require(gageData)
require(tidyverse)
require(clusterProfiler)
# By @biomedical_infromatics Edris Sharif Rahmani. Jun 3, 2023
GO_dotplot <- function(data, 
                       pval, 
                       qval,
                       top) {
  up_down <- data.frame(read_csv(data))
  up_down$entrez <- mapIds(org.Mm.eg.db,
                           keys=up_down$Genes,
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
  go_enrichment <- enrichGO(up_down$entrez,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP",
                            pvalueCutoff = pval,
                            qvalueCutoff = qval)
  
  res <- go_enrichment@result
  res$total <- as.numeric(gsub("[0-9]+/", "", res$GeneRatio))
  res$GeneRatio <- res$Count/res$total
  g <- ggplot(data = res[1:top, ], aes(reorder(Description, GeneRatio, sum),
                                       GeneRatio, color = p.adjust, size = Count)) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low = "red", high = "blue") +
    xlab("Biological Process") +
    ylab("Gene Ratio") +
    theme_classic() +
    theme(axis.text = element_text(family = "Times",size = 14 , colour = "black"),
          axis.text.x = element_text(family = "Times",colour = "black", size = 10),
          axis.text.y = element_text(family = "Times",colour = "black", size = 10),
          plot.subtitle = element_text(family = "Times",size = 16, colour = "black", hjust = 0.5),
          axis.title.y = element_text(family = "Times", size = 12, angle = 90),
          axis.title.x = element_text(family = "Times", size = 12, angle = 00),
          legend.text = element_text(size = 12, family = "Times"), 
          legend.title = element_text(size = 12, family = "Times"))
  return(g)
}

GO_dotplot(data = "DESeq2/Unique_lung_drop.csv", 
           pval = 0.05,
           qval = 0.05,
           top = 10)
