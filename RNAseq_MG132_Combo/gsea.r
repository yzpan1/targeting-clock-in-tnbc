library(org.Hs.eg.db)
library(clusterProfiler)
library(GSEABase)
library(fgsea)
library(tibble)
library(dplyr)
library(ggplot2)


#################   Prepare gmt files and gene set objects #################
# gmtfile = 'h.all.v7.4.symbols.gmt'
# pathways.hallmark <- gmtPathways(gmtfile)
cacgtg <- gmtPathways('./gmt/CACGTG_MYC_Q2.gmt')
circadianent <- gmtPathways("./gmt/GOBP_ENTRAINMENT_OF_CIRCADIAN_CLOCK.gmt")
hallmarks <- gmtPathways("./gmt/h.all.v7.4.symbols.gmt")
hiftarget <- gmtPathways("./gmt/c6.all.v7.4.symbols.gmt")

#################  Prepare ranked list objects  #################
x <- read.csv("/home/kaylab/Documents/projects/20220923_MB231_MG132/shpDESigByLFC.csv")
grpgsea <- x[,c("Symbol", "log2FoldChange")]
  # The following line can be used to get rid of NA when generating ranked list with p-value, etc
#grpgsea <- grpgsea[!is.na(shpgsea$padj),]
grpgsea
#write.table(shpgsea, file = "shpDEByP.rnk", quote = F, sep = "\t", row.names = FALSE)
ranks <- deframe(grpgsea)
grpGSEAres <- fgsea(pathways = hallmarks, stats = ranks, nperm = 1000)

grpGSEAresTidy <- grpGSEAres %>% as_tibble() %>% arrange(desc(NES))

grpGSEAresTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
ggplot(grpGSEAresTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


  # Write cls first line 
writeLines(c("6 2 1"), "shp.cls", sep = "\t")




  # export rnk files


shpUPrnk <- read.csv("/home/kaylab/Documents/projects/BRCA-PI/csv/shpSigUPByLFC.csv")
shpUPrnk <- shpUPrnk[,c("Symbol","log2FoldChange")]
write.table(shpUPrnk, file = "shpSigUPByLFC.rnk", quote = F, sep = "\t", row.names = FALSE)

shpDNrnk <- read.csv("/home/kaylab/Documents/projects/BRCA-PI/csv/shpSigDNByLFC.csv")
shpDNrnk <- shpDNrnk[,c("Symbol","log2FoldChange")]
write.table(shpDNrnk, file = "shpSigDNByLFC.rnk", quote = F, sep = "\t", row.names = FALSE)


#################  GO enrichment using "gprofiler"  #################

library(gprofiler2)

  # Get gene names of DE genes only in combo but not in shp or carf
comboOnly <- read.csv("./comboOnly_names.csv")
comboOnly1 <- as.vector(comboOnly[["MSRB1"]])
comboOnly1 <- trimws(comboOnly1)


gostres <- gost(query = comboOnly1, organism = "hsapiens", ordered_query = FALSE)
head(gostres$result,3)
GOunorder <- gostplot(gostres, capped = FALSE, interactive = TRUE)
GOunorder

shpGO <- gostres
shpGO <- gostplot(shpGO, capped = FALSE, interactive = TRUE)
carfGO <- gostres
comboGO <- GO
comboGO
shpGOPlt <- publish_gostplot(shpGO, highlight_terms = c("GO:0005201","GO:0070888", "GO:0062023", "GO:0009653", "KEGG:04710","TF:M12351_1"), width = NA,)
carfGOPlt <- publish_gostplot(GO, highlight_terms = c("TF:M00716_1","TF:M10438_1", "TF:M00695", "TF:M09894_1"), width = NA,)
comboGOPlt <- publish_gostplot(GO)
comboOnlyPlt <- publish_gostplot(GO, highlight_terms = c("TF:M00196","TF:M10072", "TF:M10071", "TF:M00931"), width = NA,)






