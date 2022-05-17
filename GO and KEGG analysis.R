rm(list=ls())
setwd("C:\\Users\\Y406\\Documents\\gene enrichment analysis")
library(openxlsx)
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)
#Load differential expression data, only need gene name
info <- read.xlsx( "C:\\Users\\Y406\\Documents\\hub基因.xlsx", rowNames = F,colNames = T)
#Species library designated for enrichment analysis
GO_database <- 'org.Hs.eg.db' 
KEGG_database <- 'hsa' 
#Gene ID conversion
gene <- bitr(info$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
#GO enrichment analysis
erich.go.BP = enrichGO(gene = gene$ENTREZID,
                       OrgDb = GO_database,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.BP)
barplot(erich.go.BP)
erich.go.CC = enrichGO(gene = gene$ENTREZID,
                       OrgDb = GO_database,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.CC)
barplot(erich.go.CC)
erich.go.MF = enrichGO(gene = gene$ENTREZID,
                       OrgDb = GO_database,
                       keyType = "ENTREZID",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.MF)
barplot(erich.go.MF)
#KEGG functional enrichment analysis
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.3,
                 qvalueCutoff = 0.3)
dotplot(KEGG)
barplot(KEGG)
hsa_kegg<-download_KEGG("hsa")
str(hsa_kegg)
length(unique(hsa_kegg$KEGGPATHID2EXTID$from)) 
length(unique(hsa_kegg$KEGGPATHID2EXTID$to))  
x<-data.frame(hsa_kegg$KEGGPATHID2EXTID)  
length(unique(hsa_kegg$KEGGPATHID2NAME$from)) 
length(unique(hsa_kegg$KEGGPATHID2NAME$to)) 
y<-data.frame(hsa_kegg$KEGGPATHID2NAME)  
info<-gene
#ID conversion
deg.id <- bitr(info$SYMBOL,fromType = "SYMBOL",toType= "ENTREZID",OrgDb = "org.Hs.eg.db") 
#Merge gene name,ENTREZID
idvec <- deg.id$ENTREZID
names(idvec) <- deg.id$SYMBOL
info$ENTREZID <- idvec[info$SYMBOL]
write.csv(info, "hub gene_ENTREZID.csv", quote = F, row.names = F)
#KEGG functional enrichment analysis
ekk <- enrichKEGG(gene = info$ENTREZID,
                  keyType = 'kegg',
                  organism = 'hsa',
                  pvalueCutoff = 0.3,
                  pAdjustMethod  = "BH")
dim(ekk)
write.csv(ekk,"enrichKEGG_test1.csv",quote = F)
dotplot(ekk,showCategory=10)
barplot(ekk)
cnetplot(ekk, categorySize="pvalue", showCategory = 5)
library(ggupset)
upsetplot(ekk)
library(pathview)
#GO/KEGG enrichment histogram and dot plot
GO<-enrichGO( gene$ENTREZID,
              OrgDb = GO_database,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)
KEGG<-enrichKEGG(gene$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 0.3,
                 qvalueCutoff = 0.3)
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dotplot(KEGG)
#Association network diagram between enriched genes and their function set/pathway set
enrichplot::cnetplot(GO,circular=TRUE,colorEdge = TRUE)
enrichplot::cnetplot(KEGG,circular=TRUE,colorEdge = TRUE)
#Heatmap form
enrichplot::heatplot(GO,showCategory = 50)#基因-通路关联热图
enrichplot::heatplot(KEGG,showCategory = 50)
#Save KEGG-enriched pathways to a local file and select pathways for display
write.table(KEGG$ID, file = "KEGG_IDs.txt", 
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
KEGG<-read.csv("C:\\Users\\Y406\\Documents\\anjiyin\\enrichKEGG_test1.csv")
browseKEGG(KEGG,"hsa04370")



