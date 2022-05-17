library(edgeR)
setwd("C:\\Users\\Y406\\Documents")
#Import gene expression value or degree value data
exprSet_all <- read.csv("C:\\Users\\Y406\\Documents\\cndmAlpha.csv")
exprSet <- exprSet_all[,-1]
group_list <- factor(c(rep("ND.Alpha",138),rep("T2D.Alpha",101)))
exprSet <- DGEList(counts = exprSet, group = group_list)
exprSet <- calcNormFactors(exprSet)
exprSet <- estimateCommonDisp(exprSet)
exprSet <- estimateTagwiseDisp(exprSet)
et <- exactTest(exprSet)
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)
write.csv(tTag,file = "NDvsT2Ddu.Alpha.csv")
