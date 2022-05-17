label<-read.csv("E:\\pdata_endocrine_cells.csv")#the id of endocrine cells
rownames(label)=label[,1]

data<-read.csv("E:\\T2Dexp_endocrine_cells.csv")#Gene expression matrix of endocrine cells
colnames(data)[1]<-'Gene'
rownames(data)=data[,1]

write.table(data, "E:\\test_counts.txt",  row.names=F, sep='\t',quote = F)
write.table(label, "E:\\test_meta.txt", row.names=F, sep='\t',quote = F) 


##Export the expression matrix and cell phenotype information to CellPhoneDB for cell communication
#conda activate cellphonedb   #activate cellphonedb
#test_counts[1:4,1:4]
#head test_meta.txt
#cut -f 2 test_meta.txt |sort |uniq -c 
#cellphonedb --help
#cellphonedb method statistical_analysis test_meta.txt test_counts.txt --counts-data=gene_name  # 如果是直接写出基因名的加最后这个参数，转化为基因ID的话不用加

#deconvoluted.txt：Average expression of genes 
#mean.txt：Average expression per receptor-ligand pair
#pvalues.txt：p-value for each receptor-ligand pair
#significant_means.txt：Average expression value for each receptor-ligand significance result

setwd('E:\\cell-cell_communication')
mypvals <- read.table("out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
mymeans <- read.table("out/means.txt",header = T,sep = "\t",stringsAsFactors = F) 

kp = grepl(pattern = "Beta", colnames(mypvals))

table(kp)
pos = (1:ncol(mypvals))[kp] 
choose_pvalues <- mypvals[,c(c(1,2,5,6,8,9),pos)]
choose_means <- mymeans[,c(c(1,2,5,6,8,9),pos)]

logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<0.05, 1, sum) 
# Only some interacting pairs with cell specificity are retained

logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]

# The same condition keeps choose_means
choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]

head(choose_means,1)
head(choose_pvalues,1)
dim(mypvals)
dim(mymeans)
dim(choose_means)
dim(choose_pvalues)

library(tidyverse)
meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = meansdf$interacting_pair,
                      CC = meansdf$variable,
                      means = meansdf$value)

pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = pvalsdf$interacting_pair,
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)

# Merge p-value and mean files
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >3.5))$means)
head(pldf)
pcc =  pldf%>% filter(means >3.5) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 1  )+ 
  theme_bw()+ 
  # scale_color_manual(values = rainbow(100))+
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8))
pcc

