rm(list=ls())
T<- read.table("D:\\T2D_beta_avr.txt",head=F)
T<-as.matrix(T)
data1<-read.csv("D:\\lie.csv",header=F)#Add the names of genes in the first column     
data2=data.frame(data1,T)
data3=read.csv("D:\\hang.csv",header=F)#Add the names of cells in the first row
data4=t(data.frame(data3,t(data2)))
colnames(data4)=data4[1,]
data4=data4[-1,]
rownames(data4)=data4[,1]
data4=data4[,-1]
data5=data.frame(data4)

data5[]=lapply(data5, as.numeric)
View(data5)
library(WGCNA)

cyt1= exportNetworkToCytoscape(data5,threshold = 0.3,
                               edgeFile="D:\\edge_T2D_beta_0.3.txt",
                               nodeFile="D:\\node_T2D_beta_0.3.txt",weighted = TRUE)#The threshold value is variable based on the network
