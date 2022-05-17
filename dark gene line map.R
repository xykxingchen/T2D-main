rm(list=ls())
setwd("C:\\Users\\Y406\\Documents\\anjiyin")
library(ggplot2)
#Set as white background theme
theme_set(theme_bw())
#Import gene expression and degree value data separately
df1<-read.csv("C:\\Users\\Y406\\Documents\\TNFAIP6jiyin.csv",header = T)
df1$stage<-factor(df1$stage,levels = c("ND","T2D"))
df2<-read.csv("C:\\Users\\Y406\\Documents\\TNFAIP6du.csv",header = T)
df2$stage<-factor(df2$stage,levels = c("ND","T2D"))
df11<-read.csv("C:\\Users\\Y406\\Documents\\TNFAIP6jiyinaverage.csv",header = T)
df11$stage<-factor(df11$stage,levels = c("ND","T2D"))
df22<-read.csv("C:\\Users\\Y406\\Documents\\TNFAIP6duaverage.csv",header = T)
df22$stage<-factor(df22$stage,levels = c("ND","T2D"))
#Scatter plus line graph
g<-ggplot(df1, aes(x=stage, y=GEM, color = stage))
g+geom_point(position = "jitter", alpha=.3)+
  geom_line(data = df11,aes(x=stage, y=average),color="black",group=df11$group,size=0.5) +
  ylim(c(0,1))+
  labs(y="EXP", 
       x="stage", 
       title="TNFAIP6")
g<-ggplot(df2, aes(x=stage, y=NDM, color = stage))
g+geom_point(position = "jitter", alpha=.3)+
  geom_line(data = df22,aes(x=stage, y=average),color="black",group=df22$group,size=0.5) +
  ylim(c(0,1))+
  labs(y="Degree", 
       x="stage", 
       title="TNFAIP6")
