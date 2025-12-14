library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
library(stringr)

set.seed(1) 

saveext = "./CREST/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#Get the dataset
D <- readRDS("Data/CREST_trophoblast_aligned.rds")

#Colours for plotting
cols<-c('CTB'='#e2c120',
'mCTb'='#e2c120',
'cCTb'='#e2c120',
'CTb'='#e2c120',
'STB'='#d68c34',
'mSTb'='#d68c34',
'cSTB'='#d68c34',
'eSTB'='#d68c34',
'STb'='#d68c34',
'EVT'='#bf3737',
'cEVT'='#bf3737',
'cEVT2'='#bf3737',
'Evt'='#bf3737',
'Tb_CS3'='#db91b3',
'Tb_CS4'='#db91b3')


#Plot the UMAP.
#This is the UMAP aligned to the composite reference (Wang et al, and Mole et al) used in Figure 1
p<-DimPlot(D, cols=cols,  pt.size = 4, reduction = "pca", split.by = "Species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Figure6D.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

D <- JoinLayers(D)

Av1 <- AverageExpression(D)

OrderList <- readRDS("./Data/TrophoblastMarkers.rds")

d1 <- Av1$RNA[OrderList,c("CTB","STB","EVT")]
d2 <- Av1$RNA[OrderList,c("cCTb","cSTB","cEVT")]
d3 <- Av1$RNA[OrderList,c("CTb","STb","Evt")]
#d4 <- Av3$RNA[OrderList,c("pVCT","pEVT","pSCT")]
d1 <- t(scale(t(d1) ))
d2 <- t(scale(t(d2) ))
d3 <- t(scale(t(d3) ))

d <- cbind(d1,d2,d3)
d <- d[OrderList,c("CTB","STB","EVT","cCTb","cSTB","cEVT","CTb","STb","Evt")]
pm6<-pheatmap(d1, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Figure6E_panel1",".pdf",sep="") ,width=2,height=10)
pm6<-pheatmap(d2, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Figure6E_panel2",".pdf",sep="") ,width=2,height=10)
pm6<-pheatmap(d3, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Figure6E_panel3",".pdf",sep="") ,width=2,height=10)
pm6<-pheatmap(d, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Figure6E",".pdf",sep="") ,width=4,height=10)

