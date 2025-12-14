library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
library(stringr)

set.seed(1) 

saveext = "~/Desktop/Data/Endometrial/InVitro/Matteo/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#Get the dataset
D <- readRDS("/Users/christopherpenfold/Downloads/CREST_epithelia.rds")

D$ID6 <- factor(D$ID6,levels=c("Ours","proliferative","early-secretory","late-secretory"))

#D$FineAno <- Idents(D)
#This contains fine anotation labels, so lets combined some labels for a broader anotation
#Idents(D,cells=WhichCells(D,idents=c("CTB","eCTB","STB","EVT","Tb_CS7","Tb_CS6","Tb_CS3","Tb_CS4","Tb_CS5"))) <- "Trophoblast"
#Idents(D,cells=WhichCells(D,idents=c("Amnion","Am_CS5","Am_CS6","Am_CS7","PGC_CS5"))) <- "Amnion"
#Idents(D,cells=WhichCells(D,idents=c("YolkSac","SYS_CS5","SYS_CS6","SYS_CS7"))) <- "Yolk Sac"
#Idents(D,cells=WhichCells(D,idents=c("Hypoblast","Hyp_CS3","VE_CS5","VE_CS6","VE_CS7","VE_CS4"))) <- "Hypoblast"
#Idents(D,cells=WhichCells(D,idents=c("ExEmbMes","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","Stalk_CS7"))) <- "Mesehncyme"
#Idents(D,cells=WhichCells(D,idents=c("EmDisc_CS5","EmbDisc","EmDisc_CS6","EmDisc_CS7","EmDiscPS_CS6","Epi_CS4","Epi_CS3"))) <- "Embryonic disc"
#Idents(D,cells=WhichCells(D,idents=c("PGC_CS6","PGC_CS7"))) <- "Embryonic disc"

#Colours for plotting based on broad anos
cols<-c("Embryonic disc"="#437abe",
        "Amnion"="#5ec7cd",
        "Hypoblast"="#f58b67",
        "Yolk Sac"="#2b8a42",
        "Mesehncyme"="#9fcb40",
        "Trophoblast"="#e1c122")
  
mycols <- c("Ciliated"="#ad3479",
"Glandular"="#ccbb44",
"SOX9"="#238a42",
"Luminal"="#65ccf1")

Idents(D) <- D$EpithelialAno
#This is the UMAP aligned to the composite reference (Wang et al, and Mole et al) used in Figure 5A
p<-DimPlot(D, cols=mycols,  pt.size = 4, reduction = "pca", split.by = "ID6", label = TRUE, repel = TRUE, dim = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/FigureS5H.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

Idents(D) <- D$ID6
Dsub <- subset(D,idents=c("Ours"))
p1<-FeaturePlot(Dsub, features = "Gland1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/SF5_H_pannel2.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(Dsub, features = "Cil1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/SF5_H_pannel3.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(Dsub, features = "Luminal1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/SF5_H_pannel4.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

#Checked up until here. Needs editing for later plots
#Find DE genes ExMes vs others
Idents(Dsub) <- Dsub$EpithelialAno
Mrk1 <- FindMarkers(Dsub, ident.2 = c("Luminal"), ident.1 = c("Glandular"), only.pos=FALSE )

newdata1 <- Mrk1[order(Mrk1$avg_log2FC,decreasing = TRUE),]

Av <- AverageExpression(Dsub)
Ae5 <- as.data.frame(Av$RNA)
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5[rownames(Mrk1),"Pval1"] <- -log2(Mrk1$p_val_adj)
Ae5[rownames(Mrk1),"FC1"] <- Mrk1$avg_log2FC
pospos1 <- which( abs(Ae5$FC1)>log(1.1) & Ae5$Pval1>-log2(0.01))
Ae5$ColInd <- 0

Ae5[ which( (Ae5$FC1)>log(1.1) & Ae5$Pval1> -log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.1) & Ae5$Pval1> -log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
Ae5$AvExp <- log(0.5*(Ae5$'Luminal' + Ae5$'Glandular') +1)

#genes.to.label1 <- c("ANK2","CER1","LINC00458","TRABD2B","GPC6","MYL4","GATA4","LRG5","SLPI","FTL","PTMA","S100A4","TMSB10","S100A6","ANAXA2","KRT18","FN1","TAGLN","NPPB","TRPC5","CDH6","PAMR1","MT1H")
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = c(genes.to.label1), repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"/DimRed/SF5J.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

