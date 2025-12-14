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
.dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#Get the dataset
D <- readRDS("./Data/CREST_aligned_to_human_rhesus_marmoset.rds")
D$FineAno <- Idents(D)
#This contains fine anotation labels, so lets combined some labels for a broader anotation
Idents(D,cells=WhichCells(D,idents=c("CTB","eCTB","STB","EVT","Tb_CS7","Tb_CS6","Tb_CS3","Tb_CS4","Tb_CS5"))) <- "Trophoblast"
Idents(D,cells=WhichCells(D,idents=c("Amnion","Am_CS5","Am_CS6","Am_CS7","PGC_CS5"))) <- "Amnion"
Idents(D,cells=WhichCells(D,idents=c("YolkSac","SYS_CS5","SYS_CS6","SYS_CS7"))) <- "Yolk Sac"
Idents(D,cells=WhichCells(D,idents=c("Hypoblast","Hyp_CS3","VE_CS5","VE_CS6","VE_CS7","VE_CS4"))) <- "Hypoblast"
Idents(D,cells=WhichCells(D,idents=c("ExEmbMes","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","Stalk_CS7"))) <- "Mesehncyme"
Idents(D,cells=WhichCells(D,idents=c("EmDisc_CS5","EmbDisc","EmDisc_CS6","EmDisc_CS7","EmDiscPS_CS6","Epi_CS4","Epi_CS3"))) <- "Embryonic disc"
Idents(D,cells=WhichCells(D,idents=c("PGC_CS6","PGC_CS7"))) <- "Embryonic disc"

#Colours for plotting based on broad anos
cols<-c("Embryonic disc"="#437abe",
        "Amnion"="#5ec7cd",
        "Hypoblast"="#f58b67",
        "Yolk Sac"="#2b8a42",
        "Mesehncyme"="#9fcb40",
        "Trophoblast"="#e1c122")
  

#This is the UMAP aligned to the composite reference (Wang et al, and Mole et al) used in Figure 5A
p<-DimPlot(D, cols=cols,  pt.size = 4, reduction = "pca", split.by = "Species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/FigureS5A.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

#Now let's look at exmese, sys and VE identities
Idents(D) <- D$FineAno
Dsub <- subset(D,idents=c("YolkSac","ExEmbMes","Hypoblast","ExMes_CS5","VE_CS5","SYS_CS5","Hyp_CS3","ExMes_CS6","SYS_CS6","VE_CS6","SYS_CS7","ExMes_CS7","VE_CS7","VE_CS4"))
Dsub <- RunPCA(Dsub, npcs = 20, verbose = FALSE)
Dsub <- RunUMAP(Dsub, reduction = "pca", dims = 1:20, n.neighbors = 50)
Dsub <- FindNeighbors(Dsub, reduction = "pca", dims = 1:20)#

#Colours for all of the fine anotations
my_cols <- c('YolkSac'='#49873b',    
             'ExEmbMes'='#a9d058',  
             'Hypoblast'='#f0926d',     
             'Hyp_CS3'='#dfbfb5',       
             'VE_CS5'='#dfbfb5',         
             'VE_CS6'='#dfbfb5',         
             'VE_CS7'='#dfbfb5',       
             'SYS_CS5'='#6ab662',        
             'SYS_CS6'='#6ab662',        
             'SYS_CS7'='#6ab662',        
             'ExMes_CS5'='#bae366',      
             'ExMes_CS6'='#bae366',    
             'ExMes_CS7'='#bae366',      
             'VE_CS5'='#f5d7ce',       
             'SYS_CS6'='#70c967',
             'SYS_CS7'='#70c967',        
             'ExMes_CS5'='#c5ed74',      
             'ExMes_CS6'='#c5ed74',     
             'ExMes_CS7'='#c5ed74',      
             'Hyp_CS3'='#fce3dc',
             'VE_CS4'='#fce3dc',  
             'VE_CS5'='#fce3dc',  
             'SYS_CS6'='#fce3dc') 

#Factor for plotting different shapes
Dsub$cellshape <- as.factor(as.double(as.factor(Dsub$Species)))
cellshape <- as.character(Dsub$cellshape)
cellshape[which(cellshape=="3")] <- "5"
cellshape[which(cellshape=="4")] <- "3"
cellshape[which(cellshape=="5")] <- "4"
Dsub$cellshape <- as.factor(cellshape)

#Supplementary Figure 5B
p<-DimPlot(Dsub, cols=my_cols, shape.by = "cellshape", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SF5B.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(Dsub, cols=my_cols, shape.by = "cellshape", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/SF5B_PC2_3.pdf",sep=""),width = 10, height = 8,p)

#Look at some markers (using RNA slot)
DefaultAssay(Dsub) <- "RNA"
p<-FeaturePlot(Dsub,split.by = "Species", reduction = "pca", pt.size = 4, features = c("LEFTY2","ANXA1","CER1","COTL1"))
ggsave(filename=paste(saveext,"/DimRed/SF5D.pdf",sep=""),width = 40, height = 30,p)

Dsub$stash <- Idents(Dsub)
Idents(Dsub) <- paste(Dsub$Species,Dsub$stash,sep="")

#Now we are going to find some DE genes SYS vs VE
Dsub <- JoinLayers(Dsub)
DefaultAssay(Dsub) <- "RNA" 
Mrk1 <- FindMarkers(Dsub, ident.1 = "HumanYolkSac", ident.2 = "HumanHypoblast")
Mrk2 <- FindMarkers(Dsub, ident.1 = c("MarmSYS_CS5","MarmSYS_CS6","MarmSYS_CS7"), ident.2 = c("MarmVE_CS5","MarmVE_CS6","MarmVE_CS7") )
Mrk3 <- FindMarkers(Dsub, ident.1 = c("CynoSYS_CS6","CynoSYS_CS7"), ident.2 = "CynoVE_CS5")

Mrk1A <- FindMarkers(Dsub, ident.1 = "HumanYolkSac", ident.2 = "HumanHypoblast",only.pos=TRUE)
Mrk2A <- FindMarkers(Dsub, ident.1 = c("MarmSYS_CS5","MarmSYS_CS6","MarmSYS_CS7"), ident.2 = c("MarmVE_CS5","MarmVE_CS6","MarmVE_CS7") ,only.pos=TRUE)
Mrk3A <- FindMarkers(Dsub, ident.1 = c("CynoSYS_CS6","CynoSYS_CS7"), ident.2 = "CynoVE_CS5",only.pos=TRUE)

Mrk1B <- FindMarkers(Dsub, ident.2 = "HumanYolkSac", ident.1 = "HumanHypoblast",only.pos=TRUE)
Mrk2B <- FindMarkers(Dsub, ident.2 = c("MarmSYS_CS5","MarmSYS_CS6","MarmSYS_CS7"), ident.1 = c("MarmVE_CS5","MarmVE_CS6","MarmVE_CS7"),only.pos=TRUE )
Mrk3B <- FindMarkers(Dsub, ident.2 = c("CynoSYS_CS6","CynoSYS_CS7"), ident.1 = "CynoVE_CS5",only.pos=TRUE)

newdata1 <- Mrk2A[order(Mrk2A$avg_log2FC,decreasing = TRUE),]
newdata2 <- Mrk2B[order(Mrk2B$avg_log2FC,decreasing = TRUE),]
newdata3 <- Mrk3A[order(Mrk3A$avg_log2FC,decreasing = TRUE),]
newdata4 <- Mrk3B[order(Mrk3B$avg_log2FC,decreasing = TRUE),]

DefaultAssay(Dsub) <- "RNA"

#Add module scores
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata1)[1:100] ), name = "SYSMarm")
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata2)[1:100] ), name = "VEMarm")
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata3)[1:100] ), name = "SYSCyno")
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata4)[1:100] ), name = "VECyno")

#Module score deltas
Dsub$SYSDelta1 <- Dsub$SYSMarm1 - Dsub$VEMarm1
Dsub$SYSDelta2 <- Dsub$SYSCyno1 - Dsub$VECyno1

#Now do some MA plots
Idents(Dsub) <- Dsub$Species
onlyOurs <- subset(Dsub,idents="Human")
Idents(onlyOurs) <- onlyOurs$FineAno

Av <- AverageExpression(onlyOurs)
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
Ae5$AvExp <- log(0.5*(Ae5$'YolkSac' + Ae5$'Hypoblast') +1)

genes.to.label <- c("VCAN","ANKRD1","ID3","TNFRSF19","S100A11","SAT1","GSTP","FTL","ANXA2","TMSB10","NEAT1","FTH1","APOE","APOA2","DIO3","IHH","NODAL","CCKBR","RBP4")

p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = c(genes.to.label), repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"/DimRed/SF5E.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)


#Find DE genes ExMes vs others
Idents(Dsub) <- paste(Dsub$Species,Dsub$stash,sep="")
Mrk2 <- FindMarkers(Dsub, ident.1 = c("MarmExMes_CS5","MarmExMes_CS6","MarmExMes_CS7"), ident.2 = c("MarmVE_CS5","MarmVE_CS6","MarmVE_CS7","MarmSYS_CS5","MarmSYS_CS6","MarmSYS_CS7") )
Mrk3 <- FindMarkers(Dsub, ident.1 = c("CynoExMes_CS5","CynoExMes_CS6","CynoExMes_CS7"), ident.2 = c("CynoVE_CS5","CynoSYS_CS6","CynoSYS_CS7") )
Mrk1 <- FindMarkers(Dsub, ident.1 = c("CynoExMes_CS5","CynoExMes_CS6","CynoExMes_CS7"), ident.2 = c("CynoVE_CS5","CynoSYS_CS6","CynoSYS_CS7") )

Mrk2A <- FindMarkers(Dsub, ident.1 = c("MarmExMes_CS5","MarmExMes_CS6","MarmExMes_CS7"), ident.2 = c("MarmVE_CS5","MarmVE_CS6","MarmVE_CS7","MarmSYS_CS5","MarmSYS_CS6","MarmSYS_CS7"), only.pos=TRUE )
Mrk3A <- FindMarkers(Dsub, ident.1 = c("CynoExMes_CS5","CynoExMes_CS6","CynoExMes_CS7"), ident.2 = c("CynoVE_CS5","CynoSYS_CS6","CynoSYS_CS7"), only.pos=TRUE )
Mrk1A <- FindMarkers(Dsub, ident.1 = c("CynoExMes_CS5","CynoExMes_CS6","CynoExMes_CS7"), ident.2 = c("CynoVE_CS5","CynoSYS_CS6","CynoSYS_CS7"), only.pos=TRUE )

Mrk2B <- FindMarkers(Dsub, ident.2 = c("MarmExMes_CS5","MarmExMes_CS6","MarmExMes_CS7"), ident.1 = c("MarmVE_CS5","MarmVE_CS6","MarmVE_CS7","MarmSYS_CS5","MarmSYS_CS6","MarmSYS_CS7"), only.pos=TRUE )
Mrk3B <- FindMarkers(Dsub, ident.2 = c("CynoExMes_CS5","CynoExMes_CS6","CynoExMes_CS7"), ident.1 = c("CynoVE_CS5","CynoSYS_CS6","CynoSYS_CS7"), only.pos=TRUE )
Mrk1B <- FindMarkers(Dsub, ident.2 = c("CynoExMes_CS5","CynoExMes_CS6","CynoExMes_CS7"), ident.1 = c("CynoVE_CS5","CynoSYS_CS6","CynoSYS_CS7"), only.pos=TRUE )

newdata1 <- Mrk2A[order(Mrk2A$avg_log2FC,decreasing = TRUE),]
newdata2 <- Mrk2B[order(Mrk2B$avg_log2FC,decreasing = TRUE),]
newdata3 <- Mrk3A[order(Mrk3A$avg_log2FC,decreasing = TRUE),]
newdata4 <- Mrk3B[order(Mrk3B$avg_log2FC,decreasing = TRUE),]

#Add module scores
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata1)[1:100] ), name = "ExMesMarm")
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata2)[1:100] ), name = "SYSExMesMarm")
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata3)[1:100] ), name = "ExMesCyno")
Dsub <- AddModuleScore(Dsub, features = list( rownames(newdata4)[1:100] ), name = "SYSExMesCyno")

#Module score deltas
Dsub$SYSDelta3 <- Dsub$ExMesMarm1 - Dsub$SYSExMesMarm1
Dsub$SYSDelta4 <- Dsub$ExMesCyno1 - Dsub$SYSExMesCyno1

Dsub$stash <- Idents(Dsub)

Idents(Dsub) <- Dsub$Species
onlyOurs <- subset(Dsub,idents="Human")

#Plots for SF5C
#p1<-FeaturePlot(onlyOurs, features = "SYSDelta1", pt.size = 4, reduction = "pca",  label = TRUE, cols =  c("blue", "red"), min.cutoff = -0.025, max.cutoff = 0.025, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/SYSvVEDeltaMarm.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(onlyOurs, features = "SYSDelta2", pt.size = 4, reduction = "pca", label = TRUE, cols =  c("blue", "red"), min.cutoff = -0.05, max.cutoff = 0.05, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SF5C_pannel1.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
#p1<-FeaturePlot(onlyOurs, features = "SYSDelta3", pt.size = 4, reduction = "pca",  label = TRUE, cols =  c("blue", "red"), min.cutoff = -0.025, max.cutoff = 0.025, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/ExMesDeltaMarm.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(onlyOurs, features = "SYSDelta4", pt.size = 4, reduction = "pca", label = TRUE, cols =  c("blue", "red"), min.cutoff = -0.05, max.cutoff = 0.05, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SF5C_pannel2.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

#Checked up until here. This needs editing.
#Now do some MA plots
Idents(onlyOurs) <- onlyOurs$FineAno

Mrk1 <- FindMarkers(onlyOurs, ident.1 = "YolkSac", ident.2 = "ExEmbMes")

Av <- AverageExpression(onlyOurs)
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
Ae5$AvExp <- log(0.5*(Ae5$'YolkSac' + Ae5$'ExEmbMes') +1)

genes.to.label1 <- c("ANK2","CER1","LINC00458","TRABD2B","GPC6","MYL4","GATA4","LRG5","SLPI","FTL","PTMA","S100A4","TMSB10","S100A6","ANAXA2","KRT18","FN1","TAGLN","NPPB","TRPC5","CDH6","PAMR1","MT1H")
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = c(genes.to.label1), repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"/DimRed/SF5F.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

