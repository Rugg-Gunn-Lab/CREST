ibrary(Seurat)
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

Data <- readRDS("Data/CREST_anotated.rds")
Idents(Data,WhichCells(Data,idents=c("Epithelia_gland","Epithelia_SOX","Epithelia_cil","Epithelia_lum"))) <- "Epithelial"
Idents(Data,WhichCells(Data,idents=c("Stroma_decid","Stroma_prolif"))) <- "Stroma"
Data <- subset(Data,idents = c("Epithelial","Stroma"))
DefaultAssay(Data) <- "RNA"

ListC3 <- FindMarkers(Data,ident.2 = c("Stroma"),ident.1 = c("Epithelial"), test.use = "MAST")


AvExp <- AverageExpression(Data)
#AvExp <- AverageExpression(mammal.combined3)
AvExp <- as.data.frame(AvExp$RNA)
AvExp$LFC <- NA
AvExp$pval <- NA
AvExp[rownames(ListC3),"LFC"] <- ListC3$avg_log2FC
AvExp[rownames(ListC3),"pval"] <- ListC3$p_val_adj
AvExp$AvExp <- log2((AvExp$Epithelial+AvExp$Stroma)/2)

genes.to.label <- c("CFAP126","C5orf49","ZBBX","DTHD1","CFAP47","CROCC2","KRT5","PIFO","MUC4","KRT7","PAEP","KRT18","CYBA","ANXA2","GSTP1","TMSB10","FTH1",
                    "SERF2","MARCKS","COL4A2","COL6A1","PTMS","CALD1","RAB31","AXL","VIM","MMP2","CNIH3","IGFBP5","TPM2","DCN","AFF3","COL6A3","LRFN5","ACTG2","THY1","COL8A1")

p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Stroma_v_Epithelia_secMrks_test.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

