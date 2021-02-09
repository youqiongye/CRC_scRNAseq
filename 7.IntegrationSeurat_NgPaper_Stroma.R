library(Seurat)
library(ggplot2)
library(magrittr)
library(Matrix)
library(wordspace)
library(Hmisc)
library(rsvd)
options(stringsAsFactors = F)
set.seed(42)
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/")
Samples <- c("p7_N","p9_N", "p10_N", "p11_N", "p12_N", "p7_T",  "p9_T",  "p10_T", "p11_T", "p12_T")

ClustersColor <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
MSCcolor <- ClustersColor[ClustersColor$SubTypes %in% "MSC",]

.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
SMC10X_Seurat <- readr::read_rds("/home/yye/Data/Public/GEO_data/scRNAseq/NG_CRC/GSE132465_SMC_scRNAseqCountSeurat.rds.gz")
KUL10X_Seurat <- readr::read_rds("/home/yye/Data/Public/GEO_data/scRNAseq/NG_CRC/GSE144735_KUL_scRNAseqCountSeurat.rds.gz")
###SMC
SMC10X_MSC <- subset(SMC10X_Seurat,cells=rownames(SMC10X_Seurat@meta.data[SMC10X_Seurat$Cell_subtype %in% 
                                                                            c("Stromal 1",'Stromal 2',"Stromal 3","Stromal 4",
                                                                              "Myofibroblasts","Pericytes","Smooth muscle cells") ,]))

SMC10X_MSC <- FindVariableFeatures(SMC10X_MSC)
SMC10X_MSC <- ScaleData(SMC10X_MSC)
SMC10X_MSC <- RunPCA(SMC10X_MSC)
SMC10X_MSC <- RunUMAP(SMC10X_MSC,dims = 1:30)
pdf("7.1Integration_NgPaper/SMC_Stromal_umap.pdf",width = 6,height = 4)
p <- DimPlot(SMC10X_MSC,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:6])
dev.off()
p2 <- FeaturePlot(SMC10X_MSC,reduction = "umap",features = "FAP")
cowplot::plot_grid(plotlist = list(p,p2))
FeaturePlot(SMC10X_MSC,reduction = "umap",features =c("ACTA2","ACTG2"))
FeaturePlot(WholeTissueList$MSC,reduction = "umap",features =c("ACTA2","ACTG2"))



##KUL

KUL10X_MSC <- subset(KUL10X_Seurat,cells=rownames(KUL10X_Seurat@meta.data[KUL10X_Seurat$Cell_subtype %in% 
                                                                            c("Stromal 1",'Stromal 2',"Stromal 3","Stromal 4","Myofibroblasts","Pericytes","Smooth muscle cells") ,]))

KUL10X_MSC <- FindVariableFeatures(KUL10X_MSC)
KUL10X_MSC <- ScaleData(KUL10X_MSC)
KUL10X_MSC <- RunPCA(KUL10X_MSC)
KUL10X_MSC <- RunUMAP(KUL10X_MSC,dims = 1:30)
pdf("7.1Integration_NgPaper/KUL10X_Stromal_umap.pdf",width = 6,height = 4)
DimPlot(KUL10X_MSC,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:5])
dev.off()
p1 <-  DimPlot(SMC10X_MSC,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:5])+NoLegend()
p2 <-  DimPlot(KUL10X_MSC,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:5])+NoLegend()
tiff("7.1Integration_NgPaper/KUL10X_SMC10X_Stromal_umap.tiff",width = 1400,height = 650,res=200)
cowplot::plot_grid(plotlist = list(p1,p2),ncol=2)
dev.off()


MSCseurat <- WholeTissueList$MSC
MSCseurat@meta.data <- MSCseurat@meta.data[,c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes")]
MSCseurat@meta.data$Class <- rep("Our_10X",times=nrow(MSCseurat@meta.data))


SMC10X_MSC@meta.data <- SMC10X_MSC@meta.data[,c("Patient","nCount_RNA", "nFeature_RNA","Class","Cell_subtype")]%>% set_colnames(c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes"))
SMC10X_MSC@meta.data$Class <- rep("SMC10X",times=nrow(SMC10X_MSC@meta.data))

KUL10X_MSC@meta.data <- KUL10X_MSC@meta.data[,c("Patient","nCount_RNA", "nFeature_RNA","Class","Cell_subtype")]%>% set_colnames(c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes"))
KUL10X_MSC@meta.data$Class <- rep("KUL10X",times=nrow(KUL10X_MSC@meta.data))



###Integration 10 and SMARTseq
SeuratList <- list(MSCseurat,SMC10X_MSC,KUL10X_MSC)
names(SeuratList) <- c("Our_10X","SMC10X","KUL10X")


Seurat.features <- SelectIntegrationFeatures(object.list = SeuratList, nfeatures = 2000)

Seurat.anchors <- FindIntegrationAnchors(object.list = SeuratList,anchor.features=Seurat.features, dims = 1:30)

Seurat.integrated <- IntegrateData(anchorset = Seurat.anchors, dims = 1:30)
DefaultAssay(Seurat.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Seurat.integrated <- ScaleData(Seurat.integrated, verbose = FALSE)
#running PCA
Seurat.integrated <- RunPCA(Seurat.integrated, npcs = 50, verbose = T)

#jackstraw
#Seurat.integrated = JackStraw(Seurat.integrated, dims = max(1:50))
#Seurat.integrated = ScoreJackStraw(Seurat.integrated, dims = 1:50)
#dims.use = Seurat.integrated@reductions$pca@jackstraw@overall.p.values
#dims.use = dims.use[dims.use[, 2] < 0.0001, 1] #taking significant PCA and assigning dims for UMAP

#running tsne and Umap 
Seurat.integrated <- RunTSNE(Seurat.integrated, reduction = "pca", dims = 1:30,#dims.use,
                             do.fast = T, k.seed = 10, check_duplicates = FALSE,
                             perplexity = 30)
Seurat.integrated <- RunUMAP(Seurat.integrated, dims = 1:30)#dims.use)
Seurat.integrated$DefineTypes2 <- ifelse(Seurat.integrated$Class %in% "Our_10X",
                                         Seurat.integrated$DefineTypes,
                                         paste(gsub("10X","",Seurat.integrated$Class),Seurat.integrated$DefineTypes,sep="_")
)
Seurat.integrated$DefineTypes2 <- factor(Seurat.integrated$DefineTypes2,
                                         levels=c(MSCcolor$DefineTypes,grep("SMC|KUL",names(table(Seurat.integrated$DefineTypes2)),value = T)))
Seurat.integrated$Class <- factor(Seurat.integrated$Class,levels=c("Our_10X","SMC10X","KUL10X"))
readr::write_rds(Seurat.integrated,path="7.1Integration_NgPaper/StromalSeurat.integrated.rds.gz",compress="gz")
Seurat.integrated <- readr::read_rds("7.1Integration_NgPaper/StromalSeurat.integrated.rds.gz")

Seurat.integrated@meta.data$DefineTypes2 <- factor(Seurat.integrated@meta.data$DefineTypes2,levels=c(MSCcolor$DefineTypes,grep("KUL|SMC",names(table(Seurat.integrated$DefineTypes2)),value=T)))
pdf("7.1Integration_NgPaper/Seurat.integrated_SplitByDataset.pdf",width = 15,height = 5)
DimPlot(Seurat.integrated,group.by="DefineTypes2",split.by = "Class",cols = c(MSCcolor$Colors,ClustersColor$Colors[1:6] ,ClustersColor$Colors[1:6]),label=F,label.size =2)
dev.off()
pdf("7.1Integration_NgPaper/Seurat.integrated_SplitByDatasetLegend.pdf",width = 15,height = 5)
DimPlot(Seurat.integrated,group.by="DefineTypes2",split.by = "Class",na.value = NA,cols = rep(NA,22),label=F,pt.size = 0)
dev.off()
tiff("7.1Integration_NgPaper/Seurat.integrated_SplitByDataset.tiff",width = 1500,height = 600,res=220)
DimPlot(Seurat.integrated,group.by="DefineTypes2",split.by = "Class",cols = c(MSCcolor$Colors,ClustersColor$Colors[1:6] ,ClustersColor$Colors[1:6]),label=F,label.size =2)+NoLegend()

dev.off()

Idents(Seurat.integrated) <- "DefineTypes2"
DefaultAssay(Seurat.integrated) <- "RNA"
tiff("7.1Integration_NgPaper/MyofibroblastFeatures.tiff",width = 2100,height = 4200,res=250)
FeaturePlot(Seurat.integrated,feature=c("MYH11","ACTG2","TAGLN","MCAM","RGS5","PDGFRB"),split.by = "Class",label=F,label.size =2)
dev.off()
FeaturePlot(Seurat.integrated,feature=c("FAP"),split.by = "Class",label=T,label.size =2)

Idents(SMC10X_MSC)  <- "Cell_subtype"
FeaturePlot(SMC10X_MSC,feature=c("MYH11","ACTG2","TAGLN"),label=T,label.size =2)

WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtype_AllFunctionTF.rds.gz")
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")


###Correlation
ClusterCols <- data.frame(unique(Seurat.integrated@meta.data[,c("Class","DefineTypes2")]),
                          Cols= c(MSCcolor$Colors,ClustersColor$Colors[1:6] ,ClustersColor$Colors[1:6]),
                          ClassCols = rep(c("red","blue","darkgreen"),times=c(10,6,6))) %>% 
  dplyr::mutate(Class=rep(c("OurDataset","SMC","KUL"),times=c(10,6,6)),DefineTypes2=gsub("SMC_|KUL_","",DefineTypes2))

ClusterID <- lapply(split(Seurat.integrated@meta.data,list(Seurat.integrated@meta.data$DefineTypes2)), function(x)rownames(x))
Data <- Seurat.integrated@assays$RNA@data %>% data.frame
#colnames(Data) <- gsub("\\.","-",colnames(Data))
ClusterMean <- t(apply(Data,1,function(x){
  lapply(ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))

ClusterMean <- ClusterMean[apply(ClusterMean,1,sum)>0,]
colnames(ClusterMean) <- ifelse(colnames(ClusterMean) %in% MSCcolor$DefineTypes,paste("OurDataset",colnames(ClusterMean),sep="_"),colnames(ClusterMean))
ClusterCols <- ClusterCols[match(colnames(ClusterMean),paste(ClusterCols$Class,ClusterCols$DefineTypes2,sep="_")),]
SampleColA <- ClusterCols %>%
  dplyr::select(c(Cols,ClassCols)) %>% as.matrix() %>% set_colnames(c("DefineTypes2","Dataset")) %>%
  set_rownames(names(ClusterID))


TypeCor <- cor(ClusterMean,method="spearman")


pdf("7.1Integration_NgPaper/DiffDataSet_StromalCor.pdf",width = 6,height = 7)
p <- heatmap.plus::heatmap.plus(TypeCor, ColSideColors = SampleColA,
                                col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100),
                                labRow = NA,cexCol = 1,margins = c(15,10))
#legend("topright",legend= c(unique(ClusterCols$Class),"", ClusterCols[ClusterCols$Class %in% c("OurDataset","SMC"),]$DefineTypes2),
#       fill=c(c("red","blue","darkgreen"),"white", ClusterCols[ClusterCols$Class %in% c("OurDataset","SMC"),]$Cols), border=FALSE, bty="n", y.intersp = 0.8, cex=0.9,xpd=T)
dev.off()


###
###Correlation
ClusterCols <- data.frame(unique(Seurat.integrated@meta.data[,c("Class","DefineTypes2")]),
                          Cols= c(MSCcolor$Colors,ClustersColor$Colors[1:6] ,ClustersColor$Colors[1:6]),
                          ClassCols = rep(c("red","blue","darkgreen"),times=c(10,6,6))) %>% 
  dplyr::mutate(Class=rep(c("OurDataset","SMC","KUL"),times=c(10,6,6)),DefineTypes2=gsub("SMC_|KUL_","",DefineTypes2))

ClusterID <- lapply(split(Seurat.integrated@meta.data,list(Seurat.integrated@meta.data$DefineTypes2)), function(x)rownames(x))
Data <- Seurat.integrated@assays$integrated@data %>% data.frame
#colnames(Data) <- gsub("\\.","-",colnames(Data))
IntClusterMean <- t(apply(Data,1,function(x){
  lapply(ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))


IntClusterMean <- IntClusterMean[apply(IntClusterMean,1,sum)>0,]
colnames(IntClusterMean) <- ifelse(colnames(IntClusterMean) %in% MSCcolor$DefineTypes,paste("OurDataset",colnames(IntClusterMean),sep="_"),colnames(IntClusterMean))
ClusterCols <- ClusterCols[match(colnames(IntClusterMean),paste(ClusterCols$Class,ClusterCols$DefineTypes2,sep="_")),]
SampleColA <- ClusterCols %>%
  dplyr::select(c(Cols,ClassCols)) %>% as.matrix() %>% set_colnames(c("DefineTypes2","Dataset")) %>%
  set_rownames(names(ClusterID))


TypeCor <- cor(IntClusterMean,method="spearman")


pdf("7.1Integration_NgPaper/IntDiffDataSet_StromalCor.pdf",width = 6,height = 7)
p <- heatmap.plus::heatmap.plus(TypeCor, ColSideColors = SampleColA,
                                col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100),
                                labRow = NA,cexCol = 1,margins = c(15,10))
#legend("topright",legend= c(unique(ClusterCols$Class),"", ClusterCols[ClusterCols$Class %in% c("OurDataset","SMC"),]$DefineTypes2),
#       fill=c(c("red","blue","darkgreen"),"white", ClusterCols[ClusterCols$Class %in% c("OurDataset","SMC"),]$Cols), border=FALSE, bty="n", y.intersp = 0.8, cex=0.9,xpd=T)
dev.off()
