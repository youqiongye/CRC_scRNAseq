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


SMC10X_Myeloid <- subset(SMC10X_Seurat,cells=rownames(SMC10X_Seurat@meta.data[SMC10X_Seurat$Cell_type %in% "Myeloids" ,]))

SMC10X_Myeloid <- FindVariableFeatures(SMC10X_Myeloid)
SMC10X_Myeloid <- ScaleData(SMC10X_Myeloid)
SMC10X_Myeloid <- RunPCA(SMC10X_Myeloid)
SMC10X_Myeloid <- RunUMAP(SMC10X_Myeloid,dims = 1:30)
pdf("7.1Integration_NgPaper/SMC_Myeloid_umap.pdf",width = 6,height = 4)
p <- DimPlot(SMC10X_Myeloid,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:6])
p
dev.off()
p2 <- FeaturePlot(SMC10X_Myeloid,reduction = "umap",features = "FAP")
cowplot::plot_grid(plotlist = list(p,p2))
FeaturePlot(SMC10X_Myeloid,reduction = "umap",features =c("ACTA2","ACTG2"))
FeaturePlot(WholeTissueList$Myeloid,reduction = "umap",features =c("ACTA2","ACTG2"))



##KUL

KUL10X_Myeloid <- subset(KUL10X_Seurat,cells=rownames(KUL10X_Seurat@meta.data[KUL10X_Seurat$Cell_type %in% "Myeloids",]))

KUL10X_Myeloid <- FindVariableFeatures(KUL10X_Myeloid)
KUL10X_Myeloid <- ScaleData(KUL10X_Myeloid)

KUL10X_Myeloid <- RunPCA(KUL10X_Myeloid)
KUL10X_Myeloid <- RunUMAP(KUL10X_Myeloid,dims = 1:30)
pdf("7.1Integration_NgPaper/KUL10X_Myeloid_umap.pdf",width = 6,height = 4)
DimPlot(KUL10X_Myeloid,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:5])
dev.off()
p1 <-  DimPlot(SMC10X_Myeloid,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:5])+NoLegend()
p2 <-  DimPlot(KUL10X_Myeloid,reduction = "umap",group.by = "Cell_subtype",cols = ClustersColor$Colors[1:5])+NoLegend()
tiff("7.1Integration_NgPaper/KUL10X_SMC10X_Myeloid_umap.tiff",width = 1400,height = 650,res=200)
cowplot::plot_grid(plotlist = list(p1,p2),ncol=2)
dev.off()


Myeloidseurat <- WholeTissueList$Myeloid
Myeloidseurat@meta.data <- Myeloidseurat@meta.data[,c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes")]
Myeloidseurat@meta.data$Class <- rep("Our_10X",times=nrow(Myeloidseurat@meta.data))


SMC10X_Myeloid@meta.data <- SMC10X_Myeloid@meta.data[,c("Patient","nCount_RNA", "nFeature_RNA","Class","Cell_subtype")]%>% set_colnames(c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes"))
SMC10X_Myeloid@meta.data$Class <- rep("SMC10X",times=nrow(SMC10X_Myeloid@meta.data))

KUL10X_Myeloid@meta.data <- KUL10X_Myeloid@meta.data[,c("Patient","nCount_RNA", "nFeature_RNA","Class","Cell_subtype")]%>% set_colnames(c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes"))
KUL10X_Myeloid@meta.data$Class <- rep("KUL10X",times=nrow(KUL10X_Myeloid@meta.data))



###Integration 10 and SMARTseq
MyeloidSeuratList <- list(Myeloidseurat,SMC10X_Myeloid,KUL10X_Myeloid)
names(MyeloidSeuratList) <- c("Our_10X","SMC10X","KUL10X")


Seurat.features <- SelectIntegrationFeatures(object.list = MyeloidSeuratList, nfeatures = 2000)

Seurat.anchors <- FindIntegrationAnchors(object.list = MyeloidSeuratList,anchor.features=Seurat.features, dims = 1:30)

MyeloidSeurat.integrated <- IntegrateData(anchorset = Seurat.anchors, dims = 1:30)
DefaultAssay(MyeloidSeurat.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
MyeloidSeurat.integrated <- ScaleData(MyeloidSeurat.integrated, verbose = FALSE)
#running PCA
MyeloidSeurat.integrated <- RunPCA(MyeloidSeurat.integrated, npcs = 50, verbose = T)


#jackstraw
#MyeloidSeurat.integrated = JackStraw(MyeloidSeurat.integrated, dims = max(1:50))
#MyeloidSeurat.integrated = ScoreJackStraw(MyeloidSeurat.integrated, dims = 1:50)
#dims.use = MyeloidSeurat.integrated@reductions$pca@jackstraw@overall.p.values
#dims.use = dims.use[dims.use[, 2] < 0.0001, 1] #taking significant PCA and assigning dims for UMAP

#running tsne and Umap 
MyeloidSeurat.integrated <- RunTSNE(MyeloidSeurat.integrated, reduction = "pca", dims = 1:30,#dims.use,
                             do.fast = T, k.seed = 10, check_duplicates = FALSE,
                             perplexity = 30)
MyeloidSeurat.integrated <- RunUMAP(MyeloidSeurat.integrated, dims = 1:50)#dims.use)
MyeloidSeurat.integrated$DefineTypes2 <- ifelse(MyeloidSeurat.integrated$Class %in% "Our_10X",
                                         MyeloidSeurat.integrated$DefineTypes,
                                         paste(gsub("10X","",MyeloidSeurat.integrated$Class),MyeloidSeurat.integrated$DefineTypes,sep="_")
)
MyeloidSeurat.integrated$DefineTypes2 <- factor(MyeloidSeurat.integrated$DefineTypes2,
                                         levels=c(Myeloidcolor$DefineTypes,grep("SMC|KUL",names(table(MyeloidSeurat.integrated$DefineTypes2)),value = T)))
MyeloidSeurat.integrated$Class <- factor(MyeloidSeurat.integrated$Class,levels=c("Our_10X","SMC10X","KUL10X"))
readr::write_rds(MyeloidSeurat.integrated,path="7.1Integration_NgPaper/MyeloidMyeloidSeurat.integrated.rds.gz",compress="gz")
MyeloidSeurat.integrated <- readr::read_rds("7.1Integration_NgPaper/MyeloidMyeloidSeurat.integrated.rds.gz")
MyeloidSeurat.integrated@meta.data$DefineTypes2 <- factor(MyeloidSeurat.integrated@meta.data$DefineTypes2,levels=c(Myeloidcolor$DefineTypes,grep("KUL|SMC",names(table(MyeloidSeurat.integrated$DefineTypes2)),value=T)))
pdf("7.1Integration_NgPaper/MyeloidSeurat.integrated_SplitByDataset.pdf",width = 15,height = 5)
DimPlot(MyeloidSeurat.integrated,group.by="DefineTypes2",split.by = "Class",cols = c(Myeloidcolor$Colors,ClustersColor$Colors[1:6] ,ClustersColor$Colors[1:6]),label=F,label.size =2)

dev.off()
tiff("7.1Integration_NgPaper/MyeloidSeurat.integrated_SplitByDataset.tiff",width = 1800,height = 800,res=200)
DimPlot(MyeloidSeurat.integrated,group.by="DefineTypes2",split.by = "Class",
        cols = c(Myeloidcolor$Colors ,ClustersColor$Colors[1:7],ClustersColor$Colors[c(2:5,7)]),label=F,label.size =2)+NoLegend()

dev.off()


###Correlation
NgMyeloidClusterCols <- data.frame(unique(MyeloidSeurat.integrated@meta.data[,c("Class","DefineTypes2")]),
                                   Cols= c(Myeloidcolor$Colors ,ClustersColor$Colors[c(2:5,7)],ClustersColor$Colors[1:7]),
                                   ClassCols = rep(c("red","blue","darkgreen"),times=c(9,7,5))) %>% 
  dplyr::mutate(Class=rep(c("OurDataset","SMC","KUL"),times=c(9,7,5)),DefineTypes2=gsub("SMC_|KUL_","",DefineTypes2))

NgMyeloidClusterID <- lapply(split(MyeloidSeurat.integrated@meta.data,list(MyeloidSeurat.integrated@meta.data$DefineTypes2)), function(x)rownames(x))
Data <- MyeloidSeurat.integrated@assays$RNA@data %>% data.frame
#colnames(Data) <- gsub("\\.","-",colnames(Data))
NgMyeloidClusterMean <- t(apply(Data,1,function(x){
  lapply(NgMyeloidClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))

NgMyeloidClusterMean <- NgMyeloidClusterMean[apply(NgMyeloidClusterMean,1,sum)>0,]
colnames(NgMyeloidClusterMean) <- ifelse(colnames(NgMyeloidClusterMean) %in% MSCcolor$DefineTypes,paste("OurDataset",colnames(NgMyeloidClusterMean),sep="_"),colnames(NgMyeloidClusterMean))


NgMyeloidClusterCols <- NgMyeloidClusterCols[match(colnames(NgMyeloidClusterMean),paste(NgMyeloidClusterCols$Class,NgMyeloidClusterCols$DefineTypes2,sep="_")),]
SampleColA <- NgMyeloidClusterCols %>%
  dplyr::select(c(Cols,ClassCols)) %>% as.matrix() %>% set_colnames(c("DefineTypes2","Dataset")) %>%
  set_rownames(names(NgMyeloidClusterID))


TypeCor <- cor(NgMyeloidClusterMean,method="spearman")


pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/7.1Integration_NgPaper/DiffDataSet_MyeloidRNACor.pdf",width = 6,height = 7)
p <- heatmap.plus::heatmap.plus(TypeCor, ColSideColors = SampleColA,
                                col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100),
                                labRow = NA,cexCol = 1,margins = c(15,10))
#legend("topright",legend= c(unique(NgMyeloidClusterCols$Class),"", NgMyeloidClusterCols[NgMyeloidClusterCols$Class %in% c("OurDataset","SMC"),]$DefineTypes2),
#       fill=c(c("red","blue","darkgreen"),"white", NgMyeloidClusterCols[NgMyeloidClusterCols$Class %in% c("OurDataset","SMC"),]$Cols), border=FALSE, bty="n", y.intersp = 0.8, cex=0.9,xpd=T)
dev.off()


###
###Correlation
NgMyeloidClusterCols <- data.frame(unique(MyeloidSeurat.integrated@meta.data[,c("Class","DefineTypes2")]),
                          Cols= c(Myeloidcolor$Colors ,ClustersColor$Colors[c(2:5,7)],ClustersColor$Colors[1:7]),
                          ClassCols = rep(c("red","blue","darkgreen"),times=c(9,7,5))) %>% 
  dplyr::mutate(Class=rep(c("OurDataset","SMC","KUL"),times=c(9,7,5)),DefineTypes2=gsub("SMC_|KUL_","",DefineTypes2))

NgMyeloidClusterID <- lapply(split(MyeloidSeurat.integrated@meta.data,list(MyeloidSeurat.integrated@meta.data$DefineTypes2)), function(x)rownames(x))
Data <- MyeloidSeurat.integrated@assays$integrated@data %>% data.frame
#colnames(Data) <- gsub("\\.","-",colnames(Data))
IntClusterMean <- t(apply(Data,1,function(x){
  lapply(NgMyeloidClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))


IntClusterMean <- IntClusterMean[apply(IntClusterMean,1,sum)>0,]
colnames(IntClusterMean) <- ifelse(colnames(IntClusterMean) %in% MSCcolor$DefineTypes,paste("OurDataset",colnames(IntClusterMean),sep="_"),colnames(IntClusterMean))
NgMyeloidClusterCols <- NgMyeloidClusterCols[match(colnames(IntClusterMean),paste(NgMyeloidClusterCols$Class,NgMyeloidClusterCols$DefineTypes2,sep="_")),]
SampleColA <- NgMyeloidClusterCols %>%
  dplyr::select(c(Cols,ClassCols)) %>% as.matrix() %>% set_colnames(c("DefineTypes2","Dataset")) %>%
  set_rownames(names(NgMyeloidClusterID))


TypeCor <- cor(IntClusterMean,method="spearman")


pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/7.1Integration_NgPaper/IntDiffDataSet_MyeloidCor.pdf",width = 6,height = 7)
p <- heatmap.plus::heatmap.plus(TypeCor, ColSideColors = SampleColA,
                                col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100),
                                labRow = NA,cexCol = 1,margins = c(15,10))
#legend("topright",legend= c(unique(ClusterCols$Class),"", ClusterCols[ClusterCols$Class %in% c("OurDataset","SMC"),]$DefineTypes2),
#       fill=c(c("red","blue","darkgreen"),"white", ClusterCols[ClusterCols$Class %in% c("OurDataset","SMC"),]$Cols), border=FALSE, bty="n", y.intersp = 0.8, cex=0.9,xpd=T)
dev.off()



###Harmony

MyeloidSeuratList <- list(Myeloidseurat,SMC10X_Myeloid,KUL10X_Myeloid)

MyeloidHarmony <- merge(x=MyeloidSeuratList[[1]],y=MyeloidSeuratList[2:3])
#MyeloidHarmony@meta.data$Class <- ifelse(is.na(MyeloidHarmony@meta.data$Class),"Our_10X",MyeloidHarmony@meta.data$Class)

MyeloidHarmony <- NormalizeData(MyeloidHarmony)
MyeloidHarmony <- FindVariableFeatures(object =MyeloidHarmony, mean.function = ExpMean, dispersion.function = LogVMR)
MyeloidHarmony <- ScaleData(object = MyeloidHarmony)
MyeloidHarmony <- RunPCA(object = MyeloidHarmony,  npcs = 50, verbose = FALSE)

MyeloidHarmony@meta.data$Sample_Tech <- paste(MyeloidHarmony$Class,MyeloidHarmony$orig.ident,sep="_")
MyeloidHarmony@meta.data$Tissue_Tech <- paste(MyeloidHarmony$Class,MyeloidHarmony$Tissues,sep="_")
MyeloidHarmony <- MyeloidHarmony %>% 
  RunHarmony("Sample_Tech", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(MyeloidHarmony, 'harmony')
MyeloidHarmony <- MyeloidHarmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
MyeloidHarmony$DefineTypes2 <- ifelse(MyeloidHarmony$Class %in% "Our_10X",
                                      MyeloidHarmony$DefineTypes,
                                      paste(gsub("10X","",MyeloidHarmony$Class),MyeloidHarmony$DefineTypes,sep="_")
)
MyeloidHarmony$DefineTypes2 <- factor(MyeloidHarmony$DefineTypes2,
                                      levels=c(Myeloidcolor$DefineTypes,grep("SMC|KUL",names(table(MyeloidHarmony$DefineTypes2)),value = T)))
MyeloidHarmony$Class <- factor(MyeloidHarmony$Class,levels=c("Our_10X","SMC10X","KUL10X"))


readr::write_rds(MyeloidHarmony,path="7.1Integration_NgPaper/Lee_MyeloidHarmony.rds.gz",compress = "gz")
MyeloidHarmony <- readr::read_rds("7.1Integration_NgPaper/Lee_MyeloidHarmony.rds.gz")

pdf("7.1Integration_NgPaper/Lee_MyeloidHarmony_SplitByDataset.pdf",width = 15,height = 5)
DimPlot(MyeloidHarmony,group.by="DefineTypes2",split.by = "Class",cols = c(Myeloidcolor$Colors,ClustersColor$Colors[1:7] ,ClustersColor$Colors[c(2:5,7)]),label=F,label.size =2)

dev.off()
tiff("7.1Integration_NgPaper/Lee_MyeloidHarmony_SplitByDataset.tiff",width = 1800,height = 800,res=200)
DimPlot(MyeloidHarmony,group.by="DefineTypes2",split.by = "Class",
        cols = c(Myeloidcolor$Colors,ClustersColor$Colors[1:7] ,ClustersColor$Colors[c(2:5,7)]),label=F,label.size =2)+NoLegend()

dev.off()


for(i in seq(0.1,0.8,by=0.1)){
  MyeloidHarmony <- FindClusters(MyeloidHarmony, resolution =i, algorithm = 1)
}

png("7.1Integration_NgPaper/Lee_MyeloidHarmony cells with Different resolution.png",width=3200,height=1400,res=200)
DiffResolution_plot <- lapply(seq(0.1,0.8,by=0.1),function(i){
  p <- DimPlot( MyeloidHarmony , reduction = "umap",label=T,group.by = paste("RNA_snn_res.",i,sep=""),cols = .cluster_cols)+
    labs(title=paste("Resolution = ",i,sep=""))
  return(p)
})

p_res <- cowplot::plot_grid(plotlist = DiffResolution_plot,ncol=4)
print(p_res)
dev.off()




