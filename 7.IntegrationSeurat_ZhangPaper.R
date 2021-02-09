library(Seurat)
library(ggplot2)
library(magrittr)
library(Matrix)
library(wordspace)
library(Hmisc)
library(rsvd)
options(stringsAsFactors = F)
set.seed(42)
setwd("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/7.Integration_ZhangPaper")
Samples <- c("p7_N","p9_N", "p10_N", "p11_N", "p12_N", "p7_T",  "p9_T",  "p10_T", "p11_T", "p12_T")


.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
Leu10X_Seurat<- readr::read_rds("/home/yye/Data/Public/GEO_data/scRNAseq/ZeminZhang_CRC/Leu10X_Seurat.rds.gz")
LeuSMART_Seurat <- readr::read_rds("/home/yye/Data/Public/GEO_data/scRNAseq/ZeminZhang_CRC/LeuSMART_Seurat.rds.gz")
LeuSMART_Myeloid <- subset(LeuSMART_Seurat,cells=rownames(LeuSMART_Seurat@meta.data[LeuSMART_Seurat$Global_Cluster %in% "Myeloid cell" & LeuSMART_Seurat$Sub_Cluster !="hM01_Mast-TPSAB1",]))
Leu10X_Myeloid <- subset(Leu10X_Seurat,cells=rownames(Leu10X_Seurat@meta.data[Leu10X_Seurat$Global_Cluster %in% "Myeloid cell" & Leu10X_Seurat$Sub_Cluster !="hM01_Mast-TPSAB1",]))
Leu10X_Myeloid <- FindVariableFeatures(Leu10X_Myeloid)
Leu10X_Myeloid <- ScaleData(Leu10X_Myeloid)
Leu10X_Myeloid <- RunPCA(Leu10X_Myeloid)

ClustersColor <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
Myeloidcolor <- ClustersColor[ClustersColor$SubTypes %in% "Myeloid",]
ggplot(LeuSMART_Myeloid@meta.data,aes(x=Sub_tSNE_1, y= Sub_tSNE_2,color=Sub_Cluster))+
  geom_point()
ggplot(Leu10X_Myeloid@meta.data,aes(x=Sub_tSNE_1, y= Sub_tSNE_2,color=Sub_Cluster))+
  geom_point()


WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtype_AllFunctionTF.rds.gz")
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")
p1 <- FeaturePlot(WholeTissueSC,feature="GJB2")
p2 <- FeaturePlot(WholeTissueList$MSC,feature="GJB2")

cowplot::plot_grid(plotlist = list(p1,p2))


Myeloidseurat <- WholeTissueList$Myeloid
Myeloidseurat@meta.data <- Myeloidseurat@meta.data[,c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes")]
Myeloidseurat@meta.data$Class <- rep("Our_10X",times=nrow(Myeloidseurat@meta.data))
LeuSMART_Myeloid@meta.data <- LeuSMART_Myeloid@meta.data[,c("Sample","nCount_RNA", "nFeature_RNA","Tissue","Sub_Cluster")] %>% set_colnames(c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes"))
LeuSMART_Myeloid@meta.data$Class <- rep("Zhang_SMART",times=nrow(LeuSMART_Myeloid@meta.data))
Leu10X_Myeloid@meta.data <- Leu10X_Myeloid@meta.data[,c("Sample","nCount_RNA", "nFeature_RNA","Tissue","Sub_Cluster")]%>% set_colnames(c("orig.ident","nCount_RNA", "nFeature_RNA","Tissues","DefineTypes"))
Leu10X_Myeloid@meta.data$Class <- rep("Zhang_10X",times=nrow(Leu10X_Myeloid@meta.data))

###Integration 10 and SMARTseq
SeuratList <- list(Myeloidseurat,Leu10X_Myeloid)
names(SeuratList) <- c("Our_10X","Zhang_10X")

Idents(Myeloidseurat) <- "DefineTypes"
FeaturePlot(Myeloidseurat,feature="FOLR2",cols=c("gray","orange","red"),label=T,label.size = 3.5)


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
Seurat.integrated$DefineTypes <- factor(Seurat.integrated$DefineTypes,levels=c(Myeloidcolor$DefineTypes,grep("hM",names(table(Seurat.integrated$DefineTypes)),value = T)))
readr::write_rds(Seurat.integrated,path="Seurat.integrated.rds.gz",compress="gz")

Seurat.integrated <- readr::read_rds("Seurat.integrated.rds.gz")
pdf("Seurat.integrated_SplitByDataset.pdf",width = 15,height = 4)
DimPlot(Seurat.integrated,group.by="DefineTypes",split.by = "Class",cols = c(Myeloidcolor$Colors,ClustersColor$Colors[1:12]),label=F,label.size =2)
dev.off()

tiff("Seurat.integrated_SplitByDataset.tiff",width = 3000,height = 1000,res=250)
DimPlot(Seurat.integrated,group.by="DefineTypes",split.by = "Class",cols = c(Myeloidcolor$Colors,ClustersColor$Colors[c(1:11,13)]),label=F,label.size =2)
dev.off()


Seurat.integrated <- FindNeighbors(Seurat.integrated)
for(i in seq(0.1,1,by=0.1)){
  Seurat.integrated <- FindClusters(Seurat.integrated, resolution = i, algorithm = 1)
}

png("Seurat.integrated with Different resolution.png",width=2000,height=1000,res=150)
DiffResolution_plot <- lapply(seq(0.1,0.6,by=0.1),function(i){
  p <- DimPlot(   Seurat.integrated, reduction = "umap",label=T,group.by = paste("integrated_snn_res.",i,sep=""),cols = .cluster_cols)+
    labs(title=paste("Resolution = ",i,sep=""))
  return(p)
})
p_res <- cowplot::plot_grid(plotlist = DiffResolution_plot,ncol=3)
print(p_res)
dev.off()
DefaultAssay(Seurat.integrated) <- "RNA"
FeaturePlot(Seurat.integrated,feature="CLEC9A",split.by = "Class")


CRC.query <- SeuratList[["Our_10X"]]
CRC.anchors <- FindTransferAnchors(reference = Leu10X_Myeloid, query = CRC.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = CRC.anchors, refdata = Leu10X_Myeloid$DefineTypes, 
                            dims = 1:30)
CRC.query <- AddMetaData(CRC.query, metadata = predictions)
CRC.query$prediction.match <- CRC.query$predicted.id == CRC.query$DefineTypes
table(CRC.query$prediction.match)

CRC.query$predicted.idCut <- ifelse(CRC.query$prediction.score.max >0.5,as.character(CRC.query$DefineTypes),CRC.query$predicted.id)



###Correlation
MyeClusterCols <- data.frame(unique(Seurat.integrated@meta.data[,c("Class","DefineTypes")]),
                          Cols= c(Myeloidcolor$Colors,ClustersColor$Colors[c(1:11,13)]),
                          ClassCols = rep(c("red","darkgreen"),times=c(9,12))) %>% 
  dplyr::mutate(Class=rep(c("OurDataset","Zhang"),times=c(9,12)))

ClusterID <- lapply(split(Seurat.integrated@meta.data,list(Seurat.integrated@meta.data$DefineTypes)), function(x)rownames(x))
Data <- Seurat.integrated@assays$RNA@data %>% data.frame
#colnames(Data) <- gsub("\\.","-",colnames(Data))
ClusterMean <- t(apply(Data,1,function(x){
  lapply(ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))

ClusterMean <- ClusterMean[apply(ClusterMean,1,sum)>0,]
colnames(ClusterMean) <- ifelse(colnames(ClusterMean) %in% MSCcolor$DefineTypes,paste("OurDataset",colnames(ClusterMean),sep="_"),colnames(ClusterMean))

SampleColA <- MyeClusterCols %>%
  dplyr::select(c(Cols,ClassCols)) %>% as.matrix() %>% set_colnames(c("DefineTypes2","Dataset")) %>%
  set_rownames(names(ClusterID))


TypeCor <- cor(ClusterMean,method="spearman")


pdf("DiffDataSet_MyeloidCor_RNA.pdf",width = 6,height = 7)
p <- heatmap.plus::heatmap.plus(TypeCor, ColSideColors = SampleColA,
                                col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100),
                                labRow = NA,cexCol = 1,margins = c(15,10))
#legend("topright",legend= c(unique(MyeClusterCols$Class),"", MyeClusterCols[MyeClusterCols$Class %in% c("OurDataset","SMC"),]$DefineTypes2),
#       fill=c(c("red","blue","darkgreen"),"white", MyeClusterCols[MyeClusterCols$Class %in% c("OurDataset","SMC"),]$Cols), border=FALSE, bty="n", y.intersp = 0.8, cex=0.9,xpd=T)
dev.off()


###
###Correlation

ClusterID <- lapply(split(Seurat.integrated@meta.data,list(Seurat.integrated@meta.data$DefineTypes)), function(x)rownames(x))
Data <- Seurat.integrated@assays$integrated@data %>% data.frame
#colnames(Data) <- gsub("\\.","-",colnames(Data))
IntClusterMean <- t(apply(Data,1,function(x){
  lapply(ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))


IntClusterMean <- IntClusterMean[apply(IntClusterMean,1,sum)>0,]
colnames(IntClusterMean) <- ifelse(colnames(IntClusterMean) %in% MSCcolor$DefineTypes,paste("OurDataset",colnames(IntClusterMean),sep="_"),colnames(IntClusterMean))

SampleColA <- MyeClusterCols %>%
  dplyr::select(c(Cols,ClassCols)) %>% as.matrix() %>% set_colnames(c("DefineTypes","Dataset")) %>%
  set_rownames(names(ClusterID))


TypeCor <- cor(IntClusterMean,method="spearman")


pdf("IntDiffDataSet_MyeloidCor.pdf",width = 6,height = 7)
p <- heatmap.plus::heatmap.plus(TypeCor, ColSideColors = SampleColA,
                                col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:10],space="rgb")(100),
                                labRow = NA,cexCol = 1,margins = c(15,10))
#legend("topright",legend= c(unique(MyeClusterCols$Class),"", MyeClusterCols[MyeClusterCols$Class %in% c("OurDataset","SMC"),]$DefineTypes2),
#       fill=c(c("red","blue","darkgreen"),"white", MyeClusterCols[MyeClusterCols$Class %in% c("OurDataset","SMC"),]$Cols), border=FALSE, bty="n", y.intersp = 0.8, cex=0.9,xpd=T)
dev.off()





