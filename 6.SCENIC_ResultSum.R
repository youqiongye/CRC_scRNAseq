library(Seurat)
library(magrittr)
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/4.SCENIC/SCENIC_Result/")
Bdata <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/B.rds")
PlasmaBdata <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/Plasma.rds")
WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypeUpdataB.rds.gz")
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")

MetabData <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/Output/KEGG_PathwayScoreBycssmillie.rds")
GOData <-  readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/Output/RBPPathwaysScoreBycssmillie.rds.gz")
ReactData <-  readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/Output/ReactomePathwaysScoreBycssmillie.rds")


####KEGG pathways
MetabolismPathways <- gmtPathways("/home/yye/Data/Public/Genelist/KEGG_metabolism.gmt")
NonMetabolismPathways <- gmtPathways("/home/yye/Data/Public/Genelist/nonMetabolic_KEGG.gmt")
KeggPathways <- c(MetabolismPathways,NonMetabolismPathways)

scores <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/Output/KEGG_PathwayScoreBycssmillie.rds")
scoresMatrix <- scores[match(gsub("Plasma","",rownames(WholeTissueSC@meta.data)),rownames(scores)),]
scoresMatrix  <- t(scoresMatrix )
scoresMatrix <- apply(scoresMatrix,2,function(x)signif(x,digits = 3))
rownames(scoresMatrix ) <- c(names(MetabolismPathways),names(NonMetabolismPathways))
colnames(scoresMatrix) <- rownames(WholeTissueSC@meta.data)
scoresMatrixSeurat <- CreateSeuratObject(counts = scoresMatrix )
scoresMatrixSeurat@meta.data <- WholeTissueSC@meta.data
scoresMatrix1 <- scoresMatrix
rownames(scoresMatrix1 ) <- c(names(MetabolismPathways),names(NonMetabolismPathways))
colnames(scoresMatrix1) <- gsub("Epithelial_|MSC_|Endo_|Glia_|^B_|^T_|PlasmaB_|ProliferatingT_|ILCs_|Myeloid_|Mast_","",rownames(WholeTissueSC@meta.data))
scoresMatrix1Seurat <- CreateSeuratObject(counts = scoresMatrix1 )


BPPathwaysscores <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/Output/RBPPathwaysScoreBycssmillie.rds.gz")
####BP
BPPathwaysscores <- BPPathwaysscores[match(gsub("Plasma","",rownames(WholeTissueSC@meta.data)),rownames(BPPathwaysscores)),]

BPPathwaysscoresMatrix <- as.matrix(BPPathwaysscores)
BPPathwaysscoresMatrix  <- t(BPPathwaysscoresMatrix )
BPPathwaysscoresMatrix <- apply(BPPathwaysscoresMatrix,2,function(x)signif(x,digits = 3))
#rownames(BPPathwaysscoresMatrix ) <- names(BPPathways)
colnames(BPPathwaysscoresMatrix) <- rownames(WholeTissueSC@meta.data)
BPPathwaysscoresMatrixSeurat <- CreateSeuratObject(counts = BPPathwaysscoresMatrix )
BPPathwaysscoresMatrixSeurat@meta.data <- WholeTissueSC@meta.data

BPPathwaysscoresMatrix1 <- BPPathwaysscoresMatrix
#rownames(BPPathwaysscoresMatrix1 ) <- names(BPPathways)
colnames(BPPathwaysscoresMatrix1) <- gsub("Epithelial_|MSC_|Endo_|Glia_|^B_|^T_|PlasmaB_|ProliferatingT_|ILCs_|Myeloid_|Mast_","",rownames(WholeTissueSC@meta.data))
BPPathwaysscoresMatrix1Seurat <- CreateSeuratObject(counts = BPPathwaysscoresMatrix1 )
####Reactome
ReactomePathwaysscores <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/Output/ReactomePathwaysScoreBycssmillie.rds")
ReactomePathwaysscores <- ReactomePathwaysscores[match(gsub("Plasma","",rownames(WholeTissueSC@meta.data)),rownames(ReactomePathwaysscores)),]

ReactomePathwaysscoresMatrix <- as.matrix(ReactomePathwaysscores)
ReactomePathwaysscoresMatrix  <- t(ReactomePathwaysscoresMatrix )
ReactomePathwaysscoresMatrix <- apply(ReactomePathwaysscoresMatrix,2,function(x)signif(x,digits = 3))
#rownames(ReactomePathwaysscoresMatrix ) <- names(ReactomePathways)
colnames(ReactomePathwaysscoresMatrix) <- rownames(WholeTissueSC@meta.data)
ReactomePathwaysscoresMatrixSeurat <- CreateSeuratObject(counts = ReactomePathwaysscoresMatrix )
ReactomePathwaysscoresMatrixSeurat@meta.data <- WholeTissueSC@meta.data
ReactomePathwaysscoresMatrix1 <- ReactomePathwaysscoresMatrix
#rownames(ReactomePathwaysscoresMatrix1 ) <-names(ReactomePathways)
colnames(ReactomePathwaysscoresMatrix1) <- gsub("Epithelial_|MSC_|Endo_|Glia_|^B_|^T_|PlasmaB_|ProliferatingT_|ILCs_|Myeloid_|Mast_","",rownames(WholeTissueSC@meta.data))
ReactomePathwaysscoresMatrix1Seurat <- CreateSeuratObject(counts = ReactomePathwaysscoresMatrix1 )

WholeTissueList1 <- lapply(names(WholeTissueList),function(x){
  sub <- WholeTissueList[[x]]
  subScore <- subset(scoresMatrix1Seurat,cells=rownames(sub@meta.data))
  subBPPathwaysscores <- subset(BPPathwaysscoresMatrix1Seurat,cells=rownames(sub@meta.data))
  subReactomePathwaysscores <- subset(ReactomePathwaysscoresMatrix1Seurat,cells=rownames(sub@meta.data))
  
  sub@assays$Metab <- subScore@assays$RNA
  sub@assays$GO <- subBPPathwaysscores@assays$RNA
  sub@assays$React <- subReactomePathwaysscores@assays$RNA
  print(x)
  return(sub)
  
})
names(WholeTissueList1) <- names(WholeTissueList)
###TF AUC
AUClist <- list.files(path="AUC/",pattern="*AUC.csv")
AUCData <- lapply(AUClist,function(x){
  AUCvalue <- read.csv(paste("AUC/",x,sep=""))
  AUCvalueMatrix <- as.matrix(AUCvalue[,2:ncol(AUCvalue)])
  AUCvalueMatrix  <- t(AUCvalueMatrix )
  AUCvalueMatrix <- apply(AUCvalueMatrix,2,function(x)signif(x,digits = 3))
  #rownames(AUCvalueMatrix ) <- names(ReactomePathways)
  colnames(AUCvalueMatrix) <- AUCvalue$Cell
  AUCvalueMatrixSeurat <- CreateSeuratObject(counts = AUCvalueMatrix )
  AUCvalueMatrixSeurat@meta.data <- WholeTissueSC@meta.data
  
  return(AUCvalueMatrixSeurat)
})
names(AUCData) <- c("Endo","Epithelial","Myeloid","MSC","T")

for(i in names(AUCData)){
  WholeTissueList1[[i]]@assays$TF <- AUCData[[i]]@assays$RNA
}
WholeTissueList1[["B"]]@assays$TF <- Bdata@assays$TF
WholeTissueList1[["PlasmaB"]]@assays$TF <- PlasmaBdata@assays$TF

readr::write_rds(WholeTissueList1,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtype_AllFunctionTF.rds.gz",compress = "gz")

WholeTissueSC@assays$Metab <- scoresMatrixSeurat@assays$RNA
WholeTissueSC@assays$BP <- BPPathwaysscoresMatrixSeurat@assays$RNA
WholeTissueSC@assays$React <- ReactomePathwaysscoresMatrixSeurat@assays$RNA
readr::write_rds(WholeTissueSC,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineSubtype_AllFunctionTF.rds.gz",compress = "gz")

CellList <- SplitObject(WholeTissueSC, split.by = "DefineTypes")
CellList1 <- lapply(CellList,function(sub){
  Idents(sub) <- sub$Tissues
  LL <- lapply(c("RNA","Metab","BP","React"),function(x){
    DefaultAssay(sub) <- x#"Metab"
    DiffFeature <- FindAllMarkers(sub,min.pct = 0.025,logfc.threshold = 0.1)
    return (DiffFeature)
  })
  return(LL)
})

CellListMAST <- lapply(CellList,function(sub){
  Idents(sub) <- sub$Tissues
  LL <- lapply(c("RNA","Metab","BP","React"),function(x){
    DefaultAssay(sub) <- x#"Metab"
    DiffFeature <- FindAllMarkers(sub,min.pct = 0.025,logfc.threshold = 0.1,test.use="MAST")
    return (DiffFeature)
  })
  return(LL)
})
readr::write_rds(CellListMAST,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.NTDiff/NT_DifferenceMast.rds.gz",compress = "gz")
CellListMAST$`FAP+ fibroblasts`[[3]] %>% summary()
CellListMAST1 <- CellListMAST
CellListMAST <- lapply(CellListMAST,function(sub){
  LL <- lapply(sub,function(x){
    DiffFeature <- x[!duplicated(x$gene),]
    return (DiffFeature)
  })
  names(LL) <- c("RNA","Metab","BP","React")
  return(LL)
})







FeaturePlot(WholeTissueList1$Myeloid,feature="IL3RA")


FeaturePlot(WholeTissueList1$T,feature="DEPDC5")
Idents(WholeTissueSC) <- WholeTissueSC$SubTypes
FeaturePlot(WholeTissueSC,feature="DEPDC5",label=T)

VlnPlot(WholeTissueSC,feature="DEPDC5")

DEPDC5 <- data.frame(DEPDC5=as.numeric(WholeTissueSC@assays$RNA@data[rownames(WholeTissueSC@assays$RNA@data) %in% "DEPDC5",]),Subtypes=WholeTissueSC$SubTypes)
DEPDC5$Class <- ifelse(DEPDC5$DEPDC5 > 0,1,0)


options(stringsAsFactors = F)
library(ggplot2)
library(magrittr)

WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypeUpdataB.rds.gz")
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")
Idents(WholeTissueSC) <- WholeTissueSC$SubTypes

MSCseurat <- WholeTissueList$MSC
ClustersColor <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
MSCcolor <- ClustersColor[ClustersColor$SubTypes %in% "MSC",]
###All samples distribution.


MSC_RegulonAUC <- read.csv("Human_Stroma_AUC.csv")
colnames(MSC_RegulonAUC) <- gsub("\\.","",colnames(MSC_RegulonAUC))

#MSC_Regulon <- read.csv("Human_Stroma_Reg.csv")[,1:10]

MSC_RegulonAUCmatrix <- as.matrix(MSC_RegulonAUC[,2:ncol(MSC_RegulonAUC)]) %>% t
rownames(MSC_RegulonAUCmatrix ) <- colnames(MSC_RegulonAUC)[2:ncol(MSC_RegulonAUC)]
colnames(MSC_RegulonAUCmatrix) <- rownames(MSCseurat@meta.data)
MSC_RegulonAUCmatrixSeurat <- CreateSeuratObject(counts = MSC_RegulonAUCmatrix )


MSCseurat@assays$RegulonAUC <- MSC_RegulonAUCmatrixSeurat@assays$RNA


MSC_Regulon_ClusterID <- lapply(split(MSCseurat@meta.data,list(MSCseurat$DefineTypes)), function(x)rownames(x))


MSC_RegulonClusterMean <- t(apply(MSC_RegulonAUCmatrix,1,function(x){
  lapply(MSC_Regulon_ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))

readr::write_rds(MSC_RegulonClusterMean,path="MSC_RegulonClusterMean.rds")


MSC_RegulonClusterMean <- MSC_RegulonClusterMean[apply(MSC_RegulonClusterMean,1,sum)>0,]

MSC_RegulonClusterMeanM <- MSC_RegulonClusterMean
MSC_RegulonClusterMeanM[,1:ncol(MSC_RegulonClusterMeanM)] <- t(apply(MSC_RegulonClusterMeanM,1,scale))
Sample_order <- Clusters


MSC_RegulonClusterMeanM[MSC_RegulonClusterMeanM< -3] <- -3
MSC_RegulonClusterMeanM[MSC_RegulonClusterMeanM> 3] <- 3


Wholetissue_right_ha = ComplexHeatmap::rowAnnotation(foo =  ComplexHeatmap::anno_mark(at = match(PlotGenes,WholeTissueSC_ClusterMarkerFC0.5M$gene), labels = PlotGenes ))



Wholetissue_barcolor <- Clusters$Colors
names(Wholetissue_barcolor) <- Clusters$DefineTypes
Wholetissue_SampleBar <-  Sample_order$DefineTypes
names(Wholetissue_SampleBar) <-  Clusters$DefineTypes

Wholetissue_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Wholetissue_SampleBar,
                                                        col = list(bar = Wholetissue_barcolor,Main))


pdf("Final_Wholetissue_Heatmap.pdf",width = 8,height = 8)
ComplexHeatmap::Heatmap(MSC_RegulonClusterMeanM,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,column_order= names(Wholetissue_SampleBar),use_raster = TRUE,  name = "mat",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),#colorRampPalette(c("blue","white","red"),space="rgb")(100),#colorRampPalette(c("magenta", "black", "yellow"))(256),#
                        show_row_names = F,show_column_names = F,right_annotation = Wholetissue_right_ha,top_annotation = Wholetissue_ha_bar)
dev.off()


ComplexHeatmap::Heatmap(MSC_RegulonClusterMeanM,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,  name = "mat",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),#colorRampPalette(c("blue","white","red"),space="rgb")(100),#colorRampPalette(c("magenta", "black", "yellow"))(256),#
                        show_row_names = F,show_column_names = F)
pheatmap::pheatmap(MSC_RegulonClusterMeanM)
