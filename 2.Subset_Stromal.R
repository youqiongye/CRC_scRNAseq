library(Seurat)
library(ggplot2)
library(magrittr)
library(Matrix)
library(wordspace)
library(Hmisc)
library(rsvd)
options(stringsAsFactors = F)
set.seed(42)
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/Stromal/")
Samples <- c("p7_N","p9_N", "p10_N", "p11_N", "p12_N", "p7_T",  "p9_T",  "p10_T", "p11_T", "p12_T")
###Features
##1. mast cells,  Expression profiling of constitutive mast cells reveals a unique identity within the immune system
##Mast: KIT (CD117) ,IL1RL1, Mast cells and basophils: trojan horses of conventional lin- stem/progenitor cell isolates
##Mast activated marker: ENPP3


.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
TissueStromal <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/StromalCells_harmony.rds.gz")

TissueStromal <- NormalizeData(TissueStromal)
TissueStromal <- FindVariableFeatures(object =TissueStromal, mean.function = ExpMean, dispersion.function = LogVMR)
TissueStromal <- ScaleData(object = TissueStromal)
TissueStromal <- RunPCA(object = TissueStromal,  npcs = 30, verbose = FALSE)

TissueStromal <- TissueStromal %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(TissueStromal, 'harmony')
TissueStromal <- TissueStromal %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) 

TissueStromal <- TissueStromal %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "umap", dims = 1:2)
for(i in seq(0.05,0.3,by=0.05)){
  TissueStromal <- FindClusters(TissueStromal, resolution =i, algorithm = 1)
}
TissueStromal@meta.data <- TissueStromal@meta.data[,c(1:6,grep("RNA_snn",colnames(TissueStromal@meta.data)))]

png("1.Cluster2 with Different resolution.png",width=2400,height=1400,res=200)
DiffResolution_plot <- lapply(seq(0.05,0.3,by=0.05),function(i){
  p <- DimPlot( TissueStromal, reduction = "umap",label=T,group.by = paste("RNA_snn_res.",i,sep=""),cols = .cluster_cols)+
    labs(title=paste("Resolution = ",i,sep=""))
  return(p)
})
p_res <- cowplot::plot_grid(plotlist = DiffResolution_plot,ncol=3)
print(p_res)
dev.off()



ref.se <- readr::read_rds("/home/yye/Project/Collaboration/2020ShenLei/Matrix/ProcessData/SingeR_ref.se.rds.gz")
bp.se <- ref.se$BP
hpca.se <- ref.se$HPCA

TissueStromal_SingleR_Main <- SingleR::SingleR(GetAssayData(TissueStromal@assays$RNA, slot = "data"),
                                               ref=ref.se,
                                               method = "cluster",
                                               clusters =TissueStromal@meta.data$RNA_snn_res.0.05,
                                               labels = list(bp.se$label.main, hpca.se$label.main))
TissueStromalSingleR_Result <- data.frame(HPCA =TissueStromal_SingleR_Main@listData$orig.results$HPCA$first.labels,
                             BP =TissueStromal_SingleR_Main@listData$orig.results$BP$first.labels,
                             clusters=sort(unique(TissueStromal@meta.data$RNA_snn_res.0.05)))

TissueStromal@meta.data$HPCA <- TissueStromalSingleR_Result$HPCA[match(TissueStromal@meta.data$RNA_snn_res.0.05,TissueStromalSingleR_Result$clusters)]
TissueStromal@meta.data$BP <- TissueStromalSingleR_Result$BP[match(TissueStromal@meta.data$RNA_snn_res.0.05,TissueStromalSingleR_Result$clusters)]

DimPlot(TissueStromal,reduction = "umap",group.by="HPCA",cols = .cluster_cols)
DimPlot(TissueStromal,reduction = "umap",group.by="BP",cols = .cluster_cols)
pdf("Feature_Count_mito.pdf",width = 6,height = 9)
VlnPlot(object = TissueStromal, features = c("nCount_RNA","nFeature_RNA","percent.mito"),group.by="RNA_snn_res.0.05",pt.size = 0, cols= .cluster_cols,ncol = 1)+NoLegend()
 dev.off() 
FeaturePlot(TissueStromal,reduction = "umap",features = c("CD3E","CAFs",""),cols = c("gray","orange","red"))

###CAFs
CAFs <- c("FAP", "THY1", "DCN", "COL1A1", "COL1A2", "COL6A1", "COL6A2", "COL6A3")
TissueStromal@meta.data$CAFs <- apply(TissueStromal@assays$RNA@data[rownames(TissueStromal@assays$RNA@data) %in% CAFs,],2,mean)
pdf("CAFs_and markers.pdf",width = 9,height = 9)
p<-FeaturePlot(TissueStromal,reduction = "umap",features = c("CAFs",CAFs ),cols = c("gray","orange","red"),combine = F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+theme(axis.title = element_blank(),axis.text = element_blank(),
                                      axis.ticks = element_blank(),axis.line = element_line(color="gray",linetype = "dashed"))
}
cowplot::plot_grid(plotlist = p,ncol = 3)
dev.off()



panCAFs <- c("VIM", "COL1A1", "COL3A1")
TissueStromal@meta.data$panCAFs <- apply(TissueStromal@assays$RNA@data[rownames(TissueStromal@assays$RNA@data) %in% panCAFs,],2,mean)
FeaturePlot(TissueStromal,reduction = "umap",features = "panCAFs",cols = c("gray","orange","red"))
pdf("PanCAFs markers.pdf",width = 9,height = 3)
p<- FeaturePlot(TissueStromal,reduction = "umap",features = panCAFs,cols = c("gray","orange","red"),combine = F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+theme(axis.title = element_blank(),axis.text = element_blank(),
                                      axis.ticks = element_blank(),axis.line = element_line(color="gray",linetype = "dashed"))
}
cowplot::plot_grid(plotlist = p,ncol = 3)
dev.off()
Adipocytes <- c("ADH1B",	"DLAT",	"GPD1",	"LBP",	"PLIN1",	"PPP1R1A","PPP2R1B",	"PTGER3",	"CILP",	"ADIPOQ",	"COL5A3",	"PNPLA2")
TissueStromal@meta.data$Adipocytes <- apply(TissueStromal@assays$RNA@data[rownames(TissueStromal@assays$RNA@data) %in% Adipocytes,],2,mean)
FeaturePlot(TissueStromal,reduction = "umap",features = "Adipocytes",cols = c("gray","orange","red"))
FeaturePlot(TissueStromal,reduction = "umap",features = "PPP1R1A",cols = c("gray","orange","red"))

### myofibroblasts(MFs) citation: James Kinchen et al, Structural Remodeling of the Human Colonic Mesenchyme in Inflammatory Bowel Disease
MFs <- c("MYH11","ACTG2","MYOCD")
TissueStromal@meta.data$MFs <- apply(TissueStromal@assays$RNA@data[rownames(TissueStromal@assays$RNA@data) %in%MFs,],2,mean)

FeaturePlot(TissueStromal,reduction = "umap",features = "MFs",cols = c("gray","orange","red"))

FeaturePlot(TissueStromal,reduction = "umap",features = MFs,cols = c("gray","orange","red"))
pdf("MF markers.pdf",width = 8,height = 3)
p<- FeaturePlot(TissueStromal,reduction = "umap",features = MFs,cols = c("gray","orange","red"),combine = F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+theme(axis.title = element_blank(),axis.text = element_blank(),
                                      axis.ticks = element_blank(),axis.line = element_line(color="gray",linetype = "dashed"))
}
cowplot::plot_grid(plotlist = p,ncol = 3)
dev.off()
##Glial S100b, citation: James Kinchen et al, Structural Remodeling of the Human Colonic Mesenchyme in Inflammatory Bowel Disease

Glial <- c("S100B")
pdf("glial_marker_c8.pdf",width = 3,height = 3)
FeaturePlot(TissueStromal,reduction = "umap",features = "S100B",cols = c("gray","orange","red"))+NoLegend()
dev.off()
Pericytes <- c("RGS5")
pdf("Pericytes marker_c8.pdf",width = 3,height = 3)
FeaturePlot(TissueStromal,reduction = "umap",features = "RGS5",cols = c("gray","orange","red"))+NoLegend()

dev.off()

EndothelialCells <- c("PECAM1","VWF","CDH5")
TissueStromal@meta.data$EndothelialCells <- apply(TissueStromal@assays$RNA@data[rownames(TissueStromal@assays$RNA@data) %in% EndothelialCells,],2,mean)
FeaturePlot(TissueStromal,reduction = "umap",features = "EndothelialCells",cols = c("gray","orange","red"))
##mv endothelial cells: PLVAP
## endothelial WWRN1+
FeaturePlot(TissueStromal,reduction = "umap",features = "PLVAP",cols = c("gray","orange","red"))
FeaturePlot(TissueStromal,reduction = "umap",features = "ACKR1",cols = c("gray","orange","red"))
FeaturePlot(TissueStromal,reduction = "umap",features = "MMRN1",cols = c("gray","orange","red"))
pdf("Endothelial markers.pdf",width = 9,height =6)
p<- FeaturePlot(TissueStromal,reduction = "umap",features = c("PECAM1","VWF","CDH5","PLVAP","ACKR1","MMRN1"),cols = c("gray","orange","red"),combine = F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+theme(axis.title = element_blank(),axis.text = element_blank(),
                                      axis.ticks = element_blank(),axis.line = element_line(color="gray",linetype = "dashed"))
}






p <-FeaturePlot(ColonData,reduction = "umap",features = EndothelialCells,cols = c("gray","orange","orange","red","red"),combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
