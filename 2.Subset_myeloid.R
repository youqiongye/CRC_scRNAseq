library(Seurat)
library(ggplot2)
library(magrittr)
library(Matrix)
library(wordspace)
library(Hmisc)
library(rsvd)
options(stringsAsFactors = F)
set.seed(42)
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/Myeloid/")
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
TissueMyeloid <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/MyeloidCells_harmony.rds.gz")
#Idents(TissueMyeloid) <- TissueMyeloid$Pat_Tissues

TissueMyeloid <- NormalizeData(TissueMyeloid)
TissueMyeloid <- FindVariableFeatures(object =TissueMyeloid, dispersion.function = LogVMR)
TissueMyeloid <- ScaleData(object = TissueMyeloid)
TissueMyeloid <- RunPCA(object = TissueMyeloid,  npcs = 30, verbose = FALSE)

TissueMyeloid <- TissueMyeloid %>% 
  RunHarmony("Pat_Tissues", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(TissueMyeloid, 'harmony')
TissueMyeloid <- TissueMyeloid %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
TissueMyeloid <- TissueMyeloid %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) 
for(i in seq(0.1,0.8,by=0.1)){
  TissueMyeloid <- FindClusters(TissueMyeloid, resolution =i, algorithm = 1)
}
TissueMyeloid@meta.data <- TissueMyeloid@meta.data[,c(1:6,grep("RNA_snn",colnames(TissueMyeloid@meta.data)))]

png("1.Cluster with Different resolution.png",width=3200,height=1400,res=200)
DiffResolution_plot <- lapply(seq(0.1,0.8,by=0.1),function(i){
  p <- DimPlot( TissueMyeloid, reduction = "umap",label=T,group.by = paste("RNA_snn_res.",i,sep=""),cols = .cluster_cols)+
    labs(title=paste("Resolution = ",i,sep=""))
  return(p)
})
p_res <- cowplot::plot_grid(plotlist = DiffResolution_plot,ncol=4)
print(p_res)
dev.off()


DimPlot( TissueMyeloid, reduction = "umap",label=T,group.by = paste("RNA_snn_res.",0.8,sep=""),cols = .cluster_cols)
  

FeaturePlot(TissueMyeloid,features = c("ENPP3","IL1RL1","KIT"))
FeaturePlot(tt,features = c("CLEC4C","CLEC9A"))

tt <- readr::read_rds("/home/jjingq/WholeTissue/myeloid/Myeloid_tmp")
DimPlot(TissueMyeloid,group.by = "RNA_snn_res.0.8")
ttSingleR_Main <- SingleR::SingleR(GetAssayData(tt@assays$RNA, slot = "data"),
                                               ref=ref.se,
                                               method = "cluster",
                                               clusters =tt@meta.data$RNA_snn_res.0.8,
                                               labels = list(bp.se$label.main, hpca.se$label.main))
ttSingleR_Result <- data.frame(HPCA =ttSingleR_Main@listData$orig.results$HPCA$first.labels,
                                          BP =ttSingleR_Main@listData$orig.results$BP$first.labels,
                                          clusters=sort(unique(tt@meta.data$RNA_snn_res.0.8)))

ref.se <- readr::read_rds("/home/yye/Project/Collaboration/2020ShenLei/Matrix/ProcessData/SingeR_ref.se.rds.gz")
bp.se <- ref.se$BP
hpca.se <- ref.se$HPCA

TissueMyeloid_SingleR_Main <- SingleR::SingleR(GetAssayData(TissueMyeloid@assays$RNA, slot = "data"),
                                            ref=ref.se,
                                            method = "cluster",
                                            clusters =TissueMyeloid@meta.data$RNA_snn_res.0.8,
                                            labels = list(bp.se$label.main, hpca.se$label.main))
TissueMyeloidSingleR_Result <- data.frame(HPCA =TissueMyeloid_SingleR_Main@listData$orig.results$HPCA$first.labels,
                             BP =TissueMyeloid_SingleR_Main@listData$orig.results$BP$first.labels,
                             clusters=sort(unique(TissueMyeloid@meta.data$RNA_snn_res.0.8)))

TissueMyeloid@meta.data$HPCA <- SingleR_Result$HPCA[match(TissueMyeloid@meta.data$RNA_snn_res.0.8,SingleR_Result$clusters)]
TissueMyeloid@meta.data$BP <- SingleR_Result$BP[match(TissueMyeloid@meta.data$RNA_snn_res.0.8,SingleR_Result$clusters)]



####
MyeloidFake <- subset(TissueMyeloid,cells=rownames(TissueMyeloid@meta.data[TissueMyeloid@meta.data$RNA_snn_res.0.8 %in% 8,]))
readr::write_rds(MyeloidFake,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/MyeloidFake_Treal.rds.gz",compress = "gz")

MyeloidUpdate <- subset(TissueMyeloid,cells=rownames(TissueMyeloid@meta.data[!(TissueMyeloid@meta.data$RNA_snn_res.0.8 %in% 8),]))
readr::write_rds(MyeloidUpdate,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/MyeloidReal.rds.gz",compress = "gz")



###Macrophage
MacrophageMarker <- c("CD14", "FCGR3A", "FCGR1A", "CD68", "TFRC", "CCR5") 

###Mast: KIT (CD117) ,IL1RL1, Mast cells and basophils: trojan horses of conventional lin- stem/progenitor cell isolates
###Mast activated marker: ENPP3
###Neutrophil:FCGR3B (CD16),CXCL8 
##macrophages "C1Q" ##citation: Dysfunctional CD8 T Cells Form a Proliferative, Dynamically Regulated Compartment within Human Melanoma
# cDC1: c("IRF8","ID2","BATF3","NFIL3")
# cDC2: c("IRF4","KLF4","ZEB2")

##monocytes: "VCAN",
##DCs: "CLEC10A","CD1C"
##macrophages "C1Q" ##citation: Dysfunctional CD8 T Cells Form a Proliferative, Dynamically Regulated Compartment within Human Melanoma
###M2: "MRC1","CD163"

FeaturePlot(TissueMyeloid,features = c("CD14"))
