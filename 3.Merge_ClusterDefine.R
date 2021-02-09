library(Seurat)
library(ggplot2)
library(magrittr)
library(harmony)
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue")
SubTypesClusters <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/SubTypesColorCode.xlsx")
Clusters <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")

TissuePlasmaB <- readr::read_rds("/home/jjingq/WholeTissue/BCell/PlasmaB/4thTissuePlasmacells.rds.gz")
TissuePlasmaB@meta.data <- TissuePlasmaB@meta.data[,c(1:6,ncol(TissuePlasmaB@meta.data))]
TissueB <- readr::read_rds("/home/jjingq/WholeTissue/BCell/BCellOnly/2ndTissueBcells.rds.gz")
TissueB@meta.data <- TissueB@meta.data[,c(1:6,ncol(TissueB@meta.data))]
ListB <- list(TissueB,TissuePlasmaB)
names(ListB) <- c("B","PlasmaB")

WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypes.rds.gz")
WholeTissueList$Mast@meta.data <- WholeTissueList$Mast@meta.data[,1:7]

WholeTissueListNew <- c(WholeTissueList[c("Epithelial","MSC","Endo","Glia")],ListB ,WholeTissueList[c("T","ProliferatingT","ILCs","Myeloid","Mast")])

readr::write_rds(WholeTissueListNew,path="~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypeUpdataB.rds.gz",compress = "gz")
openxlsx::write.xlsx(data.frame(DefineTypes=names(table(WholeTissueListNew$B$DefineTypes)),Colors=scales::hue_pal()(7)),file="tmp.xlsx")

ClustersColor <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
WholeTissueSC <- merge(x = WholeTissueListNew[[1]], y = WholeTissueListNew[names(WholeTissueListNew)[2:length(WholeTissueListNew)]],add.cell.ids=c(names(WholeTissueListNew)),merge.data=T,project = "SeuratProject")


ClusterNames <- lapply(WholeTissueListNew,function(x)names(table(x@meta.data$DefineTypes)))
MainTypes <- data.frame(Subtypes=names(WholeTissueListNew),MainTypes=c("Epithelial","Stroma","Stroma","Stroma","B","B","T","T","T","Myeloid","Myeloid"))
Clusters <- data.frame(DefineTypes=unlist(ClusterNames),SubTypes=rep(names(ClusterNames),times=unlist(lapply(ClusterNames,length))))

Clusters$MainTypes <- MainTypes$MainTypes[match(Clusters$SubTypes,MainTypes$Subtypes)]

WholeTissueSC@meta.data$SubTypes <- Clusters$SubTypes[match(WholeTissueSC@meta.data$DefineTypes,Clusters$DefineTypes)]
WholeTissueSC@meta.data$MainTypes <- MainTypes$MainTypes[match(WholeTissueSC@meta.data$SubTypes ,MainTypes$Subtypes)]
WholeTissueSC@meta.data$TopTypes <- ifelse(WholeTissueSC@meta.data$MainTypes %in% c("B","T","Myeloid"),"Immune",as.character(WholeTissueSC$MainTypes))
WholeTissueSC@meta.data$DefineTypes <- factor(WholeTissueSC@meta.data$DefineTypes,levels=Clusters$DefineTypes)
WholeTissueSC@meta.data$SubTypes <- ifelse(WholeTissueSC@meta.data$SubTypes %in% c("ILCs","ProliferatingT"),"T",as.character(WholeTissueSC@meta.data$SubTypes))
#WholeTissueSC@meta.data$SubTypes <- ifelse(WholeTissueSC@meta.data$DefineTypes %in% "Plasma","Plasma",as.character(WholeTissueSC@meta.data$SubTypes))

WholeTissueSC@meta.data$SubTypes <- factor(WholeTissueSC@meta.data$SubTypes,levels=c("Epithelial","MSC","Endo","Glia","B","PlasmaB","T","Myeloid","Mast"))

WholeTissueSC@meta.data$MainTypes <- factor(WholeTissueSC@meta.data$MainTypes,levels=c("Epithelial","Stroma","B","T","Myeloid"))
WholeTissueSC@meta.data$TopTypes <- factor(WholeTissueSC@meta.data$TopTypes,levels=c("Epithelial","Stroma","Immune"))
Patients <- data.frame(Pat_T = c("p1_N","p2_N", "p3_N", "p4_N", "p5_N", "p1_T",  "p2_T",  "p3_T", "p4_T", "p5_T"),
                       Pat_Tissues = c("p7_N","p9_N", "p10_N", "p11_N", "p12_N", "p7_T",  "p9_T",  "p10_T", "p11_T", "p12_T"))


WholeTissueSC <- NormalizeData(WholeTissueSC)
WholeTissueSC <- FindVariableFeatures(object =WholeTissueSC, mean.function = ExpMean, dispersion.function = LogVMR)
WholeTissueSC <- ScaleData(object = WholeTissueSC)
WholeTissueSC <- RunPCA(object = WholeTissueSC,  npcs = 30, verbose = FALSE)
Idents(WholeTissueSC) <- WholeTissueSC$Pat_Tissues

WholeTissueSC <- WholeTissueSC %>% 
  RunHarmony("Pat_Tissues", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(WholeTissueSC, 'harmony')
WholeTissueSC <- WholeTissueSC %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
WholeTissueSC <- WholeTissueSC %>% 
  RunTSNE(reduction = "harmony", dims = 1:10) 
readr::write_rds(WholeTissueSC,path="~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz",compress = "gz")
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")

WholeTissueSC@meta.data$SubTypes <- factor(WholeTissueSC@meta.data$SubTypes,levels=SubTypesClusters$SubTypes)

tiff("WholeTissue_MainTypesDefineUpdateB.tiff",width = 1000,height = 1000,res=280)
DimPlot(WholeTissueSC,reduction = "umap",group.by = "SubTypes",cols = SubTypesClusters$Colors,pt.size = 0.1)+NoLegend()
dev.off()

pdf("ColorDefine.pdf",width = 4,height = 10)
ggplot(Clusters,aes(x=1,y=DefineTypes))+
  geom_tile(aes(fill=DefineTypes))+
  scale_y_discrete(limit=rev(Clusters$DefineTypes),labels=rev(Clusters$DefineTypes))+
  scale_fill_manual(limit=Clusters$DefineTypes,values=Clusters$Colors,guide=F)+
  coord_fixed()+theme(legend.background = element_blank(),panel.background = element_blank(),axis.text.x = element_blank(),
                      axis.title =element_blank(),axis.ticks = element_blank() )
dev.off() 
pdf("ColorDefinePoint.pdf",width = 4,height = 10)
ggplot(Clusters,aes(x=1,y=DefineTypes))+
  geom_point(aes(color=DefineTypes))+
  scale_y_discrete(limit=rev(Clusters$DefineTypes),labels=rev(Clusters$DefineTypes))+
  scale_color_manual(limit=Clusters$DefineTypes,values=Clusters$Colors,guide=F)+
  coord_fixed()+theme(legend.background = element_blank(),panel.background = element_blank(),axis.text.x = element_blank(),
                      axis.title =element_blank(),axis.ticks = element_blank() )
dev.off()

pdf("SubTypesColorDefine1.pdf",width = 4,height = 4)
ggplot(SubTypesClusters,aes(x=1,y=SubTypes))+
  geom_point(aes(color=SubTypes))+
  scale_y_discrete(limit=rev(SubTypesClusters$SubTypes),labels=rev(SubTypesClusters$SubTypes))+
  scale_color_manual(limit=SubTypesClusters$SubTypes,values=SubTypesClusters$Colors,guide=F)+
  coord_fixed()+theme(legend.background = element_blank(),panel.background = element_blank(),axis.text.x = element_blank(),
                      axis.title =element_blank(),axis.ticks = element_blank() )
dev.off() 

####





Idents(WholeTissueSC) <- WholeTissueSC$DefineTypes
WholeTissueSC_ClusterMarker <- FindAllMarkers(object = WholeTissueSC, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)

openxlsx::write.xlsx(WholeTissueSC_ClusterMarker,file="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_SC_ClusterMarkerUpdateB.xlsx")

WholeTissueSC_ClusterMarker <- openxlsx::read.xlsx("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_SC_ClusterMarkerUpdateB.xlsx")
WholeTissueSC_ClusterMarkerF <- WholeTissueSC_ClusterMarker[WholeTissueSC_ClusterMarker$avg_logFC > 1,]
Idents(WholeTissueSC) <- WholeTissueSC$DefineTypes
WholeTissueSC_ClusterMarkerMast <- FindAllMarkers(object = WholeTissueSC, min.pct=0.1,logfc.threshold = 0.25,pseudocount.use = 0.1,only.pos = T,test.use = "MAST")

openxlsx::write.xlsx(WholeTissueSC_ClusterMarkerMast,file="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_SC_ClusterMarkerTestByMAST.xlsx")



#####Plot Figures


WholeTissueListNew <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypeUpdataB.rds.gz")

readr::write_rds(WholeTissueListNew$Myeloid,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/CRC_Myeloid.rds.gz",compress = "gz")

p1 <- FeaturePlot(WholeTissueListA$MSC,feature="PLAU",split.by = "Tissues")
p2 <- FeaturePlot(WholeTissueSC,feature="ICAM1",split.by = "Tissues")
cowplot::plot_grid(plotlist = list(p1,p2))

