setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment")
library(fgsea)
library(magrittr)
library(ggplot2)
library(nlme)
library(Seurat)
source("/home/yye/Data/Public/ToolsData/ulcerative_colitis/scores.r")
###KEGG Pathway loading
MetabolismPathways <- gmtPathways("/home/yye/Data/Public/Genelist/KEGG_metabolism.gmt")
NonMetabolismPathways <- gmtPathways("/home/yye/Data/Public/Genelist/nonMetabolic_KEGG.gmt")
KeggPathways <- c(MetabolismPathways,NonMetabolismPathways)
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")

scores <- score_cells(seur=WholeTissueSC, names=KeggPathways, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
readr::write_rds(scores,path="KEGG_PathwayScoreBycssmillie.rds")
###TORC1 signaling

###Reactome pathways
ReactomePathways <- gmtPathways("/home/yye/Data/Public/Genelist/ReactomePathways.gmt")
ReactomePathwaysscores <- score_cells(seur=WholeTissueSC, names=ReactomePathways, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
readr::write_rds(ReactomePathwaysscores ,path="ReactomePathwaysScoreBycssmillie.rds")
###Biological pathway
BPPathways <- gmtPathways("/home/yye/Data/Public/Genelist/c5.bp.v7.1.symbols.gmt")
BPPathwaysscores <- score_cells(seur=WholeTissueSC, names=BPPathways, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
readr::write_rds(BPPathwaysscores,path="RBPPathwaysScoreBycssmillie.rds.gz",compress = "gz")
###hallmark
HMPathways <- gmtPathways("/home/yye/Data/Public/Genelist/h.all.v6.1.symbols.gmt")
HMPathwaysscores <- score_cells(seur=WholeTissueSC, names=HMPathways, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
readr::write_rds(HMPathwaysscores,path="HMPathwaysScoreBycssmillie.rds.gz",compress = "gz")

###
HMPathwaysscores <- readr::read_rds("HMPathwaysScoreBycssmillie.rds.gz")
####BP
HMPathwaysscoresMatrix <- as.matrix(HMPathwaysscores)
HMPathwaysscoresMatrix  <- t(HMPathwaysscoresMatrix )
HMPathwaysscoresMatrix <- apply(HMPathwaysscoresMatrix,2,function(x)signif(x,digits = 3))
#rownames(HMPathwaysscoresMatrix ) <- names(BPPathways)
colnames(HMPathwaysscoresMatrix) <- rownames(WholeTissueSC@meta.data)
HMPathwaysscoresMatrixSeurat <- CreateSeuratObject(counts = HMPathwaysscoresMatrix )
HMPathwaysscoresMatrixSeurat@meta.data <- WholeTissueSC@meta.data

HMPathwaysscoresMatrix1 <- HMPathwaysscoresMatrix
#rownames(HMPathwaysscoresMatrix1 ) <- names(BPPathways)
colnames(HMPathwaysscoresMatrix1) <- gsub("Epithelial_|MSC_|Endo_|Glia_|B_|^T_|ProliferatingT_|ILCs_|Myeloid_|Mast_|Plasma","",rownames(WholeTissueSC@meta.data))
HMPathwaysscoresMatrix1Seurat <- CreateSeuratObject(counts = HMPathwaysscoresMatrix1 )


####Merge pathways into different subtypes of CRC as new RNA expression assay

####KEGG pathways
scores <- readr::read_rds("KEGG_PathwayScoreBycssmillie.rds")
scoresMatrix <- dplyr::bind_rows(scores)
scoresMatrix  <- t(scoresMatrix )
rownames(scoresMatrix ) <- c(names(MetabolismPathways),names(NonMetabolismPathways))
colnames(scoresMatrix) <- rownames(WholeTissueSC@meta.data)
scoresMatrixSeurat <- CreateSeuratObject(counts = scoresMatrix )
scoresMatrixSeurat@meta.data <- WholeTissueSC@meta.data
scoresMatrix1 <- scoresMatrix
rownames(scoresMatrix1 ) <- c(names(MetabolismPathways),names(NonMetabolismPathways))
colnames(scoresMatrix1) <- gsub("Epithelial_|MSC_|Endo_|Glia_|B_|^T_|ProliferatingT_|ILCs_|Myeloid_|Mast_|Plasma_","",rownames(WholeTissueSC@meta.data))
scoresMatrix1Seurat <- CreateSeuratObject(counts = scoresMatrix1 )


BPPathwaysscores <- readr::read_rds("Output/RBPPathwaysScoreBycssmillie.rds.gz")
####BP
BPPathwaysscoresMatrix <- as.matrix(BPPathwaysscores)
BPPathwaysscoresMatrix  <- t(BPPathwaysscoresMatrix )
BPPathwaysscoresMatrix <- apply(BPPathwaysscoresMatrix,2,function(x)signif(x,digits = 3))
#rownames(BPPathwaysscoresMatrix ) <- names(BPPathways)
colnames(BPPathwaysscoresMatrix) <- rownames(WholeTissueSC@meta.data)
BPPathwaysscoresMatrixSeurat <- CreateSeuratObject(counts = BPPathwaysscoresMatrix )
BPPathwaysscoresMatrixSeurat@meta.data <- WholeTissueSC@meta.data

BPPathwaysscoresMatrix1 <- BPPathwaysscoresMatrix
#rownames(BPPathwaysscoresMatrix1 ) <- names(BPPathways)
colnames(BPPathwaysscoresMatrix1) <- gsub("Epithelial_|MSC_|Endo_|Glia_|B_|^T_|ProliferatingT_|ILCs_|Myeloid_|Mast_","",rownames(WholeTissueSC@meta.data))
BPPathwaysscoresMatrix1Seurat <- CreateSeuratObject(counts = BPPathwaysscoresMatrix1 )
####Reactome
ReactomePathwaysscores <- readr::read_rds("Output/ReactomePathwaysScoreBycssmillie.rds")
ReactomePathwaysscoresMatrix <- as.matrix(ReactomePathwaysscores)
ReactomePathwaysscoresMatrix  <- t(ReactomePathwaysscoresMatrix )
ReactomePathwaysscoresMatrix <- apply(ReactomePathwaysscoresMatrix,2,function(x)signif(x,digits = 3))
#rownames(ReactomePathwaysscoresMatrix ) <- names(ReactomePathways)
colnames(ReactomePathwaysscoresMatrix) <- rownames(WholeTissueSC@meta.data)
ReactomePathwaysscoresMatrixSeurat <- CreateSeuratObject(counts = ReactomePathwaysscoresMatrix )
ReactomePathwaysscoresMatrixSeurat@meta.data <- WholeTissueSC@meta.data
ReactomePathwaysscoresMatrix1 <- ReactomePathwaysscoresMatrix
#rownames(ReactomePathwaysscoresMatrix1 ) <-names(ReactomePathways)
colnames(ReactomePathwaysscoresMatrix1) <- gsub("Epithelial_|MSC_|Endo_|Glia_|B_|^T_|ProliferatingT_|ILCs_|Myeloid_|Mast_","",rownames(WholeTissueSC@meta.data))
ReactomePathwaysscoresMatrix1Seurat <- CreateSeuratObject(counts = ReactomePathwaysscoresMatrix1 )

WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypes.rds.gz")

WholeTissueList1 <- lapply(names(WholeTissueList),function(x){
  sub <- WholeTissueList[[x]]
  subScore <- subset(scoresMatrix1Seurat,cells=rownames(sub@meta.data))
  subBPPathwaysscores <- subset(BPPathwaysscoresMatrix1Seurat,cells=rownames(sub@meta.data))
  subReactomePathwaysscores <- subset(ReactomePathwaysscoresMatrix1Seurat,cells=rownames(sub@meta.data))
  subHMPathwaysscores <- subset(HMPathwaysscoresMatrix1Seurat,cells=rownames(sub@meta.data))
  
  sub@assays$Metab <- subScore@assays$RNA
  sub@assays$GO <- subBPPathwaysscores@assays$RNA
  sub@assays$React <- subReactomePathwaysscores@assays$RNA
  sub@assays$Hallmark <- subHMPathwaysscores@assays$RNA
  return(sub)
})
names(WholeTissueList1) <- names(WholeTissueList)
###Add hallmarks GO


###
#####
WholeTissueSC@assays$PathwayScoreUC <- scoresMatrixSeurat@assays$RNA


readr::write_rds(WholeTissueListA,path="~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypes_AddScore.rds.gz")

FeaturePlot(WholeTissueList$B,feature=c("IGHG1","IGHM","IGHD","IGHA1"))

WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")
WholeTissueListA <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtype_AllFunctionTF.rds.gz")

FeaturePlot(WholeTissueSC,feature="ATF3")
FeaturePlot(WholeTissueListA$MSC,feature="ATF3",col=c("gray93","red"))

VlnPlot(WholeTissueListA$Myeloid,feature="ATF3",group.by = "DefineTypes")
VlnPlot(WholeTissueListA$Mast,feature="CD274",group.by = "DefineTypes")


MSC <- WholeTissueListA$MSC
DefaultAssay(MSC) <- "TF"
FeaturePlot(WholeTissueListA$MSC,feature="FABP4",col=c("gray93","blue"),split.by = "Tissues") 
FeaturePlot(WholeTissueListA$MSC,feature="GO-REGULATION-OF-CELLULAR-RESPONSE-TO-HYPOXIA",col=c("gray93","lightskyblue","blue"),split.by = "Tissues")
FeaturePlot(WholeTissueListA$MSC,feature="RARRES2",col=c("gray93","orange","red"),split.by = "Tissues")

FeaturePlot(WholeTissueSC,feature="RARRES2",col=c("gray93","blue"),split.by = "Tissues")

VlnPlot(WholeTissueListA$MSC,feature="RARRES2",split.by = "Tissues",group.by = "DefineTypes")
VlnPlot(WholeTissueListA$Myeloid,feature="STAT1...",group.by = "DefineTypes")

####DE features
names(WholeTissueListA)


WholeTissueListA <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtype_AllFunctionTF.rds.gz")

WholeTissueFeatureDiff <- lapply(c("Epithelial","MSC","Endo","B","PlasmaB","T","Myeloid"),function(x){
 AA <-  lapply(c("RNA","Metab","GO","React","TF"),function(AssayName){
    DefaultAssay(x) <- AssayName
    Idents(x) <- "DefineTypes"
    DEfeatures <- FindAllMarkers(x,pct=0.1,logfc.threshold = 0.25,pseudocount.use = 0.1,only.pos = T)
    print(AssayName)
    return(DEfeatures)
  })
 names(AA) <- c("RNA","Metab","GO","React","TF")
 return(AA)
})
readr::write_rds(WholeTissueFeatureDiff,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/WholeTissueFeatureDiff_EachSubTypes.rds.gz",compress = "gz")


