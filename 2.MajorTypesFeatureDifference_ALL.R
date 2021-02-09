library(magrittr)
library(ggplot2)
library(nlme)
library(Seurat)
options(stringsAsFactors = F)
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/MarjorTypes/")

###
###

WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineSubtype_AllFunctionTF.rds.gz")
###HallMark
Features <- c("RNA","Metab","BP","React", "Hallmark")
AllFeatureDiff <- lapply(Features,function(AssayName){
  DefaultAssay(WholeTissueSC) <- AssayName
  Idents(WholeTissueSC) <- "SubTypes"
  if(AssayName=="RNA"){
    DEfeatures <- FindAllMarkers(WholeTissueSC,min.pct=0.1,logfc.threshold = 0.25,pseudocount.use = 0.1,only.pos = T)
  }else{
    DEfeatures <- FindAllMarkers(WholeTissueSC,min.pct=0.1,logfc.threshold = 0.25,pseudocount.use = F,only.pos = T)
  }
  return(DEfeatures)
})
names(AllFeatureDiff ) <- Features


readr::write_rds(AllFeatureDiff,
                 path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/MarjorTypes/FeaturesDifferenceAmongSubtypes.rds.gz",compress = "gz")

SubTypesFeaturesDifference <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/MarjorTypes/SubTypesFeaturesDifference.rds.gz")

AllFeatureDiff <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/5.FunctionEnrichment/MarjorTypes/FeaturesDifferenceAmongSubtypes.rds.gz")

openxlsx::write.xlsx(AllFeatureDiff$RNA[AllFeatureDiff$RNA$avg_logFC > 1,],
                     file="DEgenesAmongSubtypes.xlsx")



