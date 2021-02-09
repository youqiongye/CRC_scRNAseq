library(corrplot)
library(survival)
library(ggpubr)
library(maxstat)
options(stringsAsFactors = F)
setwd("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputeTCGA/CIBERSORTx/")
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")

WholeTissueSC_ClusterMarker <- openxlsx::read.xlsx("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_SC_ClusterMarkerUpdateB.xlsx")
WholeTissueSC_ClusterMarkerFC0.5 <- WholeTissueSC_ClusterMarker[WholeTissueSC_ClusterMarker$avg_logFC > 0.5,]

NT_SCexpNorm <- WholeTissueSC@assays$RNA@data %>% as.matrix
T_SCexpNorm <- NT_SCexpNorm[,rownames(WholeTissueSC@meta.data[WholeTissueSC@meta.data$Tissues %in% "T" ,])]
T_SCexpNorm <- T_SCexpNorm[apply(T_SCexpNorm,1,sum) >1,]
T_SCexpNorm <- apply(T_SCexpNorm,2,function(x){signif(x,digits = 3)})
T_SCexpNorm <- T_SCexpNorm[rownames(T_SCexpNorm) %in% WholeTissueSC_ClusterMarker$gene,]

colnames(T_SCexpNorm) <- WholeTissueSC@meta.data$SubTypes[match(colnames(T_SCexpNorm),rownames(WholeTissueSC@meta.data))]
write.table(T_SCexpNorm,file="~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputeTCGA/Cibersortx_T_Major_Input.txt",quote = F,row.names = T,col.names = T,sep="\t")


#T_SCexpNorm <- read.delim("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputeTCGA/Cibersortx_Input_2020.txt")
NT_SCexpNorm <- WholeTissueSC@assays$RNA@data %>% as.matrix
N_SCexpNorm <- NT_SCexpNorm[,rownames(WholeTissueSC@meta.data[WholeTissueSC@meta.data$Tissues %in% "N" ,])]
N_SCexpNorm <- N_SCexpNorm[apply(N_SCexpNorm,1,sum) >1,]
N_SCexpNorm <- apply(N_SCexpNorm,2,function(x){signif(x,digits = 3)})
N_SCexpNorm <- N_SCexpNorm[rownames(N_SCexpNorm) %in% WholeTissueSC_ClusterMarker$gene,]
colnames(N_SCexpNorm) <- WholeTissueSC@meta.data$SubTypes[match(colnames(N_SCexpNorm),rownames(WholeTissueSC@meta.data))]
write.table(N_SCexpNorm,file="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputeTCGA/Cibersortx_N_Major_Input.txt",quote = F,row.names = T,col.names = T,sep="\t")
