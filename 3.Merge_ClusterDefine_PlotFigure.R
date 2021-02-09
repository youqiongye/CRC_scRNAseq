library(Seurat)
library(ggplot2)
library(magrittr)
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue")
PatientColor <- data.frame(Patients = c("p1","p2","p3","p4","p5"),PatColor=c("#D53E4F", "purple2", "#FDAE61", "#66C2A5", "#3288BD"))
##LRRC15
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")

WholeTissueSC_ClusterMarker <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_SC_ClusterMarkerUpdateB.xlsx")
WholeTissueSC_ClusterMarkerFC0.5 <- WholeTissueSC_ClusterMarker[WholeTissueSC_ClusterMarker$avg_logFC > 0.5,]
WholeTissueSC_ClusterMarkerFC0.5$DefineTypes <- WholeTissueSC_ClusterMarkerFC0.5$cluster

Clusters <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
SubClusters <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/SubTypesColorCode.xlsx")
SubClusters[,"SubTypes"] <- ifelse(SubClusters$SubTypes %in% "Plasma","PlasmaB",SubClusters$SubTypes)
DefineMarkers <- openxlsx::read.xlsx("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/AllCellTypeMarkers-0513.xlsx",sheet="DefineTypesMarkers")
PanMarkers <-  openxlsx::read.xlsx("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/AllCellTypeMarkers-0513.xlsx",sheet="PanMarkers")
AddGenes <- setdiff(c(DefineMarkers$GeneSymbol,PanMarkers$GeneSymbol),WholeTissueSC_ClusterMarkerFC0.5$gene)
AddGenes <- DefineMarkers[DefineMarkers$GeneSymbol %in% AddGenes,] %>% dplyr::select(c(DefineTypes,GeneSymbol)) %>% set_colnames(c("DefineTypes","gene")) %>%
  dplyr::mutate(avg_logFC=rep(1,times=nrow(.)))
DefineMarkers$DefineTypes <- factor(DefineMarkers$DefineTypes,levels=Clusters$DefineTypes)
PlotGenes <- unique(c(DefineMarkers$GeneSymbol,PanMarkers$GeneSymbol))

PlotGenesLeft <- unique(c(DefineMarkers[DefineMarkers$MainTypes %in% c("Epithelial","Stroma"),]$GeneSymbol,PanMarkers[PanMarkers$TopTypes %in% c("Epithelial","Stroma"),]$GeneSymbol))
PlotGenesRight <- unique(c(DefineMarkers[!(DefineMarkers$MainTypes %in% c("Epithelial","Stroma")),]$GeneSymbol,
                           PanMarkers[!(PanMarkers$TopTypes %in% c("Epithelial","Stroma")),]$GeneSymbol))


WholeTissueSC_ClusterMarkerFC0.5M <- WholeTissueSC_ClusterMarkerFC0.5 %>% dplyr::select(c(DefineTypes,gene,avg_logFC)) %>%  
  dplyr::bind_rows(AddGenes) %>% tidyr::spread(DefineTypes,avg_logFC,fill=0)
WholeTissueSC_ClusterMarkerFC0.5M <- WholeTissueSC_ClusterMarkerFC0.5M[,c("gene",Clusters$DefineTypes)]
#WholeTissueSC_ClusterMarkerFC0.5M <- WholeTissueSC_ClusterMarkerFC0.5M[!is.na(WholeTissueSC_ClusterMarkerFC0.5M$gene),Clusters$DefineTypes]
WholeTissueSC_ClusterMarkerFC0.5M <- WholeTissueSC_ClusterMarkerFC0.5M[do.call(order, c(WholeTissueSC_ClusterMarkerFC0.5M[,2:ncol(WholeTissueSC_ClusterMarkerFC0.5M)], list(decreasing=TRUE))),]


WholeTissue_ClusterID <- lapply(split(WholeTissueSC@meta.data,list(WholeTissueSC$DefineTypes)), function(x)rownames(x))

WholeTissueSCm <- WholeTissueSC@assays$RNA@data%>% as.matrix
colnames(WholeTissueSCm) <- gsub("\\.","-",colnames(WholeTissueSCm))

WholeTissueSCm <- WholeTissueSCm[apply(WholeTissueSCm,1,sum) > 0,]

#WholeTissueSCm <- WholeTissueSCm[rownames(WholetissueCom) %in% VariableFeatures(WholetissueHarmony),]

ClusterMean <- t(apply(WholeTissueSCm,1,function(x){
  lapply(WholeTissue_ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))

readr::write_rds(ClusterMean,path="EachClusterGeneMeanExpression.rds")
ClusterMean <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/EachClusterGeneMeanExpression.rds")
ClusterMean <- ClusterMean[apply(ClusterMean,1,sum)>0,]

ClusterMeanM <- ClusterMean[match(WholeTissueSC_ClusterMarkerFC0.5M$gene,rownames(ClusterMean)),]
ClusterMeanM[,1:ncol(ClusterMeanM)] <- t(apply(ClusterMeanM,1,scale))

Sample_order <- Clusters


ClusterMeanM[ClusterMeanM< -3] <- -3
ClusterMeanM[ClusterMeanM> 3] <- 3


Wholetissue_right_ha = ComplexHeatmap::rowAnnotation(foo =  ComplexHeatmap::anno_mark(at = match(PlotGenes,WholeTissueSC_ClusterMarkerFC0.5M$gene), labels = PlotGenes ))

Wholetissue_RightPanel_ha = ComplexHeatmap::rowAnnotation(foo =  ComplexHeatmap::anno_mark(at = match(PlotGenesRight,WholeTissueSC_ClusterMarkerFC0.5M$gene), labels = PlotGenesRight))

Wholetissue_LeftPanel_ha = ComplexHeatmap::rowAnnotation(foo =  ComplexHeatmap::anno_mark(at = match(PlotGenesLeft,WholeTissueSC_ClusterMarkerFC0.5M$gene), labels = PlotGenesLeft))

Wholetissue_barcolor <- Clusters$Colors
names(Wholetissue_barcolor) <- Clusters$DefineTypes
Wholetissue_SampleBar <-  Sample_order$DefineTypes
names(Wholetissue_SampleBar) <-  Clusters$DefineTypes

Wholetissue_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Wholetissue_SampleBar,
                                                     col = list(bar = Wholetissue_barcolor))


pdf("Final_Wholetissue_HeatmapUpdateB.pdf",width = 8,height = 8)
ComplexHeatmap::Heatmap(ClusterMeanM,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,  name = "mat",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),#colorRampPalette(c("blue","white","red"),space="rgb")(100),#colorRampPalette(c("magenta", "black", "yellow"))(256),#
                        show_row_names = F,show_column_names = F,right_annotation = Wholetissue_right_ha,top_annotation = Wholetissue_ha_bar)
dev.off()



pdf("Final_Wholetissue_HeatmapUpdateB_Twoside.pdf",width = 8,height = 8)
ComplexHeatmap::Heatmap(ClusterMeanM,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,column_order= names(Wholetissue_SampleBar),use_raster = TRUE,  name = "mat",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),#colorRampPalette(c("blue","white","red"),space="rgb")(100),#colorRampPalette(c("magenta", "black", "yellow"))(256),#
                        show_row_names = F,show_column_names = F,
                        right_annotation = Wholetissue_RightPanel_ha,left_annotation =Wholetissue_LeftPanel_ha ,top_annotation = Wholetissue_ha_bar)
dev.off()



pdf("Final_Wholetissue_HeatmapMagentaYellow.pdf",width = 8,height = 8)
ComplexHeatmap::Heatmap(ClusterMeanM,cluster_rows = FALSE, cluster_columns = FALSE, cluster_column_slices=F,use_raster = TRUE,  name = "mat",
                        col=colorRampPalette(c("magenta", "black", "yellow"))(256),#colorRampPalette(c("blue","white","red"),space="rgb")(100),#
                        show_row_names = F,show_column_names = F,right_annotation = Wholetissue_right_ha,top_annotation = Wholetissue_ha_bar)
dev.off()

ClustersColor <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
Tcolor <- ClustersColor[ClustersColor$SubTypes %in% "T",]

p1 <- DimPlot(WholeTissueList$T,reduction = "umap",group.by="DefineTypes",cols =Tcolor$Colors)+
  coord_fixed()

p2 <- FeaturePlot(WholeTissueList$T,feature="GPR171",cols=c("gray93","blue"))+
  coord_fixed()
p3 <- VlnPlot(WholeTissueList$T,feature="GPR171",group.by="DefineTypes",cols =Tcolor$Colors,pt.size = 0)+NoLegend()
pdf("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/GPR171_InCRC_Tcells.pdf",width = 6,height = 10)
cowplot::plot_grid(plotlist = list(p1,p2,p3),ncol=1,rel_heights = c(3.2,3.2,3.6))
dev.off()
SubTypesClusters <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/SubTypesColorCode.xlsx")

PP1 <- FeaturePlot(WholeTissueSC,feature="GPR171",cols=c("gray93","blue"))+
  coord_fixed()
p1 <- DimPlot(WholeTissueSC,reduction = "umap",group.by="SubTypes",cols =SubTypesClusters$Colors)+
  coord_fixed()
p3 <- VlnPlot(WholeTissueSC,feature="GPR171",group.by="SubTypes",cols =SubTypesClusters$Colors,pt.size = 0)+NoLegend()
pdf("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/GPR171_InCRC_Allcells.pdf",width = 5,height = 10)
cowplot::plot_grid(plotlist = list(p1,PP1,p3),ncol=1,rel_heights = c(3.2,3.2,3.6))
dev.off()

####
WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypeUpdataB.rds.gz")
WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")
ClustersColor <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
ClustersColor$MainTypes <- unique(WholeTissueSC@meta.data[,c("DefineTypes","MainTypes")])$MainTypes[match(ClustersColor$DefineTypes,unique(WholeTissueSC@meta.data[,c("DefineTypes","MainTypes")])$DefineTypes)]
#ClustersColor$DefineTypes <- factor(ClustersColor$DefineTypes,levels=names(table(WholeTissueSC$DefineTypes)))
#ClustersColor <- ClustersColor[order(ClustersColor$DefineTypes),]

DistributionPlot <- lapply(names(table(WholeTissueSC$MainTypes)),function(x){
  SubSeurat <- subset(WholeTissueSC,cells=rownames(WholeTissueSC@meta.data[WholeTissueSC$MainTypes %in% x,]))
  MSCcolor <- ClustersColor[ClustersColor$MainTypes %in% x,]
  ###All samples distribution.
  
  SubSeurat@meta.data$DefineTypes <- factor(SubSeurat@meta.data$DefineTypes,levels=names(table(WholeTissueSC@meta.data[WholeTissueSC$MainTypes %in% x,]$DefineTypes)[table(WholeTissueSC@meta.data[WholeTissueSC$MainTypes %in% x,]$DefineTypes)>0]))
  
  SubCellsDisA <- table(SubSeurat@meta.data[,c("Pat_Tissues","DefineTypes")]) %>% 
    data.frame %>% set_colnames(c("Sample","CellTypes","Number"))
  
  SubCellsDisA <- lapply(split(SubCellsDisA,SubCellsDisA$Sample),function(X){
    X%>% dplyr::mutate(Per=100*Number/sum(Number))
  }) %>% dplyr::bind_rows(.)
  
  SubCellsDisA[,c("PatientID","Tissues")] <- data.frame(do.call(rbind,strsplit(as.character(SubCellsDisA$Sample),"_")))
  SubCellsDisA$Tissues <- factor(SubCellsDisA$Tissues,levels=c("N","T"),labels=c("Normal","Tumor"))
  p <- ggplot(data = SubCellsDisA, aes(x = PatientID, fill=CellTypes,y=Per)) + #log2(as.numeric(CopyNumber)))
    geom_bar(stat = "identity",width = 0.6)+
    scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+
    facet_wrap(~Tissues)+
    scale_x_discrete(limits=c(paste("p",1:5,sep="")))+
    scale_fill_manual(limit=MSCcolor$DefineTypes,
                      values= MSCcolor$Colors,name="",guide=F)+
    theme(axis.text.x=element_text(color = "black",hjust=0.5,size=9),
          panel.background = element_blank(),panel.grid=element_blank(),strip.background = element_blank(),strip.text = element_text(size=12),
          legend.title=element_blank(), axis.text.y=element_text(color = "black"),
          axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=8))
  return(p)
  
})



pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/DistributionPlot.pdf",width = 4,height = 20)

cowplot::plot_grid(plotlist = DistributionPlot,ncol=1)
dev.off()
Idents(WholeTissueSC) <- "SubTypes"
DistributionPlotData <- lapply(names(table(WholeTissueSC$SubTypes))[c(1:3,5,7:9)],function(x){
  grep("Plasma",Cells(WholeTissueSC),value=T)
  SubSeurat <- subset(WholeTissueSC,idents= x) #grep("Plasma",colnames(WholeTissueSC),value=T)
  MSCcolor <- ClustersColor[ClustersColor$MainTypes %in% x,]
  ###All samples distribution.
  
  SubSeurat@meta.data$DefineTypes <- factor(SubSeurat@meta.data$DefineTypes,
                                            levels=names(table(WholeTissueSC@meta.data[WholeTissueSC$SubTypes %in% x,]$DefineTypes)[table(WholeTissueSC@meta.data[WholeTissueSC$SubTypes %in% x,]$DefineTypes)>0]))
  
  SubCellsDisA <- table(SubSeurat@meta.data[,c("Pat_Tissues","DefineTypes")]) %>% 
    data.frame %>% set_colnames(c("Sample","CellTypes","Number"))
  
  SubCellsDisA <- lapply(split(SubCellsDisA,SubCellsDisA$Sample),function(X){
    X%>% dplyr::mutate(Per=100*Number/sum(Number))
  }) %>% dplyr::bind_rows(.)
  
  SubCellsDisA[,c("PatientID","Tissues")] <- data.frame(do.call(rbind,strsplit(as.character(SubCellsDisA$Sample),"_")))
  SubCellsDisA$Tissues <- factor(SubCellsDisA$Tissues,levels=c("N","T"),labels=c("Normal","Tumor"))
   return(SubCellsDisA)
})
names(DistributionPlotData) <- names(table(WholeTissueSC$SubTypes))[c(1:3,5,7:9)]

DiffCellTypes <- lapply(DistributionPlotData ,function(sub){
  sub <- sub %>% dplyr::select(-c(Sample,Number) ) %>% tidyr::spread(Tissues,Per)
  sapply(split(sub,sub$CellTypes),function(x){
    aa <- data.frame(TumorMed=mean(x$Tumor),NormalMed=mean(x$Normal),Diff=mean(x$Tumor)-mean(x$Normal),
                     pval=t.test(x$Tumor,x$Normal,paired = T)$p.value)
    return(aa)
  }) %>% t %>% data.frame %>% dplyr::mutate(DefineTypes=rownames(.))
}) %>% dplyr::bind_rows()
DiffCellTypes$Diff <- as.numeric(DiffCellTypes$Diff)
DiffCellTypes <- merge(DiffCellTypes,unique(WholeTissueSC@meta.data[,c("SubTypes","DefineTypes")]),by="DefineTypes")
DiffCellTypes <- DiffCellTypes[order(DiffCellTypes$SubTypes,DiffCellTypes$Diff),]
openxlsx::write.xlsx(DiffCellTypes,file="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_InfiltrationDiff.xlsx")

DiffCellTypesSig <- DiffCellTypes[DiffCellTypes$Diff >0 & DiffCellTypes$pval < 0.05,]

WholeTissueSC_ClusterMarker <- openxlsx::read.xlsx("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_SC_ClusterMarkerUpdateB.xlsx")
TumorUpMarkers <- WholeTissueSC_ClusterMarker[WholeTissueSC_ClusterMarker$avg_logFC > 0.75 & WholeTissueSC_ClusterMarker$cluster %in% DiffCellTypesSig$DefineTypes[grep("Prolifer",DiffCellTypesSig$DefineTypes,invert = T)],]
TumorUpMarkers$gene %>% unique
write.table(unique(TumorUpMarkers$gene),file="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/TumorUpCellTypesMarkers.txt",quote = F,row.names = F,sep="\t")

DiffCellTypes$Class <- ifelse(DiffCellTypes$Diff >0 & DiffCellTypes$pval < 0.05,"Tumor",ifelse(DiffCellTypes$Diff < 0 & DiffCellTypes$pval < 0.05,"Normal","Both"))
DiffCellTypes$Diff <- signif(as.numeric(DiffCellTypes$Diff),digits = 3)
DiffCellTypes$Per <- -log2(as.numeric(DiffCellTypes$TumorMed)/as.numeric(DiffCellTypes$NormalMed))
DiffCellTypes$Class <- factor(DiffCellTypes$Class,levels=c("Normal","Both","Tumor"))
pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/InfiltrationDifference.pdf",width = 4,height = 12)

ggplot(DiffCellTypes,aes(x=1,y=DefineTypes))+
  geom_tile(aes(fill=Diff))+
  scale_fill_gradientn(limit=c(-max(abs(DiffCellTypes$Diff)),max(abs(DiffCellTypes$Diff))),colours = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(8,"RdYlBu")),space="rgb")(100))+
  geom_text(aes(label=Diff))+
  geom_tile(data=DiffCellTypes[DiffCellTypes$Class %in% "Tumor",],aes(x=1,y=DefineTypes),color="red",fill=NA)+
  geom_tile(data=DiffCellTypes[DiffCellTypes$Class %in% "Normal",],aes(x=1,y=DefineTypes),color="blue",fill=NA)+
  geom_tile(data=DiffCellTypes[DiffCellTypes$Class %in% "Both",],aes(x=1,y=DefineTypes),color="gray",fill=NA)+
  scale_y_discrete(limit=DiffCellTypes[order(DiffCellTypes$Class,DiffCellTypes$Diff),]$DefineTypes)+
  theme(axis.title = element_blank(),axis.ticks = element_blank(),panel.background = element_blank(),axis.text.x = element_blank())
dev.off() 

openxlsx::write.xlsx(DiffCellTypes,file="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/DiffCellTypes_Normal_Tumor.xlxs")
###Patients colors
PatientUmapPlot <- lapply(c("Epithelial","MSC","Endo","T","Myeloid","Mast","B","PlasmaB"),function(x){
  p <- DimPlot(WholeTissueList[[x]],reduction = "umap",group.by="orig.ident",cols = as.character(PatientColor$PatColor))+NoLegend()
  return(p)
})

PatientUmapPlotN <- lapply(c("Epithelial","MSC","Endo","T","Myeloid","Mast","B","PlasmaB"),function(x){
  sub <- WholeTissueList[[x]]
  sub <- subset(sub,cells = rownames(sub@meta.data[sub@meta.data$Tissues %in% "N",]))
  p <- DimPlot(sub,reduction = "umap",group.by="orig.ident",cols = as.character(PatientColor$PatColor))+NoLegend()
  return(p)
})
PatientUmapPlotT <- lapply(c("Epithelial","MSC","Endo","T","Myeloid","Mast","B","PlasmaB"),function(x){
  sub <- WholeTissueList[[x]]
  sub <- subset(sub,cells = rownames(sub@meta.data[sub@meta.data$Tissues %in% "T",]))
  p <- DimPlot(sub,reduction = "umap",group.by="orig.ident",cols = as.character(PatientColor$PatColor))+NoLegend()
  return(p)
})


tiff("PatientUmapPlotSubTypes Plot.tiff",width = 2800,height = 1400,res=200)
cowplot::plot_grid(plotlist = PatientUmapPlot,ncol = 4)
dev.off()


tiff("PatientUmapPlotSubTypes PlotNT.tiff",width = 5600,height = 1400,res=200)
cowplot::plot_grid(plotlist = c(PatientUmapPlotT,PatientUmapPlotN),ncol = 8)
dev.off()

###Subtypes distribution
###MainTypes distribution



PatientUmapPlotNT <- lapply(c("MSC","Endo","T","Myeloid","Mast","B","PlasmaB","Epithelial"),function(x){
 
  #sub <- subset(sub,cells = rownames(sub@meta.data[sub@meta.data$Tissues %in% "T",]))
  p <- DimPlot(WholeTissueList[[x]],reduction = "umap",group.by="orig.ident",cols = as.character(PatientColor$PatColor),split.by = "Tissues")+NoLegend()
  return(p)
})
names(PatientUmapPlotNT) <- c("MSC","Endo","T","Myeloid","Mast","B","PlasmaB","Epithelial")
tiff("PatientUmapPlotSubTypes PlotNT_.tiff",width = 1400,height = 4200,res=200)
cowplot::plot_grid(plotlist = PatientUmapPlotNT[c("Endo","T","Mast","B","PlasmaB","Epithelial")],ncol = 1)
dev.off()

tiff("PatientUmapPlotNT_.tiff",width = 800,height = 800,res=200)
DimPlot(WholeTissueSC,reduction = "umap",group.by="Tissues",cols = c("blue2","tomato"),pt.size = 0.1)+NoLegend()
dev.off()

WholeTissueList
WholeTissueSC_ClusterMarker <- openxlsx::read.xlsx("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/WholeTissue_SC_ClusterMarkerUpdateB.xlsx")
WholeTissueSC_ClusterMarkerF <- WholeTissueSC_ClusterMarker[WholeTissueSC_ClusterMarker$avg_logFC > 1,]

###Pan marker
PanMarkers <- c("EPCAM","COL1A1","COL3A1","CD19","MZB1","PTPRC","CD3E","CD14")
p <- FeaturePlot(WholeTissueSC,reduction = "umap",features = PanMarkers,cols = c("gray","orange","red"),pt.size = 0.5,combine = F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+theme(axis.title = element_blank(),axis.text = element_blank(),
                                      axis.ticks = element_blank(),axis.line = element_line(color="gray",linetype = "dashed"))
  
}

tiff("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/PanMarker.tiff",width = 4200,height = 600,res=200)
cowplot::plot_grid(plotlist = p,ncol = 8)
dev.off()
