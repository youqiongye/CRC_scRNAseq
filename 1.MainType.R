WholeTissueDataN <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Process/WholeTissueData_Filter_nFeatureMito.rds.gz")

EpithelialCells <-  c("EPCAM","KRT8","KRT18")
StromalCells <-  c("COL1A1","COL1A2","COL6A1","COL6A2","VWF","PLVAP","CDH5","S100B")
ImmuneCells <- c("CD52","CD2","CD3D","CD3G","CD3E","CD79A","CD79B","CD14","FCGR3A","CD68","CD83","CSF1R","FCER1G")
MyeloidCells <- c("CD68", "XCR1", "CLEC9A", "CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM")
TCells <- c("NKG7", "KLRC1", "CCR7", "FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D")
BCells<- c("MZB1", "IGHA1", "SELL", "CD19", "AICDA")

Samples <- c("p7_N","p9_N", "p10_N", "p11_N", "p12_N", "p7_T",  "p9_T",  "p10_T", "p11_T", "p12_T")
WholeTissueDataN <- NormalizeData(WholeTissueDataN)
WholeTissueDataN <- FindVariableFeatures(object =WholeTissueDataN, mean.function = ExpMean, dispersion.function = LogVMR)
WholeTissueDataN <- ScaleData(object = WholeTissueDataN)
WholeTissueDataN <- RunPCA(object = WholeTissueDataN,  npcs = 30, verbose = FALSE)
Idents(WholeTissueDataN) <- WholeTissueDataN$Pat_Tissues


ColonData <- WholeTissueDataN %>% 
  RunHarmony("Pat_Tissues", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(ColonData, 'harmony')
ColonData <- ColonData %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
ColonData <- ColonData %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) 


for(i in seq(0.1,0.8,by=0.1)){
  ColonData <- FindClusters(ColonData, resolution =i, algorithm = 1)
}
ColonData@meta.data$MyeloidCells<- apply(ColonData@assays$RNA@data[rownames(ColonData@assays$RNA@data) %in% MyeloidCells,],2,mean)
ColonData@meta.data$TCells <- apply(ColonData@assays$RNA@data[rownames(ColonData@assays$RNA@data) %in% TCells,],2,mean)
ColonData@meta.data$BCells <- apply(ColonData@assays$RNA@data[rownames(ColonData@assays$RNA@data) %in% BCells,],2,mean)
ColonData@meta.data$EpithelialCells <- apply(ColonData@assays$RNA@data[rownames(ColonData@assays$RNA@data) %in% EpithelialCells,],2,mean)
ColonData@meta.data$StromalCells<- apply(ColonData@assays$RNA@data[rownames(ColonData@assays$RNA@data) %in% StromalCells,],2,mean)
ColonData@meta.data$ImmuneCellsScore <- apply(ColonData@assays$RNA@data[rownames(ColonData@assays$RNA@data) %in% ImmuneCells,],2,mean)
ColonData@meta.data$ImmuneCells <- apply(ColonData@meta.data[,c("ImmuneCellsScore","MyeloidCells","TCells","BCells")],1,max)
readr::write_rds(ColonData,path="/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColonData_Harmony.rds.gz",compress = "gz")
Types <- c("MyeloidCells","TCells", "BCells","ImmuneCellsScore","EpithelialCells", "StromalCells")
ThreeTypes <- c("ImmuneCells","EpithelialCells", "StromalCells")
ImmTypes <-  c("MyeloidCells","TCells", "BCells")


p1 <- DimPlot( ColonData, reduction = "umap",label=T,group.by = "RNA_snn_res.0.8",
               cols = colorRampPalette(.cluster_cols,space="rgb")(length(unique(ColonData$RNA_snn_res.0.8))))+NoLegend()
p <-FeaturePlot(ColonData,reduction = "umap",features =  Types,cols = c("gray","orange","orange","red","red"),combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

# pdf(paste("Fig/",Sample,"Cluster_Score.pdf",sep=""),width = 8,height =8)
# cowplot::plot_grid(plotlist = p)
# dev.off()
# #cowplot::plot_grid(plotlist = p)
# pdf(paste("Fig/",Sample,"Cluster_VlnPlot.pdf",sep=""),width = 10,height =3)
ColonDataS$RNA_snn_res.0.8 <- factor(ColonDataS$RNA_snn_res.0.8,levels = 0:(length(unique(ColonDataS$RNA_snn_res.0.8))-1))
Vln_p <- VlnPlot(ColonData, features= Types,group.by = "RNA_snn_res.0.8",
                 pt.size = 0,cols=.cluster_cols[1:length(unique(ColonData@meta.data$RNA_snn_res.0.8))],combine = F)
# dev.off()
ClusterSum <- sapply(split(ColonData@meta.data[,ThreeTypes],ColonData@meta.data$RNA_snn_res.0.8),function(x){apply(x,2,mean)}) %>% t %>% data.frame %>%
  dplyr::mutate(Cluster=rownames(.))
ClusterSum[,1:3] <- apply(ClusterSum[,1:3],2,function(x){(x)})
ClusterSum$ClusterBelong <- apply(ClusterSum[,1:3],1,function(x)match(max(x),x))
ClusterSum$ClusterDefine <- ifelse(ClusterSum$ClusterBelong==1,"ImmuneCells",ifelse(ClusterSum$ClusterBelong==2,"EpithelialCells","StromalCells"))
ClusterSum <- ClusterSum[order(ClusterSum$ClusterDefine),]
ClusterSum.m <- ClusterSum %>% tidyr::gather(key=ThreeTypes,value=ClusterScore,-c(Cluster,ClusterBelong,ClusterDefine))
ClusterDefineP  <- ggplot(ClusterSum.m,aes(y=Cluster,x=ThreeTypes))+
  geom_tile(aes(fill=ClusterScore))+
  scale_x_discrete(limit=ThreeTypes)+
  scale_y_discrete(limit=ClusterSum$Cluster)+
  #coord_fixed(ratio = 2)+
  scale_fill_gradientn(colors = colorRampPalette(c("gray","orange","red"),space = "rgb")(20),guide = F)+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=8,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.ticks=element_blank())
##Immune types
ImmClusterSum <- sapply(split(ColonData@meta.data[,ImmTypes],ColonData@meta.data$RNA_snn_res.0.8),function(x){apply(x,2,mean)}) %>% t %>% data.frame %>%
  dplyr::mutate(Cluster=rownames(.))
ImmClusterSum[,1:3] <- apply(ImmClusterSum[,1:3],2,function(x){(x)})
ImmClusterSum$ClusterBelong <- apply(ImmClusterSum[,1:3],1,function(x)match(max(x),x))
ImmClusterSum$ClusterDefine <- ClusterSum$ClusterDefine[match(ImmClusterSum$Cluster,ClusterSum$Cluster)]
ImmClusterSum$ClusterDefine <- ifelse(ImmClusterSum$ClusterDefine %in% c("EpithelialCells","StromalCells"),
                                      ImmClusterSum$ClusterDefine,
                                      ifelse(ImmClusterSum$ClusterBelong==1,"MyeloidCells",ifelse(ImmClusterSum$ClusterBelong==2,"TCells","BCells")))
ImmClusterSumM <- ImmClusterSum
ImmClusterSum <- ImmClusterSum[order(ImmClusterSum$ClusterDefine),]
ImmClusterSum <- ImmClusterSum[ImmClusterSum$ClusterDefine %in% ImmTypes,]
ImmClusterSum.m <- ImmClusterSum %>% tidyr::gather(key=ImmTypes,value=ClusterScore,-c(Cluster,ClusterBelong,ClusterDefine))
ImmTypesP  <- ggplot(ImmClusterSum.m,aes(y=Cluster,x=ImmTypes))+
  geom_tile(aes(fill=ClusterScore))+
  scale_x_discrete(limit=ImmTypes)+
  scale_y_discrete(limit=ImmClusterSum$Cluster)+
  #coord_fixed(ratio = 2)+
  scale_fill_gradientn(colors = colorRampPalette(c("gray","orange","red"),space = "rgb")(20),guide = F)+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=8,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.ticks=element_blank())
Pre_l <- length(p)
for(i in 1:length(Vln_p)){
  p[[(i+Pre_l)]] <- Vln_p[[i]] + NoLegend()+theme(axis.text.x = element_text(size=7))
}
p[[(length(p)+1)]] <- p1+ NoLegend() + NoAxes()
p[[(length(p)+1)]] <- ClusterDefineP
p[[(length(p)+1)]] <- ImmTypesP


###SingleR

  ref.se <- readr::read_rds("/home/yye/Project/Collaboration/2020ShenLei/Matrix/ProcessData/SingeR_ref.se.rds.gz")
  bp.se <- ref.se$BP
  hpca.se <- ref.se$HPCA
  ColonDataS <- ColonData
  ###PBMC
  ColonDataS_SingleR_Main <- SingleR::SingleR(GetAssayData(ColonDataS@assays$RNA, slot = "data"),
                                                  ref=ref.se,
                                                  method = "cluster",
                                                  clusters =ColonDataS@meta.data$RNA_snn_res.0.8,
                                                  labels = list(bp.se$label.main, hpca.se$label.main))
  SingleR_Result <- data.frame(HPCA =ColonDataS_SingleR_Main@listData$orig.results$HPCA$first.labels,
                               BP =ColonDataS_SingleR_Main@listData$orig.results$BP$first.labels,
                               clusters=sort(unique(ColonDataS@meta.data$RNA_snn_res.0.8)))
  

ColonDataS@meta.data$HPCA <- SingleR_Result$HPCA[match(ColonDataS@meta.data$RNA_snn_res.0.8,SingleR_Result$clusters)]
ColonDataS@meta.data$BP <- SingleR_Result$BP[match(ColonDataS@meta.data$RNA_snn_res.0.8,SingleR_Result$clusters)]

HPCA_p <- DimPlot( ColonDataS, reduction = "umap",label=T,group.by = "HPCA",cols = .cluster_cols)+NoLegend()
BP_p <- DimPlot( ColonDataS, reduction = "umap",label=T,group.by = "BP",cols = .cluster_cols)+NoLegend()
p[[(length(p)+1)]] <- HPCA_p
p[[(length(p)+1)]] <- BP_p
#####


ColonData@meta.data$MainTypes <- ClusterSum$ClusterDefine[match(ColonData@meta.data$RNA_snn_res.0.8,ClusterSum$Cluster)]
ColonData@meta.data$AllTypes <- ImmClusterSumM$ClusterDefine[match(ColonData@meta.data$RNA_snn_res.0.8,ImmClusterSumM$Cluster)]

MainType_p <- DimPlot( ColonData, reduction = "umap",label=T,group.by = "MainTypes",cols = .cluster_cols)+NoLegend()
AllType_p <- DimPlot( ColonData, reduction = "umap",label=T,group.by = "AllTypes",cols = .cluster_cols[c(9,1,12,3,8)])+NoLegend()

p[[(length(p)+1)]] <- MainType_p
p[[(length(p)+1)]] <- AllType_p

png("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/res0.8_Cluster.png",width = 4200,height =2100,res=200)
ggp <-  cowplot::plot_grid(plotlist = p,ncol = 6)
print(ggp)
dev.off()

readr::write_rds(ColonData,path="~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColonData_MaintypeDefine.rds.gz",compress = "gz")


ColonData$AllTypes <- ifelse(ColonData$RNA_snn_res.0.8 %in% 13,"MyeloidCells",ColonData$AllTypes)
ColonData$AllTypes <- ifelse(ColonData$RNA_snn_res.0.8 %in% c(19,24),"TCells",ColonData$AllTypes)

for(i in unique(ColonData$AllTypes)){
  sub <- subset(ColonData,cells=rownames(ColonData@meta.data[ColonData@meta.data$AllTypes %in% i,]))
  readr::write_rds(sub,path=paste("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/",i,"_harmony.rds.gz",sep=""),compress = "gz")
}

####Cluster distribution

EachSampleClusterDis <- lapply(unique(ColonDataS@meta.data$Pat_Tissues),function(Sample){
  seur <- subset(ColonDataS,cells=rownames(ColonDataS@meta.data[ColonDataS@meta.data$Pat_Tissues %in% Sample,]))
  seur@meta.data$AllTypes  %>% table %>% data.frame %>% set_colnames(c("CellTypes","Number")) %>% dplyr::mutate(Per=100*Number/sum(Number))
})

names(EachSampleClusterDis) <- unique(ColonDataS@meta.data$Pat_Tissues)
EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample=rep(names(EachSampleClusterDis),times=unlist(lapply(EachSampleClusterDis,nrow))))
pdf(paste("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/Cluster_Distribution.pdf",sep=""),width = 6,height = 4)
ggplot(data = EachSampleClusterDisA, aes(x = Sample, fill=CellTypes,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+
  scale_x_discrete(limits=Samples)+
  scale_fill_manual(values=.cluster_cols[c(9,1,12,3,8)],name="")+
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black"),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))
dev.off()

EachSampleClusterDisA[,c("Patients","Tissue")] <- data.frame(do.call(rbind,strsplit(as.character(EachSampleClusterDisA$Sample),"_")))

pdf("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/PercentageDefinedDistribution_Line.pdf",
    width = 7,height = 3)
ggplot(data = EachSampleClusterDisA, aes(x = Tissue, color=Patients,y=Per,group=Patients)) + #log2(as.numeric(CopyNumber)))
  geom_point(size=2)+
  facet_wrap(~CellTypes,nrow=1,scale="free")+
  scale_x_discrete(limit=c("N","T"),labels=c("N","T"))+
  scale_y_continuous(expand=c(0.2,0))+
  scale_color_manual(limit=c("p7","p9","p10","p11","p12"),values=RColorBrewer::brewer.pal(5,"Set2"),name="")+
  geom_line(linetype = "dashed")+#stat_compare_means(label = "p.format")+
  ylab("Percentage (%)")+
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black"),
        strip.background = element_blank(),strip.text = element_text(size=12),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))+
  ggpubr::stat_compare_means(comparisons = list(v1=c("N","T")),label="p.format",method = "t.test",paired =T )
dev.off() 


pdf("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/Cluster_Distribution_Number.pdf",width = 8,height = 4)
ggplot(data = EachSampleClusterDisA, aes(x = Sample, y= CellTypes)) + #log2(as.numeric(CopyNumber)))
  geom_tile(aes(fill=CellTypes),color="black")+
  geom_text(aes(label=Number))+
  coord_fixed(ratio = 0.8)+
  scale_y_discrete(limits = c("MyeloidCells","TCells", "BCells","EpithelialCells", "StromalCells"))+labs(y="percentage(%)")+
  scale_x_discrete(limits=names(EachSampleClusterDis))+
  scale_fill_manual(values=.cluster_cols[c(9,1,12,3,8)],name="",guide=F)+
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        axis.text.y=element_text(color = "black"),axis.ticks = element_blank(),
        axis.title = element_blank(),legend.text = element_text(size=10))
dev.off()
###


####
#plasma cells CD38+; SDC1:CD138+
FeaturePlot(ColonData,reduction = "umap",features = c("CD19"),cols = c("gray","orange","red"))

