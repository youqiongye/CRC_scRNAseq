library(corrplot)
library(survival)
library(ggpubr)
library(maxstat)
library(magrittr)
options(stringsAsFactors = F)
my.cor.test <- function(...) {
  obj<-try(cor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
setwd("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/MajorTypes/Tumor/")
ClinicalData <- read.delim("/home/public/data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt",stringsAsFactors = F)

GEO_CIBERSORTID <- openxlsx::read.xlsx("TumorMajorTypesImputionOutputList.xlsx")
GEO_CIBERSORTList <- list.files(pattern = "*_Results.xlsx")
#GEO_CIBERSORTList <- GEO_CIBERSORTList[ gsub("CIBERSORTx_Job|_Results.xlsx",'',GEO_CIBERSORTList) %in% GEO_CIBERSORTID$JobID]
GEO_CIBERSORTResult <- lapply(GEO_CIBERSORTList,function(x){
  openxlsx::read.xlsx(x)
})
names(GEO_CIBERSORTResult) <- gsub("_InputCIBERSORT",'',GEO_CIBERSORTID$DatasetID)[match(gsub("CIBERSORTx_Job|_Results.xlsx","",GEO_CIBERSORTList),GEO_CIBERSORTID$JobID)]

GEO_CIBERSORT_Tumor <- tibble::tibble(Dataset = names(GEO_CIBERSORTResult),
                                      Infiltation = GEO_CIBERSORTResult)
TCGA_CRC <- tibble::tibble(Dataset = c("COAD","READ"),
                           Infiltation = lapply(c("COAD","READ"),function(x){
                             GEO_CIBERSORT_Tumor$Infiltation[[1]][substr(GEO_CIBERSORT_Tumor$Infiltation[[1]]$Mixture,1,12) %in% gsub("-",".",ClinicalData[ClinicalData$type %in% x,]$bcr_patient_barcode),]
                           }))
GEO_CIBERSORT_Tumor <- rbind(GEO_CIBERSORT_Tumor,TCGA_CRC)
GEO_CIBERSORT_Tumor  <- GEO_CIBERSORT_Tumor%>% dplyr::mutate(MutualCorr = purrr::map(.x = Infiltation, function(.x){
  .x <- .x[,1:10]
  .x %>% tidyr::gather(key=DefineTypes,value=Infiltration,-Mixture) %>% data.frame() -> .dx
  lapply(split(.dx,.dx$DefineTypes), function(y){
    y %>% dplyr::inner_join(.x,by="Mixture") %>% data.frame -> dy
    colnames(dy) <- c(colnames(y),colnames(.x)[2:10])
    Corr = data.frame(DefineTypes =colnames(.x)[2:10],Marker2 = rep(unique(y$DefineTypes),times=9))
    Corr[,c("estimate.rho","p.value")] <- t(apply(dy[,colnames(.x)[2:10]],2,function(x)unlist(my.cor.test(as.numeric(x),as.numeric(dy$Infiltration),method="spearman"))))
    Corr$FDR <- p.adjust(Corr$p.value,method="fdr")
    Corr$Class <- ifelse(Corr$estimate.rho > 0.3 & Corr$p.value < 0.05,1,ifelse(Corr$estimate.rho < -0.3 & Corr$p.value < 0.05,-1,0))
    return(Corr)
  }) %>% dplyr::bind_rows() -> AllCorr
  AllCorr$DefineTypes <- gsub("\\."," ",AllCorr$DefineTypes)
  AllCorr$Marker2 <- gsub("\\."," ",AllCorr$Marker2)
  return(AllCorr)
})) 
GEO_CIBERSORT_Tumor <- GEO_CIBERSORT_Tumor[unlist(lapply(GEO_CIBERSORT_Tumor$Infiltation,nrow))>=35,]
GEO_CIBERSORT_Tumor <- GEO_CIBERSORT_Tumor[GEO_CIBERSORT_Tumor$Dataset != "GSE37364",]

GEO_CIBERSORT_Tumor <- GEO_CIBERSORT_Tumor[!duplicated(GEO_CIBERSORT_Tumor$Dataset),]

GEO_CIBERSORT_MSC <- lapply(GEO_CIBERSORT_Tumor$MutualCorr,function(x){
  x <- x[x$DefineTypes %in% "MSC",]
  x <- x[order(x$estimate.rho),]
  # match("MSC",x$Marker2) - match("Myeloid",x$Marker2)
})

GEO_CIBERSORT_Myeloid <- lapply(GEO_CIBERSORT_Tumor$MutualCorr,function(x){
  x <- x[x$DefineTypes %in% "Myeloid",]
  x <- x[order(x$estimate.rho),]
  #match("Myeloid",x$Marker2) - match("MSC",x$Marker2)
}) 

lapply(GEO_CIBERSORT_Tumor$MutualCorr,function(x){
  x <- x[x$DefineTypes %in% "Myeloid",]
  x <- x[order(x$estimate.rho),]
  match("Myeloid",x$Marker2) - match("MSC",x$Marker2)
}) %>% unlist 

lapply(GEO_CIBERSORT_Tumor$MutualCorr,function(x){
  x[x$DefineTypes %in% "Myeloid" & x$Marker2 %in% "MSC",]
}) 

MSC_Myeloid_Cor <- lapply(GEO_CIBERSORT_Tumor$Dataset[2:length(GEO_CIBERSORT_Tumor$Dataset)],function(i){
  sub <- GEO_CIBERSORT_Tumor[GEO_CIBERSORT_Tumor$Dataset %in% i,]$Infiltation[[1]]
  colnames(sub) <- gsub("\\."," ",colnames(sub))
  sub_MSC_MACRO_p  <- signif(as.numeric(unlist(cor.test(sub$`MSC`,sub$`Myeloid`,method="spearman"))[c("estimate.rho","p.value")]),digits = 2) 
  
  p <- ggplot(  sub,aes(x=`MSC`,y=`Myeloid`))+
    geom_point(color="#FC8D62")+
    geom_smooth(method="lm",color="red")+labs(y="Myeloid Cells",x="Fibroblast")+
    labs(title=paste(i,"\n","Rs = ", sub_MSC_MACRO_p[1],"; p = ", sub_MSC_MACRO_p[2],sep=""))+
    theme(panel.background = element_rect(fill = NA, color = "black"),panel.grid = element_blank(),
          axis.text = element_text(color="black"))
  return(p)
})

pdf("../MSC_Myeloid_Cor_Scatterplot.pdf",width = 14,height = 9.2)
cowplot::plot_grid(plotlist =MSC_Myeloid_Cor ,ncol=5)
dev.off()
readr::write_rds(GEO_CIBERSORT_Tumor,path="../GEO_CIBERSORT_Tumor.rds.gz",compress = "gz")
GEO_CIBERSORT_Tumor <- readr::read_rds("../GEO_CIBERSORT_Tumor.rds.gz")
GEO_CIBERSORT_Tumor <- GEO_CIBERSORT_Tumor[rev(order(unlist(lapply(GEO_CIBERSORT_Tumor$Infiltation,nrow)))),]
paste(paste(GEO_CIBERSORT_Tumor$Dataset," (n = ",unlist(lapply(GEO_CIBERSORT_Tumor$Infiltation,nrow)),")",sep=""),collapse = "; ")

##
TT <- GEO_CIBERSORT_Tumor$MutualCorr[2:15]
GEO_CIBERSORT_TumorCorrCountAll <- dplyr::bind_rows(TT) %>% 
  dplyr::mutate(Dataset = rep(TT$Dataset,times=unlist(lapply(TT$MutualCorr,nrow))) )
GEO_CIBERSORT_TumorCorrCountAll$Class <- ifelse(is.na(GEO_CIBERSORT_TumorCorrCountAll$Class),0,GEO_CIBERSORT_TumorCorrCountAll$Class)
GEO_CIBERSORT_TumorCorrCountAll$Class <- ifelse(GEO_CIBERSORT_TumorCorrCountAll$DefineTypes == GEO_CIBERSORT_TumorCorrCountAll$Marker2,1,GEO_CIBERSORT_TumorCorrCountAll$Class)

GEO_CIBERSORT_TumorCorrCount <- t(sapply(split(GEO_CIBERSORT_TumorCorrCountAll[,"Class"],list(GEO_CIBERSORT_TumorCorrCountAll$DefineTypes,GEO_CIBERSORT_TumorCorrCountAll$Marker2)),function(x)table(factor(x,levels=c(-1,0,1))))) %>% 
  data.frame() %>% set_colnames(c("NegCor","NonCor","PosCor"))

GEO_CIBERSORT_TumorCorrCount[,c("Marker1","Marker2")] <- apply(data.frame(do.call(rbind,strsplit(rownames(GEO_CIBERSORT_TumorCorrCount),"\\."))),2,function(x){gsub("_"," ",gsub("__"," ",x))})
GEO_CIBERSORT_TumorCorrCountm <- reshape2::melt(GEO_CIBERSORT_TumorCorrCount,id.vars=c("Marker1","Marker2"),measure.vars=c("NegCor","NonCor","PosCor"))
##
LABEL=c("Epithelial", "Fibroblast", "Endothelial", "Glia Cells", "Myeloid Cells", "Mast Cells", "T Cells", "B Cells","Plasma Cells")
Order_MajorTypes <- c("Epithelial","MSC","Endo",'Glia',"Myeloid","Mast","T","B","PlasmaB")
GEO_CIBERSORT_TumorCorrCountm$Marker1 <- factor(GEO_CIBERSORT_TumorCorrCountm$Marker1,levels=Order_MajorTypes)
GEO_CIBERSORT_TumorCorrCountm$Marker2 <- factor(GEO_CIBERSORT_TumorCorrCountm$Marker2,levels=Order_MajorTypes)

CPpairs <- data.frame(Marker1=Order_MajorTypes,Marker2=Order_MajorTypes)
for(i in Order_MajorTypes){
  for(j in Order_MajorTypes){
    aa <- data.frame(Marker1=i,Marker2=j)
    if(!(paste(aa$Marker1,aa$Marker2,sep="_") %in% c(paste(CPpairs$Marker1,CPpairs$Marker2,sep="_"),paste(CPpairs$Marker2,CPpairs$Marker1,sep="_")))){
      CPpairs=rbind(CPpairs,aa)
    }}}
#CPpairs <- CPpairs[CPpairs$Marker1 != CPpairs$Marker2,]
GEO_CIBERSORT_TumorCorrCountm <- merge(CPpairs,GEO_CIBERSORT_TumorCorrCountm,by=c("Marker1","Marker2"))
GEO_CIBERSORT_TumorCorrCountm$Marker1 <- factor(GEO_CIBERSORT_TumorCorrCountm$Marker1,levels=Order_MajorTypes,labels = LABEL)
GEO_CIBERSORT_TumorCorrCountm$Marker2 <- factor(GEO_CIBERSORT_TumorCorrCountm$Marker2,levels=rev(Order_MajorTypes),labels = rev(LABEL))

pdf("../MajorTypes Correaltion_UP2.pdf",width=5,height = 5)
ggplot(data = GEO_CIBERSORT_TumorCorrCountm, aes(x = factor(1), y = value, fill = factor(variable))) +
  geom_bar(stat="identity",position = "stack", color = NA) +
  scale_y_continuous(limits = c(0,14),expand=c(0,0)) +
  coord_polar("y") +
  facet_grid(Marker2~Marker1) +
  theme(axis.text=element_blank(),axis.title = element_blank(), 
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle =0,hjust=0,color="black",size=8),strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,angle = 90,vjust = 0),legend.text = element_text(size=10),
        legend.position = "bottom",panel.spacing  = unit(0.01, "lines")) +
  scale_fill_manual(limits=c("NegCor","NonCor","PosCor"),values=c( "#377EB8","lightgray","#E41A1C"))
dev.off()


#### Tumor and normal cell types difference
GEO_CIBERSORTResultA <- GEO_CIBERSORTResult
GEO_CIBERSORTResultA$COAD <- CIBERSORTx_TCGAsum$Infiltation[[1]]
GEO_CIBERSORTResultA$READ <- CIBERSORTx_TCGAsum$Infiltation[[2]]

DatasetWithNormal <- intersect(GEO_CIBERSORT_Tumor$Dataset,gsub("_Normal","",grep("_Normal",names(GEO_CIBERSORTResult),value=T)))
DatasetWithNT <- lapply(DatasetWithNormal,function(DataName){
  sub <- GEO_CIBERSORTResultA[grep(DataName,names(GEO_CIBERSORTResultA),value=T)]
  subA <- dplyr::bind_rows(sub)
  subA <- subA %>% dplyr::mutate(Class=rep(names(sub),times=unlist(lapply(sub,nrow))))
  subA$Class <- ifelse(subA$Class %in% grep("_Normal",names(sub),value=T),"Normal","Tumor")
  return(subA)
})
names(DatasetWithNT) <- DatasetWithNormal
DatasetWithNT <- tibble::tibble(Dataset=names(DatasetWithNT),NTexp = DatasetWithNT)

COAD_Inf <- DatasetWithNT$NTexp[[6]]
colnames(COAD_Inf) <- gsub("\\."," ",colnames(COAD_Inf))
COAD_Inf_Plot <- lapply(DiffCellTypes[DiffCellTypes$Class != "Both",]$DefineTypes,function(CellName){
  sub <- COAD_Inf[,c("Mixture",CellName,"Class")] %>% set_colnames(c("Mixture","Infiltration","Class")) %>% 
    dplyr::mutate(Mixture = substr(Mixture,1,12)) %>% tidyr::spread(Class,Infiltration) %>% 
    dplyr::filter(!is.na(Normal)) %>% tidyr::gather(key=Class,value=Infiltration,-Mixture)
  subm <- COAD_Inf[,c("Mixture",CellName,"Class")] %>% set_colnames(c("Mixture","Infiltration","Class")) %>% 
    dplyr::mutate(Mixture = substr(Mixture,1,12)) %>% tidyr::spread(Class,Infiltration) %>% 
    dplyr::filter(!is.na(Normal))
  p <- ggplot(sub, aes(x = Class, y=Infiltration,group= Mixture)) + #log2(as.numeric(CopyNumber)))
    geom_point(aes(color=Class),size=2)+
    scale_x_discrete(limit=c("Normal","Tumor"),labels=c("N","T"))+
    scale_y_continuous(expand=c(0.1,0))+
    scale_color_manual(limit=c("Normal","Tumor"),values=c("blue","red"),name="",guide=F)+ #l(7,"YlOrBr")[3:7]
    geom_segment(data=subm,aes(x=1,xend=2,y=Normal,yend=Tumor),linetype = "dashed",color="gray")+#stat_compare_means(label = "p.format")+
    labs(title=CellName)+
    theme(axis.text=element_text(color = "black"),
          panel.background = element_blank(),panel.grid=element_blank(),
          legend.title=element_blank(), axis.text.y=element_text(color = "black"),
          strip.background = element_blank(),strip.text = element_text(size=12),
          axis.title = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))+
    ggpubr::stat_compare_means(comparisons =list(V1=c("Normal","Tumor")),label= "p.format",method = "wilcox.test",paired =T )
  return(p)
})
names(COAD_Inf_Plot) <- DiffCellTypes[DiffCellTypes$Class != "Both",]$DefineTypes
pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputedGEO/PairedUpCellTypes.pdf",width=12,height = 5)
cowplot::plot_grid(plotlist = COAD_Inf_Plot[DiffCellTypes[DiffCellTypes$Class %in% "Tumor",]$DefineTypes],ncol=6)
dev.off()
pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputedGEO/PairedDownCellTypes.pdf",width=10,height = 5)
cowplot::plot_grid(plotlist = COAD_Inf_Plot[DiffCellTypes[DiffCellTypes$Class %in% "Normal",]$DefineTypes],ncol=5)
dev.off()


COAD_Inf_Plot_Box <- lapply(DiffCellTypes[DiffCellTypes$Class != "Both",]$DefineTypes,function(CellName){
  sub <- COAD_Inf[,c("Mixture",CellName,"Class")] %>% set_colnames(c("Mixture","Infiltration","Class")) %>% 
    dplyr::mutate(Mixture = substr(Mixture,1,12)) 
  p <- ggplot(sub, aes(x = Class, y=Infiltration)) + #log2(as.numeric(CopyNumber)))
    geom_boxplot(aes(color=Class),outlier.shape = NA)+
    scale_x_discrete(limit=c("Normal","Tumor"),labels=c("N","T"))+
    scale_y_continuous(expand=c(0.1,0))+
    scale_color_manual(limit=c("Normal","Tumor"),values=c("blue","red"),name="",guide=F)+ #l(7,"YlOrBr")[3:7]
    labs(title=CellName)+
    theme(axis.text=element_text(color = "black"),
          panel.background = element_blank(),panel.grid=element_blank(),
          legend.title=element_blank(), axis.text.y=element_text(color = "black"),
          strip.background = element_blank(),strip.text = element_text(size=12),
          axis.title = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))+
    ggpubr::stat_compare_means(comparisons =list(V1=c("Normal","Tumor")),label= "p.format",method = "wilcox.test" )
  return(p)
})
names(COAD_Inf_Plot_Box) <- DiffCellTypes[DiffCellTypes$Class != "Both",]$DefineTypes
pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputedGEO/UpCellTypes_Box.pdf",width=12,height = 5)
cowplot::plot_grid(plotlist = COAD_Inf_Plot_Box[DiffCellTypes[DiffCellTypes$Class %in% "Tumor",]$DefineTypes],ncol=6)
dev.off()
pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputedGEO/DownCellTypes_Box.pdf",width=10,height = 5)
cowplot::plot_grid(plotlist = COAD_Inf_Plot_Box[DiffCellTypes[DiffCellTypes$Class %in% "Normal",]$DefineTypes],ncol=5)
dev.off()


DiffCellTypes <- openxlsx::read.xlsx("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/DiffCellTypes_Normal_Tumor.xlxs")

READ_Inf <- DatasetWithNT$NTexp[[7]]
colnames(READ_Inf) <- gsub("\\."," ",colnames(READ_Inf))
READ_Inf_Plot <- lapply(DiffCellTypes[DiffCellTypes$Class != "Both",]$DefineTypes,function(CellName){
  sub <- READ_Inf[,c("Mixture",CellName,"Class")] %>% set_colnames(c("Mixture","Infiltration","Class")) %>% 
    dplyr::mutate(Mixture = substr(Mixture,1,12)) %>% tidyr::spread(Class,Infiltration) %>% 
    dplyr::filter(!is.na(Normal)) %>% tidyr::gather(key=Class,value=Infiltration,-Mixture)
  subm <- READ_Inf[,c("Mixture",CellName,"Class")] %>% set_colnames(c("Mixture","Infiltration","Class")) %>% 
    dplyr::mutate(Mixture = substr(Mixture,1,12)) %>% tidyr::spread(Class,Infiltration) %>% 
    dplyr::filter(!is.na(Normal))
  p <- ggplot(sub, aes(x = Class, y=Infiltration,group= Mixture)) + #log2(as.numeric(CopyNumber)))
    geom_point(aes(color=Class),size=2)+
    scale_x_discrete(limit=c("Normal","Tumor"),labels=c("N","T"))+
    scale_y_continuous(expand=c(0.1,0))+
    scale_color_manual(limit=c("Normal","Tumor"),values=c("blue","red"),name="",guide=F)+ #l(7,"YlOrBr")[3:7]
    geom_segment(data=subm,aes(x=1,xend=2,y=Normal,yend=Tumor),linetype = "dashed",color="gray")+#stat_compare_means(label = "p.format")+
    labs(title=CellName)+
    theme(axis.text=element_text(color = "black"),
          panel.background = element_blank(),panel.grid=element_blank(),
          legend.title=element_blank(), axis.text.y=element_text(color = "black"),
          strip.background = element_blank(),strip.text = element_text(size=12),
          axis.title = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))+
    ggpubr::stat_compare_means(comparisons =list(V1=c("Normal","Tumor")),label= "p.format",method = "wilcox.test",paired =T )
  return(p)
})
names(READ_Inf_Plot) <- DiffCellTypes[DiffCellTypes$Class != "Both",]$DefineTypes
pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputedGEO/READ_PairedUpCellTypes.pdf",width=12,height = 5)
cowplot::plot_grid(plotlist = READ_Inf_Plot[DiffCellTypes[DiffCellTypes$Class %in% "Tumor",]$DefineTypes],ncol=6)
dev.off()
pdf("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputedGEO/READ_PairedDownCellTypes.pdf",width=10,height = 5)
cowplot::plot_grid(plotlist = READ_Inf_Plot[DiffCellTypes[DiffCellTypes$Class %in% "Normal",]$DefineTypes],ncol=5)
dev.off()
