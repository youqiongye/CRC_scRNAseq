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
setwd("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/MajorTypes/")
ClinicalData <- read.delim("/home/public/data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt",stringsAsFactors = F)
#WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")

COAD_READ_Clinic <- ClinicalData[ClinicalData$type %in% c("COAD","READ"),]
COAD_Clinic <- ClinicalData[ClinicalData$type %in% c("COAD"),]
READ_Clinic <- ClinicalData[ClinicalData$type %in% c("READ"),]

GEO_CIBERSORT_Normal<- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/MajorTypes//GEO_CIBERSORT_Normal.rds.gz")
GEO_CIBERSORT_Tumor <- readr::read_rds("/home/yye/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/MajorTypes//GEO_CIBERSORT_Tumor.rds.gz")
###TCGA
TCGA_Normal <- GEO_CIBERSORT_Normal[GEO_CIBERSORT_Normal$Dataset %in% "COAD_READ",]$Infiltation[[1]]
TCGA_Tumor <- GEO_CIBERSORT_Tumor[GEO_CIBERSORT_Tumor$Dataset %in% "COAD_READ",]$Infiltation[[1]]

NT_Inf <- tibble::tibble(Cancer_types=c("COAD_READ","COAD","READ"),Clinic=list(COAD_READ_Clinic,COAD_Clinic,READ_Clinic))
NT_Inf <- NT_Inf %>% dplyr::mutate(PairedNT_Inf=purrr::map(.x=Clinic,function(.x){
  Normal <- TCGA_Normal %>% dplyr::mutate(Mixture=substr(gsub("\\.","-",Mixture),1,12)) %>%
    dplyr::filter(Mixture %in% .x$bcr_patient_barcode)
  Tumor <- TCGA_Tumor %>% dplyr::mutate(Mixture=substr(gsub("\\.","-",Mixture),1,12)) %>%
    dplyr::filter(Mixture %in% Normal$Mixture )
  PairedNT <- dplyr::bind_rows(Normal,Tumor) %>% dplyr::mutate(Class=rep(c("N","T"),times=c(nrow(Normal),nrow(Tumor)))) %>%
    dplyr::filter(Mixture %in% intersect(Normal$Mixture,Tumor$Mixture))
  return(PairedNT)
}))
Names <- colnames(NT_Inf$PairedNT_Inf[[1]])[2:10]

NT_Inf <- NT_Inf %>% dplyr::mutate(PairNT_Plot=purrr::map(.x=PairedNT_Inf,function(.x){
  pp <- lapply(Names,function(CellName){
    sub <- .x[,c("Mixture",CellName,"Class")] %>% set_colnames(c("Mixture","Infiltration","Class")) 
    sub <- sub[!duplicated(sub[,c(1,3)]),]
    #sub <- sub %>% dplyr::mutate(Mixture = substr(Mixture,1,12)) %>% tidyr::spread(Class,Infiltration) %>% 
    #  dplyr::filter(!is.na(N)) %>% tidyr::gather(key=Class,value=Infiltration,-Mixture)
    subm <- sub %>% 
      dplyr::mutate(Mixture = substr(Mixture,1,12)) %>% tidyr::spread(Class,Infiltration) %>% 
      dplyr::filter(!is.na(N))
    p <- ggplot(sub, aes(x = Class, y=Infiltration,group= Mixture)) + #log2(as.numeric(CopyNumber)))
      geom_point(aes(color=Class),size=2)+
      scale_x_discrete(limit=c("N","T"),labels=c("N","T"))+
      scale_y_continuous(expand=c(0.1,0))+
      scale_color_manual(limit=c("N","T"),values=c("blue","red"),name="",guide=F)+ #l(7,"YlOrBr")[3:7]
      geom_segment(data=subm,aes(x=1,xend=2,y=N,yend=T),linetype = "dashed",color="gray")+#stat_compare_means(label = "p.format")+
      labs(title=CellName)+
      theme(axis.text=element_text(color = "black"),
            panel.background = element_blank(),panel.grid=element_blank(),
            legend.title=element_blank(), axis.text.y=element_text(color = "black"),
            strip.background = element_blank(),strip.text = element_text(size=12),
            axis.title = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))+
      ggpubr::stat_compare_means(comparisons =list(V1=c("N","T")),label= "p.format",method = "wilcox.test",paired =T )
    return(p)
  })
  ppp <- cowplot::plot_grid(plotlist = pp,ncol=3)
  return(ppp)
}))
lapply(NT_Inf$Cancer_types,function(x){
  pdf(paste(x,"DifferentInfiltration.pdf",sep=""),height = 9,width=7.5)
  print(NT_Inf[NT_Inf$Cancer_types %in% x,]$PairNT_Plot[[1]])
  dev.off()
})



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