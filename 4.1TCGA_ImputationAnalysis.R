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
#WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")

COAD_READ_Clinic <- ClinicalData[ClinicalData$type %in% c("COAD","READ"),]
COAD_Clinic <- ClinicalData[ClinicalData$type %in% c("COAD"),]
READ_Clinic <- ClinicalData[ClinicalData$type %in% c("READ"),]
ClustersColor <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
ClustersSubTypes <- unique(WholeTissueSC@meta.data[,c("DefineTypes","SubTypes","MainTypes","TopTypes")])
#COAD_Clinic <- read.delim("/home/yye/Data/Public/TCGA/TCGA_clinical/COAD_clinical_clean.txt")
#READ_Clinic <- read.delim("/home/yye/Data/Public/TCGA/TCGA_clinical/READ_clinical_clean.txt")
MSI_Colon <- rbind(COAD_Clinic[,c("barcode","Subtype_MSI_status")],READ_Clinic[,c("barcode","Subtype_MSI_status")])
GEO_CIBERSORT_Tumor<- readr::read_rds("../GEO_CIBERSORT_Tumor.rds.gz")


CIBERSORTx_Absolute <- GEO_CIBERSORT_Tumor[GEO_CIBERSORT_Tumor$Dataset %in% "COAD_READ",]$Infiltation[[1]]#openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/3.WholeTissue/ImputeTCGA/CIBERSORTx2/CIBERSORTx_COAD.xlsx")
colnames(CIBERSORTx_Absolute) <- gsub("\\."," ",colnames(CIBERSORTx_Absolute))

Names <- colnames(CIBERSORTx_Absolute)[2:10]
ClinicOut <- tibble::tibble(Cancer_types=c("COAD_READ","COAD","READ"),Clinic=list(COAD_READ_Clinic,COAD_Clinic,READ_Clinic))
ClinicOut <- ClinicOut %>% dplyr::mutate(gs_Best=purrr::map(.x=Clinic,function(.x){
  gs_Best <- lapply(Names ,function(i){
    Sub <- CIBERSORTx_Absolute %>% data.frame %>%  set_colnames(colnames(CIBERSORTx_Absolute)) %>% dplyr::select(c("Mixture",i)) %>%
      set_colnames(c("barcode","Por"))
    Sub$barcode <- substr(gsub("\\.","-",Sub$barcode),1,12)
    Sub_PhaseI <- .x[,c("bcr_patient_barcode","PFI","PFI.time")] %>% set_colnames(c("barcode","PFI","PFIT"))
    Sub_PhaseI <- merge(Sub_PhaseI,Sub[,c("barcode","Por")],by="barcode")
    Sub_PhaseI$PFIT <- as.numeric(Sub_PhaseI$PFIT)/30
    Sub_PhaseI <- Sub_PhaseI[!is.na(Sub_PhaseI$PFIT ),]
    Sub_PhaseI$PFI <- as.numeric(Sub_PhaseI$PFI)
    #Sub_PhaseI$Por <- as.numeric(Sub_PhaseI$Por)*100
    #aa_result <- try(aa <- maxstat.test(Surv(PFIT, PFI) ~ Por, data=Sub_PhaseI,smethod="LogRank"),silent=TRUE)
    tempresult <- try(model1 <- coxph(Surv(PFIT, PFI) ~ Por, data=Sub_PhaseI, na.action=na.exclude),silent=TRUE)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=2)
    HR_detail = summary(model1)
    CI =  paste(signif(HR_detail$conf.int[,"lower .95"],3),'-',signif(HR_detail$conf.int[,"upper .95"],3),sep=" ") #Result[i1, c('CI_95%_for_HR')]
    #Sub_PhaseI$PorGroup <- ifelse(Sub_PhaseI$Por >  0.04538792,"High","Low")
    aa_try <- try(aa <- maxstat.test(Surv(PFIT, PFI) ~ Por, data=Sub_PhaseI,smethod="LogRank"),silent=TRUE)
    if(!is(aa_try, "try-error") & median(Sub_PhaseI$Por)>0 & length(table(Sub_PhaseI$Por >= aa$estimate))==2){
      
      print(i)
      Sub_PhaseI$PorGroup <- ifelse(Sub_PhaseI$Por >= aa_try$estimate,"High","Low")
      #Sub_PhaseI$PorGroup <- ifelse(Sub_PhaseI$Por >= median(Sub_PhaseI$Por),"High","Low")
      # Sub_PhaseI$PorGroup <- ifelse(Sub_PhaseI$Por >= quantile(Sub_PhaseI$Por)[4],"High",ifelse(Sub_PhaseI$Por <= quantile(Sub_PhaseI$Por)[2],"Low","Intermediate"))
      #Sub_PhaseI <-  Sub_PhaseI[ Sub_PhaseI$Por !="Intermediate",]
      Sub_PhaseI$PorGroup <- factor(Sub_PhaseI$PorGroup,levels=c("High","Low"))
      model1 <- survdiff(Surv(PFIT, PFI) ~ PorGroup, data=Sub_PhaseI, na.action=na.exclude)
      fit <- survival::survfit(survival::Surv(PFIT, PFI) ~  PorGroup, data=Sub_PhaseI, na.action=na.exclude)
      KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(Sub_PhaseI$PorGroup)))-1)
      p <-  survminer::ggsurvplot(fit, data = Sub_PhaseI,
                                  # surv.median.line = "hv", # Add medians survival
                                  
                                  # Change legends: title & labels
                                  legend.title = paste( i,paste("HR ",HR, " [95% CI:",CI,"]",sep=""),sep="\n"),
                                  legend.labs = c(paste("High (n=",table(Sub_PhaseI$PorGroup)[1],")",sep=""), paste("Low (n=",table(Sub_PhaseI$PorGroup)[2],")",sep="")),
                                  
                                  legend = c(0.8,0.8),
                                  # Add p-value and tervals
                                  pval = TRUE,
                                  xlim = c(0,155),
                                  axes.offset = T,
                                  #conf.int = TRUE,
                                  # Add risk table
                                  risk.table = F,
                                  #tables.height = 0.3,
                                  #tables.theme = survminer::theme_cleantable(),
                                  # 
                                  # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                                  # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                                  palette = c("red","blue"),
                                  xlab="Months",ylab="Survival rate",
                                  ggtheme = theme(panel.background=element_rect(colour=NA,fill=NA,size=0.5),
                                                  panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
                                                  panel.spacing = unit(0.1,"line"),
                                                  legend.background = element_blank(),legend.key = element_blank(),
                                                  axis.title=element_blank(),
                                                  axis.text=element_text(size=10,color="black"), axis.line=element_line(),
                                                  axis.ticks.length = unit(.15, "cm"))) # Change ggplot
      return(p$plot)
    }else{
      return(NA)
    }
    
  })
  return(gs_Best)
}))

ClinicOut <- ClinicOut %>% dplyr::mutate(gs_HR_Best=purrr::map(.x=Clinic,function(.x){
  gs_HR_Best <- lapply(Names ,function(i){
    Sub <- CIBERSORTx_Absolute %>% data.frame %>%  set_colnames(colnames(CIBERSORTx_Absolute)) %>% dplyr::select(c("Mixture",i)) %>%
      set_colnames(c("barcode","Por"))
    Sub$barcode <- substr(gsub("\\.","-",Sub$barcode),1,12)
    Sub_PhaseI <- .x[,c("bcr_patient_barcode","PFI","PFI.time")] %>% set_colnames(c("barcode","PFI","PFIT"))
    Sub_PhaseI <- merge(Sub_PhaseI,Sub[,c("barcode","Por")],by="barcode")
    Sub_PhaseI$PFIT <- as.numeric(Sub_PhaseI$PFIT)/30
    Sub_PhaseI <- Sub_PhaseI[!is.na(Sub_PhaseI$PFIT ),]
    Sub_PhaseI$PFI <- as.numeric(Sub_PhaseI$PFI)
    Sub_PhaseI$Por <- as.numeric(Sub_PhaseI$Por)
    
    #aa_result <- try(aa <- maxstat.test(Surv(PFIT, PFI) ~ Por, data=Sub_PhaseI,smethod="LogRank"),silent=TRUE)
    tempresult <- try(model1 <- coxph(Surv(PFIT, PFI) ~ Por, data=Sub_PhaseI, na.action=na.exclude),silent=TRUE)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
    HR_detail = summary(model1)
    CILow =  HR_detail$conf.int[,"lower .95"]
    CIHigh =  HR_detail$conf.int[,"upper .95"]
    
    CI =  paste(signif(HR_detail$conf.int[,"lower .95"],3),'-',signif(HR_detail$conf.int[,"upper .95"],3),sep=" ")   #aa_result <- try(aa <- maxstat.test(Surv(PFIT, PFI) ~ Por, data=Sub_PhaseI,smethod="LogRank"),silent=TRUE)
    #print(CI)
    aa_try <- try(aa <- maxstat.test(Surv(PFIT, PFI) ~ Por, data=Sub_PhaseI,smethod="LogRank"),silent=TRUE)
    if(!is(aa_try, "try-error") & median(Sub_PhaseI$Por)>0 & length(table(Sub_PhaseI$Por >= aa$estimate))==2){
      print(i)
      Sub_PhaseI$PorGroup <- ifelse(Sub_PhaseI$Por >= aa_try$estimate,"High","Low")
      #Sub_PhaseI$PorGroup <- ifelse(Sub_PhaseI$Por >= median(Sub_PhaseI$Por),"High","Low")
      # Sub_PhaseI$PorGroup <- ifelse(Sub_PhaseI$Por >= quantile(Sub_PhaseI$Por)[4],"High",ifelse(Sub_PhaseI$Por <= quantile(Sub_PhaseI$Por)[2],"Low","Intermediate"))
      #Sub_PhaseI <-  Sub_PhaseI[ Sub_PhaseI$Por !="Intermediate",]
      Sub_PhaseI$PorGroup <- factor(Sub_PhaseI$PorGroup,levels=c("High","Low"))
      model1 <- survdiff(Surv(PFIT, PFI) ~ PorGroup, data=Sub_PhaseI, na.action=na.exclude)
      fit <- survival::survfit(survival::Surv(PFIT, PFI) ~  PorGroup, data=Sub_PhaseI, na.action=na.exclude)
      KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(Sub_PhaseI$PorGroup)))-1)
      aa <- data.frame(DefineTypes=i,HR,CILow,CIHigh,CI,Coxp=HR_detail$coefficients[5],KMP)
      aa$Group <- ifelse(aa$HR >=1,"High","Low")
      return(aa)
      
    }else{
      return(NA)
    }
  })
  
  return(gs_HR_Best)
}))
readr::write_rds(ClinicOut,path="../TCGA_ImputationPFI.rds.gz",compress = "gz")                                         


gs_Best1 <- ClinicOut$gs_Best[[1]][unlist(lapply(ClinicOut$gs_HR_Best[[1]] ,function(x)is.data.frame(x)))]
names(gs_Best1) <- Names[unlist(lapply(ClinicOut$gs_HR_Best[[1]] ,function(x)is.data.frame(x)))]
pdf("../TCGA_MajorTypesPFIBestP.pdf",width = 9,height = 6)
cowplot::plot_grid(plotlist = gs_Best1,ncol = 3)
dev.off()

ClinicOut$gs_HR_Best[[1]] %>% dplyr::bind_rows(.)

###Other clinical features


COAD_Clinic <- read.delim("/home/yye/Data/Public/TCGA/TCGA_clinical/COAD_clinical_clean.txt")
READ_Clinic <- read.delim("/home/yye/Data/Public/TCGA/TCGA_clinical/READ_clinical_clean.txt")
MSI_Colon <- rbind(COAD_Clinic[,c("barcode","Subtype_MSI_status")],READ_Clinic[,c("barcode","Subtype_MSI_status")])

COAD_READ_ClinicF <- COAD_READ_Clinic[,c("bcr_patient_barcode","type","gender","ajcc_pathologic_tumor_stage")] %>% set_colnames(c("barcode","type","gender","stage"))
#COAD_READ_ClinicF <- dplyr::inner_join(COAD_READ_ClinicF,MSI_Clinic[,c("Sample.Name","")],by=c("barcode"="Sample.Name"))



PropTumorSubPor <- CIBERSORTx_Absolute  %>% 
  dplyr::mutate(barcode=gsub("\\.","-",substr(Mixture,1,12))) %>%dplyr::select(c("barcode",Names)) %>% 
  dplyr::inner_join(COAD_READ_ClinicF,by="barcode") %>% dplyr::mutate(stage=gsub("[ABC]$","",stage)) %>%
  dplyr::mutate(stage = ifelse(!(stage %in% c("Stage I","Stage II","Stage III","Stage IV")),"No Avail",stage)) %>%
  dplyr::filter(!(stage %in% "No Avail"))
PropTumorSubPor$class <- ifelse(PropTumorSubPor$stage %in% c("Stage I","Stage II"),"Early","Late")

PropTumorSubPorM <- PropTumorSubPor %>% tidyr::gather(key=MajorTypes,value=Por,Names)
StageCandidateExpStat <- sapply(split(PropTumorSubPorM,PropTumorSubPorM$MajorTypes),function(x){
  Med <- sapply(split(x[,"Por"],x$stage),median)
  MaxMed <- names(Med[Med==max(Med)])[1]
  MinMed <- names(Med[Med==min(Med)])[1]
  # MaxFC <- mean(log2(x[x$stage==MaxMed,"value"]+1)) - mean(log2(x[x$stage==MinMed,"value"]+1))
  Result <- c(Stage_pval = signif(oneway.test(log2(Por+1)~stage,data = x)$p.value,digits = 2),
              EarlyLate_pval = signif(t.test(log2(Por+1)~class,data = x)$p.value,digits = 2),MaxMedValue=max(Med),MaxStage = MaxMed,
              Max_FC = mean(log2(x[x$stage==MaxMed,"Por"]+1)) - mean(log2(x[x$stage==MinMed,"Por"]+1)),
              EarlyLate_FC = mean(log2(x[which(x$class=="Late"),"Por"]+1)) - mean(log2(x[which(x$class=="Early"),"Por"]+1))) 
  return(Result) 
}) %>% t() %>% data.frame %>% dplyr::mutate(MajorTypes = row.names(.))
StageCandidateExpStat$Stage_pval <- as.numeric(StageCandidateExpStat$Stage_pval)
StageCandidateExpStat$Stage_FDR <- p.adjust(StageCandidateExpStat$Stage_pval,method="fdr")
StageCandidateExpStat$EarlyLate_FDR <- p.adjust(StageCandidateExpStat$EarlyLate_pval,method="fdr")
StageCandidateExpStatSig <- StageCandidateExpStat[StageCandidateExpStat$Stage_pval <0.05 & StageCandidateExpStat$MaxMedValue > 0,]
StageCandidateExpStatSig  <- StageCandidateExpStatSig [!is.na(StageCandidateExpStatSig$Stage_pval),]
StageCandidateExpStatSig <- StageCandidateExpStatSig[order(StageCandidateExpStatSig$Stage_pval),]
PropTumorSubPorMSig <- PropTumorSubPorM[PropTumorSubPorM$MajorTypes %in% StageCandidateExpStatSig$MajorTypes,]
PropTumorSubPorMSig$MajorTypes <- factor(PropTumorSubPorMSig$MajorTypes,levels=StageCandidateExpStatSig$MajorTypes)


pdf("../2.MajorTypes Absolute proportion in stage.pdf",width = 9,height = 5)
ggplot(PropTumorSubPorMSig,aes(x=stage,y=Por))+
  geom_boxplot(aes(fill=stage),outlier.shape = NA,width=0.5)+
  facet_wrap(~MajorTypes,scale="free",nrow=2)+ stat_compare_means(label="p.format")+
  scale_x_discrete(labels=c("I","II","III","IV"))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Set2"))+labs(x="Stage",y="Percentage(%)")+
  theme(panel.background=element_rect(fill = "white",color="black"),strip.background = element_blank(),
        panel.grid =element_blank(),axis.text = element_text(color="black"))
dev.off()



PropTumorSubPorM_MSI <- merge(PropTumorSubPorM,MSI_Colon,by="barcode")
PropTumorSubPorM_MSI <- PropTumorSubPorM_MSI[!is.na(PropTumorSubPorM_MSI$Subtype_MSI_status),]
PropTumorSubPorM_MSI$MajorTypes <- factor(PropTumorSubPorM_MSI$MajorTypes,levels=Names)
pdf("../2.MajorTypes expression in MSI status.pdf",width = 12,height = 3)
ggplot(PropTumorSubPorM_MSI,aes(x=Subtype_MSI_status,y=Por*100))+
  geom_boxplot(aes(fill=Subtype_MSI_status),outlier.shape = NA,width=0.5)+
  facet_wrap(~MajorTypes,nrow=1,scales = "free")+stat_compare_means(label="p.format")+
  scale_x_discrete(labels=c("MSI-H", "MSI-L","MSS"))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Set1"))+labs(x="MSS/MSI status",y="Infiltration")+
  theme(panel.background=element_rect(fill = "white",color="black"),strip.background = element_blank(),
        panel.grid =element_blank(),axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black",angle=90))
dev.off()

###Cross

Myeloid_MSC <- CIBERSORTx_Absolute %>% data.frame %>%  set_colnames(colnames(CIBERSORTx_Absolute)) %>% 
  dplyr::select(c("Mixture", "MSC","Myeloid")) %>%
  set_colnames(c("barcode","MSC","Myeloid"))

Myeloid_MSC$barcode <- substr(gsub("\\.","-",Myeloid_MSC$barcode),1,12)
Myeloid_MSC_PhaseI <- COAD_READ_Clinic[,c("bcr_patient_barcode","PFI","PFI.time")] %>% set_colnames(c("barcode","PFI","PFIT"))
Myeloid_MSC_PhaseI <- merge(Myeloid_MSC_PhaseI,Myeloid_MSC,by="barcode")

Myeloid_MSC_PhaseI$PFIT <- as.numeric(Myeloid_MSC_PhaseI$PFIT)/30
Myeloid_MSC_PhaseI <- Myeloid_MSC_PhaseI[!is.na(Myeloid_MSC_PhaseI$PFIT ),]
Myeloid_MSC_PhaseI$PFI <- as.numeric(Myeloid_MSC_PhaseI$PFI)
Myeloid_MSC_PhaseI$MSC_G <- ifelse(Myeloid_MSC_PhaseI$MSC > median(Myeloid_MSC_PhaseI$MSC),"High","Low")
Myeloid_MSC_PhaseI$Myeloid_G <- ifelse(Myeloid_MSC_PhaseI$Myeloid > median(Myeloid_MSC_PhaseI$Myeloid),"High","Low")
Myeloid_MSC_PhaseI$Group <- paste("MSC.",Myeloid_MSC_PhaseI$MSC_G,"_","Myeloid.",Myeloid_MSC_PhaseI$Myeloid_G,sep="")
model1 <- survival::coxph(survival::Surv(PFIT, PFI) ~  Group, data=Myeloid_MSC_PhaseI, na.action=na.exclude)
fit <- survival::survfit(survival::Surv(PFIT, PFI) ~  Group, data=Myeloid_MSC_PhaseI, na.action=na.exclude)
KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(Myeloid_MSC_PhaseI$Group)))-1)
pdf("MSC_Myeloid cell median infiltration PFS.pdf",width = 5,height = 5)
p <- survminer::ggsurvplot(fit, data = Myeloid_MSC_PhaseI,
                           legend.title="",
                           legend.labs = paste(names(table(Myeloid_MSC_PhaseI$Group)),"(n = ",table(Myeloid_MSC_PhaseI$Group),")",sep=""),
                           legend = c(0.7,0.8),
                           pval = TRUE,
                           xlim = c(0,155),
                           axes.offset = T,
                           #conf.int = TRUE,
                           # Add risk table
                           risk.table = F,
                           palette = c("red","orange","cyan","blue"),
                           xlab="Months",ylab="Survival rate",
                           ggtheme = theme(panel.background=element_rect(colour=NA,fill=NA,size=0.5),
                                           panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
                                           panel.spacing = unit(0.1,"line"),
                                           legend.background = element_blank(),legend.key = element_blank(),
                                           axis.title=element_blank(),
                                           axis.text=element_text(size=10,color="black"), axis.line=element_line(),
                                           axis.ticks.length = unit(.15, "cm"))) # Change ggplot
print(p$plot)
dev.off()


aa_MSC <- try(aa <- maxstat.test(Surv(PFIT, PFI) ~ MSC, data=Myeloid_MSC_PhaseI,smethod="LogRank"),silent=TRUE)
aa_Myeloid <- try(aa <- maxstat.test(Surv(PFIT, PFI) ~ Myeloid, data=Myeloid_MSC_PhaseI,smethod="LogRank"),silent=TRUE)
Myeloid_MSC_PhaseI$MSC_BG <- ifelse(Myeloid_MSC_PhaseI$MSC > aa_MSC$estimate,"High","Low")
Myeloid_MSC_PhaseI$Myeloid_BG <- ifelse(Myeloid_MSC_PhaseI$Myeloid > aa_Myeloid$estimate,"High","Low")
Myeloid_MSC_PhaseI$GroupB <- paste("MSC.",Myeloid_MSC_PhaseI$MSC_BG,"_","Myeloid.",Myeloid_MSC_PhaseI$Myeloid_BG,sep="")
model1 <- survival::coxph(survival::Surv(PFIT, PFI) ~  GroupB, data=Myeloid_MSC_PhaseI, na.action=na.exclude)

fit <- survival::survfit(survival::Surv(PFIT, PFI) ~  GroupB, data=Myeloid_MSC_PhaseI, na.action=na.exclude)
KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(Myeloid_MSC_PhaseI$GroupB)))-1)
pdf("MSC_Myeloid cell best Cutoff infiltration PFS.pdf",width = 5,height = 5)
p <- survminer::ggsurvplot(fit, data = Myeloid_MSC_PhaseI,
                           legend.title="",
                           legend.labs = paste(names(table(Myeloid_MSC_PhaseI$GroupB)),"(n = ",table(Myeloid_MSC_PhaseI$GroupB),")",sep=""),
                           legend = c(0.7,0.8),
                           pval = TRUE,
                           xlim = c(0,155),
                           axes.offset = T,
                           #conf.int = TRUE,
                           # Add risk table
                           risk.table = F,
                           palette = c("red","orange","cyan","blue"),
                           xlab="Months",ylab="Survival rate",
                           ggtheme = theme(panel.background=element_rect(colour=NA,fill=NA,size=0.5),
                                           panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
                                           panel.spacing = unit(0.1,"line"),
                                           legend.background = element_blank(),legend.key = element_blank(),
                                           axis.title=element_blank(),
                                           axis.text=element_text(size=10,color="black"), axis.line=element_line(),
                                           axis.ticks.length = unit(.15, "cm"))) # Change ggplot
print(p$plot)
dev.off()

