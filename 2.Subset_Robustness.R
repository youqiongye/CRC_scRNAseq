library(multiROC)
library(Seurat)
library(randomForest)
library(ggplot2)
library(caTools)
library(magrittr)
options(expressions = 5e5)
options(stringsAsFactors = F)
setwd("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/Fig/Robustness")
Terms <- c("Endo","Epithelial","Myeloid","MSC","T","B")#

WholeTissueSC <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineUpdateB.rds.gz")
WholeTissueList <- readr::read_rds("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/WholeTissue_DefineEachSubtypeUpdataB.rds.gz")
Clusters <- openxlsx::read.xlsx("~/Project/Collaboration/2020ILC/WholeTissue/Harmony/ProcessData/ColorCode_Order.xlsx")
Terms <- c("Endo","Epithelial","Myeloid","MSC","T","B")#
plotp <- list()
for(term in Terms){
  Tg <- WholeTissueList[[term]]
  TgCluster <- Clusters[Clusters$SubTypes %in% term,]
  TgCluster$Cluster <- 1:nrow(TgCluster )
  Tg <- ScaleData(Tg)
  Tg_ScaleData <- Tg@assays$RNA@scale.data %>% t %>% data.frame
  Tg_ScaleData$Cluster <- TgCluster$Cluster[match(Tg$DefineTypes,TgCluster$DefineTypes)]
  
  #2 Get train and test data
  set.seed(42)
  split <- sample.split(Tg_ScaleData$Cluster, SplitRatio = 0.5)
  
  train_df <- subset(Tg_ScaleData, split == TRUE)
  test_df <- subset(Tg_ScaleData, split == FALSE)
  #3 random frest
  rf_res <- randomForest::randomForest(factor(Cluster) ~., data = train_df, ntree = 100)
  rf_pred <- predict(rf_res, test_df,  type = 'prob') 
  rf_pred <- data.frame(rf_pred)
  colnames(rf_pred) <- paste("Cluster.",gsub("X","",colnames(rf_pred)), " _pred_RF",sep="")
  
  
  #4 Merge true labels and predicted values
  true_label <- dummies::dummy(test_df$Cluster, sep = ".")
  true_label <- data.frame(true_label)
  colnames(true_label) <-  colnames(true_label)
  colnames(true_label) <- paste(colnames(true_label), "_true")
  final_df <- cbind(true_label,  rf_pred)
  #5 multiROC and multiPR
  roc_res <- multi_roc(final_df, force_diag=T) #This function calculates the Specificity, Sensitivity and AUC of multi-class classifications.
  pr_res <- multi_pr(final_df, force_diag=T) #his function calculates the Precision, Recall and AUC of multi-class classifications.
  #6 Plot
  plot_roc_df <- plot_roc_data(roc_res)
  plot_pr_df <- plot_pr_data(pr_res)
  
  plot_roc_df <- plot_roc_df[plot_roc_df$Group %in% paste("Cluster.",unique(test_df$Cluster)," ",sep=""), ]
  plot_roc_df_AUC <- roc_res$AUC$RF %>% unlist %>% signif(digits = 3)
  names(plot_roc_df_AUC) <- gsub(" ","",names(plot_roc_df_AUC))
  plot_roc_df$Group <- gsub(" ","",plot_roc_df$Group)
  plot_roc_df_AUC1 <- data.frame(  plot_roc_df_AUC )
  plot_roc_df_AUC1$Cluster <- rownames(  plot_roc_df_AUC1)
  plot_roc_df_AUC1 <- plot_roc_df_AUC1[!(plot_roc_df_AUC1$Cluster %in% c("macro","micro")),] %>% dplyr::mutate(x=1:nrow(.),y=1:nrow(.))
  
  p <- ggplot(plot_roc_df_AUC1, aes(x = x, y=y,color=Cluster)) +
   geom_point()+
       labs(x="False Positive Rate",y="True Positive Rate")+
    scale_color_manual(limit=paste("Cluster.",TgCluster$Cluster,sep=""),
                       values=as.character(TgCluster$Color),
                       label= paste(TgCluster$DefineTypes," (",plot_roc_df_AUC[match(paste("Cluster.",TgCluster$Cluster,sep=""),names(plot_roc_df_AUC))],")",sep=""),
                       name="AUC")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                       legend.justification=c(1, 0), legend.position=c(.95, .05), 
                       legend.background = element_rect(fill=NULL, size=0.5, 
                                                        linetype="solid", colour ="black"))
  plotp[[term]] <- p
  
  
}

pdf("ClusterRobustness.pdf",width = 12,height = 8)
p <- cowplot::plot_grid(plotlist = plotp,ncol = 3)
print(p)
dev.off()
