# Script for analyzing TCGA BRCA data
# Created by Yuanzhong Pan on Apr 3rd, 2023
#
# Data Source: UCSC Xena TGGA Breast Cancer (BRCA) Cohort :https://xenabrowser.net/datapages/?cohort=TCGA%20Breast%20Cancer%20(BRCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#      Gene Expression RNAseq: IlluminaHiSeq, renamed as "Xena_BRCA_HiSeqV2"
#      Survival: Phenotype -> Curated survival data, saved as "
#      Phenotype: Phenotype -> Phenotypes

library("ggplot2")
library("ggpubr")
library('pheatmap')
library('paletteer')
library('corrplot')
library('RColorBrewer')
library("magrittr")
library("devtools")
library("stringr")
library("tidyr")
library("survival")
library("survminer")
library("edgeR")
library("glmnet")

#1 Data preparation: read in data, initial data cleaning, and define gene sets
brca <- read.csv("./data/Xena_BRCA_HiSeqV2", row.names = 1, sep = "\t")
brca <- t(brca)
  ## Fix row names (colnames cannot take "-" and they are replaced by "." automatically. After transpose, the names can be fixed to have dash)
rownames(brca) <- gsub("\\.", "-", rownames(brca))
  ## Calculate z scores 
brca_z <- brca
for (i in 1:length(colnames(brca_z))){
  brca_z[,i] <- (brca_z[,i] - mean(brca_z[,i]))/sd(brca_z[,i])
}
  ## Separate tumor and normal tissue expression data and order the rows by sample ID
brca_norm <- brca[substring(rownames(brca),14,15) == "11",]
brca_norm <- brca_norm[order(row.names(brca_norm)),]
brca <- brca[substring(rownames(brca),14,15) == "01",]
brca <- brca[order(row.names(brca)),]

brca_norm_z <- brca_z[rownames(brca_norm),]
brca_z <- brca_z[rownames(brca),]

  ## Read survival data and only keep samples that contain expression data. This data is already ordered.
surv_brca <- read.csv("./data/survival_BRCA.txt", row.names = 1, sep = "\t")
surv_brca <- surv_brca[rownames(surv_brca) %in% rownames(brca),]
  ## subset BRCA expression data that have matched survival data. n: 1097 -> 1095
brca_surv <- brca[rownames(brca) %in% rownames(surv_brca),]
  ## Subset survival data of deceased patients 
brca_dec <- brca[surv_brca[,"OS"]==1,]
OS_dec <- surv_brca[surv_brca[,"OS"]==1,"OS.time"]

# BRCA Subtypes (Histological and PAM50)
  ## Read phenotype data and only keep sample with gene expression data
pheno <- read.csv("./data/clinicalMatrix_BRCA", row.names = 1, sep = "\t")
pheno <- pheno[rownames(brca),]
  ## Find Her2 status: if IHC is not defined or equivocal, Her2 status is defined by FISH; if IHC is indeterminate and FISH is recorded, her2 status is defined by FISH
  ## Some samples has indeterminate IHC but empty FISH. These samples' Her2 status is not defined.
her2_stat <- pheno[,c("lab_proc_her2_neu_immunohistochemistry_receptor_status",
"lab_procedure_her2_neu_in_situ_hybrid_outcome_type")]
for (i in 1:1097){
  if (her2_stat[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"][i] == "" | her2_stat[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"][i] == "Equivocal"){
    her2_stat[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"][i] <- her2_stat[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"][i]
  }
  else if (her2_stat[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"][i] == "Indeterminate" & her2_stat[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"][i] != ""){
    her2_stat[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"][i] <- her2_stat[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"][i]
  }
}
subtypes <- data.frame(pheno[,"breast_carcinoma_estrogen_receptor_status"], pheno[,"breast_carcinoma_progesterone_receptor_status"], her2_stat)
colnames(subtypes) <- c("ER","PR","HER2","Hist")
for (i in 1:1097){
  subtypes[,"Hist"][i] <- ""
  if (sum(subtypes[i,] == c("Negative","Negative","Negative")) ==3) {subtypes[,"Hist"][i] <- "TNBC"}
  else if (sum(subtypes[i,] == c("Positive","Negative","Negative")) ==3) {subtypes[,"Hist"][i] <- "HR+"}
  else if (sum(subtypes[i,] == c("Negative","Positive","Negative")) ==3) {subtypes[,"Hist"][i] <- "HR+"}
  else if (sum(subtypes[i,] == c("Positive","Positive","Negative")) ==3) {subtypes[,"Hist"][i] <- "HR+"}
  else if (sum(subtypes[i,] == c("Negative","Negative","Positive")) ==3) {subtypes[,"Hist"][i] <- "HER2+"}
  else if (sum(subtypes[i,] == c("Negative","Positive","Positive")) ==3) {subtypes[,"Hist"][i] <- "HR+HER2+"}
  else if (sum(subtypes[i,] == c("Positive","Negative","Positive")) ==3) {subtypes[,"Hist"][i] <- "HR+HER2+"}
  else if (sum(subtypes[i,] == c("Positive","Positive","Positive")) ==3) {subtypes[,"Hist"][i] <- "HR+HER2+"}
}
  ## Manually fix individual cases
subtypes["TCGA-AN-A0FW-01","Hist"] <- "HR+"
subtypes["TCGA-D8-A1X8-01","Hist"] <- "HR+"
subtypes$pam50 <- pheno[rownames(subtypes),"PAM50Call_RNAseq"]

 ## Expression by subtypes
brca_hr <- brca[subtypes$Hist =="HR+",]
brca_her2 <- brca[subtypes$Hist == "HER2+",]
brca_tnbc <- brca[subtypes$Hist == "TNBC",]
  ## Define gene sets
ebtf <- c('MYC','ARNTL','HIF1A','MAX', 'CLOCK','ARNT')
eboxco <- c('SP1', 'SP2', 'YY1', 'NFYA', 'NFYB', 'E2F1', 'E2F4','CTCF', 'GABPA')
clock <- c('ARNTL','CLOCK','CRY1', 'CRY2', 'PER1', 'PER2', 'PER3', 'NR1D1', 'NR1D2')
clock_extended <- c('ARNTL','CLOCK','RORA','RORB','RORC','CRY1', 'CRY2', 'PER1', 'PER2', 'PER3', 'NR1D1', 'NR1D2')

 ## Utility line to check the distribution of single genes across the cohort
brcac %>%
  as.data.frame() %>%
  ggplot(aes(x=ARNTL)) +
  geom_density()

#2 Tumor-normal comparison and DEG analysis
 ## Creating a matrix that concatenates paired tumor and normal samples. 
 ## This will be used for conditional logistic model, thus named brca_clogit
brca_clogit <- brca[substring(rownames(brca),1,12) %in% substring(rownames(brca_norm),1,12),]
brca_clogit <- brca_clogit[order(row.names(brca_clogit)),]
brca_clogit <- brca_clogit[,c("ARNTL","CLOCK","CRY1","CRY2","PER1","PER2","PER3","NR1D1","NR1D2")]
brca_clogit <- as.data.frame(brca_clogit)
brca_clogit$id <- substring(rownames(brca_clogit),1,12)
brca_clogit$case <- rep(1,length(114))
brca_normclogit <-brca_norm[,clock]
brca_normclogit <- data.frame(brca_normclogit, id=substring(rownames(brca_clogit),1,12), case = rep(0,length(114)))
brca_clogit <- rbind(brca_clogit,brca_normclogit)

Tissue_type <- c(rep("tumor",114),rep("normal",114)) %>% factor()
names(Tissue_type) <- rownames(brca_clogit)
pheatmap(brca_clogit[,1:9], annotation_row = as.data.frame(Tissue_type), cluster_rows = F, show_rownames = F)
brca_norm_paired_tumor <- brca[substring(rownames(brca),1,12) %in% substring(rownames(brca_norm),1,12),]
plot_tumor_vs_normal("PER1")

## DEG testing using Wilcoxon (from recommendation by https://doi.org/10.1186/s13059-022-02648-4, 
## Pipeline adapted from https://rpubs.com/LiYumei/806213)
## I later tried edgeR as well. The results seems to be pretty consistent, with around 11,000 DEGs
## This is quite surprising given the total number of genes is only ~15,000, meaning most of the genes
## are actually DEG genes. This weakens the conclusion that "circadian genes expression are often
## dis-regulated in cancer. Most genes are, while the expression of circadian genes like ARNTL and CLOCK 
## are actually not significantly changed. This is very likely due to the ubiquity of circadian genes.
brca_deg <- data.frame(t(rbind(brca_norm_paired_tumor,brca_norm)))
brca_deg <- 2^brca_deg -1
brca_deg <- DGEList(brca_deg, group = Tissue_type)
keep <- filterByExpr(brca_deg)
brca_deg <- brca_deg[keep,keep.lib.sizes=FALSE]
genenames <- rownames(brca_deg)
brca_deg <- as.data.frame(brca_deg)
row.names(brca_deg)<- genenames
wilcoxon_p <- sapply(1:nrow(brca_deg),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(brca_deg[i,])),Tissue_type)
  p=wilcox.test(gene~Tissue_type, data)$p.value
  return(p)
})
fdr <- p.adjust(wilcoxon_p,method = "fdr")
dataCon1=brca_deg[,c(which(Tissue_type=="normal"))]
dataCon2=brca_deg[,c(which(Tissue_type=="tumor"))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
brca_wilcox <- data.frame(cbind(foldChanges,fdr))

## comparing normal tissue with all tumor tissues
conditions <- c(rep("tumor",1097),rep("normal",114)) %>% factor()
brca_deg <- data.frame(t(rbind(brca,brca_norm)))
brca_deg <- 2^brca_deg -1
brca_deg <- DGEList(brca_deg, group = conditions)
keep <- filterByExpr(brca_deg)
brca_deg <- brca_deg[keep,keep.lib.sizes=FALSE]
brca_deg <- as.data.frame(brca_deg)

wilcoxon_p <- sapply(1:nrow(brca_deg),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(brca_deg[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr <- p.adjust(wilcoxon_p,method = "fdr")
dataCon1=brca_deg[,c(which(conditions=="normal"))]
dataCon2=brca_deg[,c(which(conditions=="tumor"))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
brca_wilcox <- data.frame(cbind(foldChanges,fdr))

 ## DEG analysis using edgeR (according to edgeR manual, it can take RSEM estimated counts)
brca_edgeR <- data.frame(t(rbind(brca,brca_norm)))
brca_edgeR <- 2^brca_edgeR -1
brca_edgeR <- DGEList(brca_edgeR, group = conditions)
keep <- filterByExpr(brca_edgeR)
brca_edgeR <- brca_edgeR[keep,keep.lib.sizes=FALSE]
brca_edgeR <- estimateDisp(brca_edgeR)
brca_edgeR <- exactTest(brca_edgeR)
brca_edgeR[clock,]
 ## Plot bar plot of DE genes
clc_bar_plot <- brca_wilcox[clock,]
clc_bar_plot$gene <- rownames(clc_bar_plot)
clc_bar_plot$qsig <- c("*","n.s.","n.s.","****","****","****","****","****","****")
clc_bar_plot$sign <- as.factor(sign(clc_bar_plot$foldChanges))
ggbarplot(clc_bar_plot, "gene","foldChanges",
          label = clc_bar_plot$qsig, color = "sign",lab.size = 6,
          fill = "sign", palette = c("#00AFBB","#FC4E07"), 
          width = 0.8, xlab = "Circadian Genes", ylab = "log2FoldChange",
          title = "TCGA_BRCA_RNAseq Tumor v.s. Normal") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  guides(fill="none",color="none")

# Conditional logistic regression
clc_clogit <- clogit(case ~ ARNTL+CLOCK+strata(id),brca_clogit)
## Testing survival difference of conditional logistic regression model
brca_clscore <- clock_score_1(brca)
BC_Score <- brca_clscore[names(brca_clscore) %in% row.names(surv_brca)]
BC_Score <- ifelse(BC_Score > median(BC_Score), "high", "low")
fit_clst = survfit(Surv(OS.time,OS) ~ BC_Score, data = surv_brca)
print(fit_clst)
ggsurvplot(fit_clst, data = surv_brca, pval = T)
## by subtype
brca_clscore <- clock_score_5(brca_her2)
BC_Score <- brca_clscore[names(brca_clscore) %in% row.names(surv_brca)]
BC_Score <- ifelse(BC_Score > median(BC_Score), "high", "low")
fit_clst = survfit(Surv(OS.time,OS) ~ BC_Score, data = surv_brca[rownames(brca_her2),])
print(fit_clst)
ggsurvplot(fit_clst, data = surv_brca[rownames(brca_her2),], pval = T)


# Logistic regression
clc_logistic <- glm(case ~ ARNTL+CLOCK+CRY1+CRY2+PER1+PER2+PER3+NR1D1+NR1D2, family = binomial(link = 'logit'), data = brca_clogit)
clc_logit_verify <- predict(clc_logistic, as.data.frame(brca_clogit), type = "response")
clc_logit_verify <- ifelse(clc_logit_verify>0.5, 1,0)
clc_logit_train_err <- 1 - sum(clc_logit_verify == brca_clogit$case)/228
brca_wo_normal <- brca_clc[!rownames(brca_clc) %in% rownames(brca_clogit),]
clc_logit_test <- predict(clc_logistic, as.data.frame(brca_wo_normal), type = "response")
clc_logit_test <- ifelse(clc_logit_test >= 0.5, 1, 0)
clc_logit_test_err <- 1-sum(clc_logit_test)/982
# Testing ligistic model on METABRIC data set. This is not a valid test since METABRIC is from micro-array
#meta_clc <- metabric[,c("ARNTL","CLOCK","CRY1","CRY2","PER1","PER2","PER3","NR1D1","NR1D2")]
#meta_logit_test <- predict(clc_logistic, as.data.frame(meta_clc), type = "response")
#meta_logit_test <- ifelse(meta_logit_test >= 0.5, 1, 0)
#meta_logit_test_err <- 1-sum(meta_logit_test)/1980

brca_logistic <- clock_score_3(brca)
Logistic_score <- brca_logistic[names(brca_logistic) %in% rownames(surv_brca)]
Logistic_score <- ifelse(Logistic_score > median(Logistic_score), "high", "low")
fit_clst = survfit(Surv(OS.time,OS) ~ Logistic_score, data = surv_brca)
print(fit_clst)
ggsurvplot(fit_clst, data = surv_brca, pval = T, conf.int=TRUE,xlab="OS Time/Days")

subtypes_pam50 <- subtypes[rownames(subtypes) %in% rownames(surv_brca),]
surv_pam50 <- surv_brca[subtypes_pam50[,"pam50"]!="",]
surv_pam50 <- data.frame(surv_pam50, subtypes_pam50[subtypes_pam50[,"pam50"]!="",])
fit_clst = survfit(Surv(OS.time,OS) ~ surv_pam50$pam50, data = surv_pam50)
print(fit_clst)
ggsurvplot(fit_clst, data = surv_pam50, pval = T)

subtypes_hist <- subtypes[rownames(subtypes) %in% rownames(surv_brca),]
surv_hist <- surv_brca[subtypes_hist[,"Hist"]!="",]
surv_hist <- data.frame(surv_hist, subtypes_hist[subtypes_hist[,"Hist"]!="",])
fit_clst = survfit(Surv(OS.time,OS) ~ surv_hist$Hist, data = surv_hist)
print(fit_clst)
ggsurvplot(fit_clst, data = surv_hist, pval = T)


## Run same logistic regression analysis with z-score
brca_z_clogit <- as.data.frame(brca_z[rownames(brca_z) %in% rownames(brca_clogit),])
brca_z_clogit$case <- rep(1,length(114))
brca_norm_z_clogit <- as.data.frame(brca_norm_z)
brca_norm_z_clogit$case <- rep(2, length(114))
brca_z_clogit <- rbind(brca_z_clogit,brca_norm_z_clogit)
brca_z_clogit$id <- rep(1:114,2)
clogit(case ~ ARNTL+CLOCK+CRY1+CRY2+PER1+PER2+PER3+NR1D1+NR1D2+strata(id),brca_z_clogit)

brca_z_clscore <- clc_clogit_score_calc(brca_z)
cl_z_cluster <- as.data.frame(brca_z_clscore[rownames(cl_z_cluster) %in% rownames(surv_brca)])
cl_z_cluster$clst <- ifelse(cl_z_cluster[,1] > 4.55113, rep(1,1095), rep(2,1095))

fit_clst = survfit(Surv(OS.time,OS) ~ cl_z_cluster$clst, data = surv_brca)
print(fit_clst)
ggsurvplot(fit_clst, data = surv_brca, pval = T)


#3 Correlation analysis
  ## Plot correlation of any two genes in a scatter plot with regression line
scatter_x <- "ARNTL"
scatter_y <- "CLOCK"
ggscatter(as.data.frame(brca), x = scatter_x, y = scatter_y,
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = scatter_x, ylab = scatter_y)

  ## Correlation between EBTFs and CoTFs in tumor samples
brca_pearson <- as.data.frame(matrix(nrow = 6, ncol = 9))
rownames(brca_pearson) <- ebtf
colnames(brca_pearson) <- eboxco
brca_pearson_p <- brca_pearson

for (eb in ebtf){
  for (eco in eboxco){
    cor_res <- cor.test(brca[,eb], brca[,eco], method = 'pearson')
    brca_pearson[eb,eco] <- cor_res$estimate['cor']
    brca_pearson_p[eb,eco] <- cor_res$p.value
  }
}
  ## Correlation between Clock genes in tumor samples
clc_brca_pearson<- as.data.frame(matrix(nrow = 9, ncol = 9))
rownames(clc_brca_pearson) <- clock
colnames(clc_brca_pearson) <- clock
clc_brca_pearson_p <- clc_brca_pearson
for (eb in clock){
  for (eb2 in clock){
    clc_cor_res <- cor.test(brca_norm[,eb], brca_norm[,eb2], method = 'pearson')
    clc_brca_pearson[eb,eb2] <- clc_cor_res$estimate['cor']
    clc_brca_pearson_p[eb,eb2] <- clc_cor_res$p.value
  }
}
my_color = rev(paletteer_d('RColorBrewer::RdYlBu')[-1])
my_color = colorRampPalette(my_color)(10)
corrplot(as.matrix(clc_brca_pearson), 
         col = COL2('RdBu', 10),
         p.mat = as.matrix(clc_brca_pearson_p),
         sig.level = 0.01,
         insig = 'blank',
         addCoef.col = 'black',
         tl.col = 'black',
         tl.srt = 45,
         cl.pos = 'r')

  ## Plot heat map of clock genes to visually test if there are patterns
brca_clc <- brca_surv[,clock]
brca_norm_clc <- brca_norm[,clock]
heat <- pheatmap(t(brca_clc), cutree_cols = 3,
                 show_colnames = F)
normal_heat <- pheatmap(t(brca_norm_clc), show_colnames = F)
extended_clc <- brca_surv[,clock_extended]
extended_clc_heat <- pheatmap(t(extended_clc), show_colnames = F)
    ### Observed a B/C-low cluster
heat_cluster <- cutree(heat$tree_col,k=2)
clusplot(brca_clc, heat_cluster)

fit_clst = survfit(Surv(OS.time,OS) ~ heat_cluster, data = surv_brca)
print(fit_clst)
ggsurvplot(fit_clst, data = surv_brca, pval = T)
    ### Results: B/C low cluster shows better survival over time although not significant, 
    ###          perhaps because cross of survival curve in the beginning and small number 
    ###          of the first cluster.
    ###          This  has no clear relationship with histological subtypes of BRCA

  ## Run logistic regression on these two clusters to find parameters of clock gene expressions that determined this separation 
brca_heat_reg <- data.frame(brca_clc, heat_cluster)
   ### heat_cluster 2 is B/C low, assigned as 0 in logistic regression
brca_heat_reg[,"heat_cluster"] <- ifelse(brca_heat_reg[,"heat_cluster"]==2,FALSE,TRUE)
heat_reg <- glm(formula = heat_cluster ~ ARNTL+CLOCK+CRY1+CRY2+PER1+PER2+PER3+NR1D1+NR1D2, family = binomial(link = 'logit'), data = brca_heat_reg)
# heat_reg_plot <- heat_reg_calc(brca_clc)
# heat_reg_plot <- 1/(1+exp(-1 * heat_reg_plot))
heat_cluster_3 <- cutree(heat$tree_col,k=3)
clusplot(brca_clc,heat_cluster_3)
names(heat_cluster) == rownames(surv_brca)
## cluster 1 is the right cluster in heatmap and cluster 2 is the left cluster in heatmap
get_heat_cluster_1_3 <- heat_cluster_3[heat_cluster_3==1 | heat_cluster_3==2]
surv_heat_cluster_3 <- surv_brca[rownames(surv_brca) %in% names(get_heat_cluster_1_3),]
fit_clst = survfit(Surv(OS.time,OS) ~ get_heat_cluster_1_3, data = surv_heat_cluster_3)
print(fit_clst)
ggsurvplot(fit_clst, data = surv_heat_cluster_3, pval = T)

brca_clc_heat_cluster_1_3 <- as.data.frame(brca_clc[names(get_heat_cluster_1_3),])
brca_clc_heat_cluster_1_3$heat_cluster <- ifelse(get_heat_cluster_1_3==2,FALSE,TRUE)
heat_reg_3 <- glm(formula = heat_cluster ~ ARNTL+CLOCK+CRY1+CRY2+PER1+PER2+PER3+NR1D1+NR1D2, family = binomial(link = 'logit'), data = brca_clc_heat_cluster_1_3)
summary(heat_reg_3)


#4 Cox regression
os_brca <- surv_brca[,c('OS.time','OS')]
os_brca <- os_brca[!is.na(os_brca$OS.time),]
os_brca <- os_brca[os_brca$OS.time != 0,]
os_brca <- as.matrix(os_brca)
colnames(os_brca) <- c("time","status")
cox_clc <- as.matrix(brca_clc)
cox_clc <- cox_clc[row.names(brca_clc) %in% row.names(os_brca),]

coxreg <- glmnet(cox_clc, os_brca, family = 'cox')
plot(coxreg, xvar = "dev", label = TRUE)
coef(coxreg, s = 0.0)

survival::survfit(coxreg, s = 0.01, x = cox_clc, y = os_brca)
clc_os <- cbind(cox_clc, os_brca)
clc_os <- as.data.frame(clc_os)
multi_cox <- coxph(Surv(time,status) ~ ARNTL+CLOCK+CRY1+CRY2+PER1+PER2+PER3+NR1D1+NR1D2, data = clc_os)

summary(multi_cox)
plot(unicox)

cox_pheno <- pheno[,c("age_at_initial_pathologic_diagnosis","pathologic_stage","pathologic_T","pathologic_N","pathologic_M")]
cox_pheno$age_at_initial_pathologic_diagnosis <- cox_pheno$age_at_initial_pathologic_diagnosis > 60
cox_pheno$pathologic_stage <- cox_pheno$pathologic_stage == "Stage IV"
cox_pheno$pathologic_T <- substring(cox_pheno$pathologic_T, 1, 2) %in% c("T3","T4")
cox_pheno$pathologic_N <- substring(cox_pheno$pathologic_N, 1, 2) %in% c("N2","N3")
cox_pheno$pathologic_M <- substring(cox_pheno$pathologic_M, 1, 2) == "M1"
brca_coxscore <- clock_score_5(brca)
cox_pheno$clc_risk <- ifelse(brca_coxscore > median(brca_coxscore), "TRUE","FALSE")
cox_pheno$ER <- subtypes$ER == "Positive"
cox_pheno$PR <- subtypes$PR == "Positive"
cox_pheno$HER2 <- subtypes$HER2 == "Positive"
cox_pheno <- cox_pheno[rownames(surv_brca),]
cox_pheno$OS <- surv_brca$OS
cox_pheno$OS.time <- surv_brca$OS.time
multi_cox <- coxph(Surv(OS.time,OS) ~ ., data = cox_pheno)
summary(multi_cox)

univ_formulas <- sapply(clock,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = clc_os)})

brca_coxscore <- clock_score_5(brca_surv)
cox_cluster <- as.data.frame(brca_coxscore, row.names = rownames(brca_coxscore))
  ### median(brca_coxscore) = 0.9419143
cox_cluster <- ifelse(brca_coxscore > 0.9419143, 1, 0)

fit_clst = survfit(Surv(OS.time,OS) ~ cox_cluster, data = surv_brca)
print(fit_clst)
ggsurvplot( fit_clst, data = surv_brca, pval = T)






## 
clock_scores <- data.frame(
  clock_score_1(brca_clc),
  clock_score_2(brca_clc),
  clock_score_3(brca_clc),
  clock_score_4(brca_clc),
  clock_score_5(brca_clc))
clock_scores$hist <- subtypes[rownames(clock_scores),"Hist"]
clock_scores$pam50 <- subtypes[rownames(clock_scores),"pam50"]
colnames(clock_scores)[1:5] <- c("score_1","score_2","score_3","score_4","score_5")
clc_plots <- list(rep(0,5))
for (i in 1:5){
  over_median <- ifelse(clock_scores[,i]>= median(clock_scores[,i]),1,0)
  fit_clst = survfit(Surv(OS.time,OS) ~ over_median, data = surv_brca)
  print(fit_clst)
  clc_plots[i] <- ggsurvplot(fit_clst, data = surv_brca, pval = T)
}
clc_plots

clock_scores_z <- data.frame(
  clock_score_1(brca_z),
  clock_score_2(brca_z),
  clock_score_3(brca_z),
  clock_score_5(brca_z))
colnames(clock_scores_z)[1:4] <- c("score_1","score_2","score_3","score_5")
clock_scores_z$hist <- subtypes$Hist
clock_scores_z$pam50 <- clock_scores$pam50
clc_plots_z <- list(rep(0,4))
ggviolin(clock_scores, x = "hist", y = "score_5",
         color = "hist", pallette = "npg",
         add = c("jitter","mean_sd"), remove = "",
         label = "",
         xlab = "Subtype")+
  stat_compare_means(
    comparisons = list(c("TNBC", "HR+"),c("TNBC","HR+HER2+"), c("TNBC", "HER2+")),
    label = "p.signif"
  )


draft_df1 <- as.data.frame(brca_clc)
draft_df1$case <- rep("tumor",1095)
draft_df2 <- as.data.frame(brca_norm_clc)
draft_df2$case <- rep("normal",114)
all_sample_clc <- rbind(draft_df1,draft_df2)
clock_scores_all_sample <- data.frame(
  clock_score_1(all_sample_clc),
  clock_score_2(all_sample_clc),
  clock_score_3(all_sample_clc),
  clock_score_4(all_sample_clc))
colnames(clock_scores_all_sample) <- c("score_1","score_2","score_3","score_4")
clock_scores_all_sample$case <- all_sample_clc[,"case"]
ggviolin(clock_scores_all_sample, x = "case", y = "score_4",
         color = "case", pallette = "npg",
         add = c("jitter","mean_sd"), remove = "",
         label = "",
         xlab = "")

    ### Plot single genes by hist and pam50 subtype
brca_clc_plot<- data.frame(brca_clc,clock_scores$hist,clock_scores$pam50)
colnames(brca_clc_plot)[10:11] <- c("hist","pam50")
ggviolin(brca_clc_plot,x = "pam50",y = "ARNTL",
         color = "pam50", pallette = "npg",
         add = c("jitter","mean_sd"), remove = "",
         xlab = "Histological Subtype")

  ## KM_plot of single circadian genes in histo subtypes devided by median
  ## For reproduction of plots, simply replace the variable name of the subtype in the for loop,
  ## assign i back to 1 and rerun the for loop. This is not written in a function because
  ## it seems ggsurvplot does not take local variable of high_low. I might be wrong tho...
brca_clc_plot_hr <- brca_clc_plot[brca_clc_plot$hist == "HR+",]
brca_clc_plot_her2 <- brca_clc_plot[brca_clc_plot$hist == "HER2+",]
brca_clc_plot_tnbc <- brca_clc_plot[brca_clc_plot$hist == "TNBC",]
km_genes <- list(rep(0,9))


i <- 1
for (gene in clock){
  high_low <- ifelse(brca_clc_plot[,gene]>= median(brca_clc_plot[,gene]),"high","low")
  surv_data <- surv_brca[rownames(brca_clc_plot),]
  surv_fit = survfit(Surv(OS.time,OS) ~ high_low, data = surv_data)
  print(surv_fit)
  km_genes[i] <- ggsurvplot(surv_fit, data = surv_data, pval = T)
  i <- i + 1
}



##
brca_eb <- as.data.frame(brca[,colnames(brca) %in% ebtf])
clock_scores %>% as.data.frame() %>%
  pivot_longer(cols = everything()) %>%
  ggviolin(x = "name", y = "value")

