# Script for analyzing METABRIC BRCA data
# Created by Yuanzhong Pan on Jan 8th, 2024 from copied TCGA script of the same author
#
# Data Source: METABRIC Cohort: 
#      Gene Expression: Illumina microarray 
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

# Data preparation: read in data, clean data, and define gene sets
metabric <- read.csv("./data/data_mrna_illumina_microarray.txt", sep = "\t")
metabric <- metabric[order(metabric$Hugo_Symbol),]
metabric <- metabric[!duplicated(metabric$Hugo_Symbol),]
row.names(metabric) <- metabric$Hugo_Symbol
metabric <- metabric[,colnames(metabric)!= "Entrez_Gene_Id" & colnames(metabric)!= "Hugo_Symbol"]

metabric <- t(metabric)
  
  ## Fix row names (colnames cannot take "-" and they are replaced by "." automatically. After transpose, the names can be fixed to have dash)
rownames(metabric) <- gsub("\\.", "-", rownames(metabric))
metabric <- metabric[sort(rownames(metabric)),]
  ## Calculate z scores 
meta_z <- metabric
for (i in 1:length(colnames(meta_z))){
  meta_z[,i] <- (meta_z[,i] - mean(meta_z[,i]))/sd(meta_z[,i])
}

  ## Read patien data and only keep samples that contain expression data.
patient <- read.csv("./data/data_clinical_patient.txt", row.names = 1, sep = "\t", skip = 4)
patient <- patient[rownames(patient) %in% rownames(metabric),]
sample <- read.csv("./data/data_clinical_sample.txt", row.names = 1, sep = "\t", skip = 4)
sample <- sample[rownames(sample) %in% rownames(metabric),]
  ## All expression-available patients also have matched survival data.
mb_os <- patient[,c("OS_MONTHS","OS_STATUS")]
mb_os$OS_STATUS <- substr(as.character(mb_os$OS_STATUS), 0, 1)
mb_os$OS_STATUS <- as.numeric(mb_os$OS_STATUS)

mb_rfs <- patient[,c("RFS_MONTHS","RFS_STATUS")]
mb_rfs$RFS_STATUS <- substr(as.character(mb_rfs$RFS_STATUS), 0, 1)
mb_rfs$RFS_STATUS <- as.numeric(mb_rfs$RFS_STATUS)
  ## check that expression data and os data aligns
unique(rownames(metabric) == rownames(mb_os))

# BRCA Subtypes (Histological and PAM50)
subtypes <- data.frame(sample[,"ER_STATUS"], sample[,"PR_STATUS"], sample[,"HER2_STATUS"])
subtypes$Hist <- patient$HISTOLOGICAL_SUBTYPE
rownames(subtypes) <- rownames(sample)
for (i in 1:1980){
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
subtypes$pam50 <- patient$CLAUDIN_SUBTYPE

meta_HR <- metabric[subtypes_hist=="HR+",]
meta_HER2 <- metabric[subtypes_hist=="HER2+",]
meta_TNBC <- metabric[subtypes_hist == "TNBC",]

  ## Define gene sets
ebtf <- c('MYC','ARNTL','HIF1A','MAX', 'CLOCK','ARNT')
eboxco <- c('SP1', 'SP2', 'YY1', 'NFYA', 'NFYB', 'E2F1', 'E2F4','CTCF', 'GABPA')
clock <- c('ARNTL','CLOCK','CRY1', 'CRY2', 'PER1', 'PER2', 'PER3', 'NR1D1', 'NR1D2')
clock_extended <- c('ARNTL','CLOCK','RORA','RORB','RORC','CRY1', 'CRY2', 'PER1', 'PER2', 'PER3', 'NR1D1', 'NR1D2')


# Correlation analysis
  ## Plot correlation of any two genes in a scatter plot with regression line
scatter_x <- "SP1"
scatter_y <- "CLOCK"
ggscatter(as.data.frame(metabric), x = scatter_x, y = scatter_y,
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = scatter_x, ylab = scatter_y)

  ## Correlation between EBTFs and CoTFs in tumor samples
meta_pearson <- as.data.frame(matrix(nrow = 6, ncol = 9))
rownames(meta_pearson) <- ebtf
colnames(meta_pearson) <- eboxco
meta_pearson_p <- meta_pearson

for (eb in ebtf){
  for (eco in eboxco){
    cor_res <- cor.test(metabric[,eb], metabric[,eco], method = 'pearson')
    meta_pearson[eb,eco] <- cor_res$estimate['cor']
    meta_pearson_p[eb,eco] <- cor_res$p.value
  }
}
  ## Correlation between Clock genes in tumor samples
clc_meta_pearson<- as.data.frame(matrix(nrow = 9, ncol = 9))
rownames(clc_meta_pearson) <- clock
colnames(clc_meta_pearson) <- clock
clc_meta_pearson_p <- clc_meta_pearson
for (eb in clock){
  for (eb2 in clock){
    clc_cor_res <- cor.test(metabric[,eb], metabric[,eb2], method = 'pearson')
    clc_meta_pearson[eb,eb2] <- clc_cor_res$estimate['cor']
    clc_meta_pearson_p[eb,eb2] <- clc_cor_res$p.value
  }
}
my_color = rev(paletteer_d('RColorBrewer::RdYlBu')[-1])
my_color = colorRampPalette(my_color)(10)
corrplot(as.matrix(meta_pearson), 
         col = COL2('RdBu', 10),
         p.mat = as.matrix(meta_pearson_p),
         sig.level = 0.01,
         insig = 'blank',
         addCoef.col = 'black',
         tl.col = 'black',
         tl.srt = 45,
         cl.pos = 'r')

  ## Plot heat map of clock genes to visually test if there are patterns
meta_clc <- metabric[,clock]
heat <- pheatmap(t(meta_clc), cutree_cols = 2,
                 show_colnames = F)

extended_clc <- metabric[,clock_extended]
extended_clc_heat <- pheatmap(t(extended_clc), show_colnames = F)
heat_cluster <- cutree(heat$tree_col,k=2)
clusplot(meta_clc, heat_cluster)
      ### Observed a PER2-CRY2 high part, but overall no obvious clustering

fit_clst = survfit(Surv(OS_MONTHS,OS_STATUS) ~ heat_cluster, data = mb_os)
print(fit_clst)
ggsurvplot(fit_clst, data = mb_os, pval = T)
      ### Results: PER2-CRY2 cluster has significant indication on survival 

# Cox regression
clc_os <- cbind(metabric[,clock], mb_os)
multi_cox <- coxph(Surv(OS_MONTHS,OS_STATUS) ~ ARNTL+CLOCK+CRY1+CRY2+PER1+PER2+PER3+NR1D1+NR1D2, data = clc_os)
summary(multi_cox)

meta_coxscore <- clock_score_3(metabric)
cox_cluster <- as.data.frame(meta_coxscore, row.names = rownames(meta_coxscore))
cox_cluster <- ifelse(meta_coxscore > median(meta_coxscore), "High", "Low")
fit_clst = survfit(Surv(OS_MONTHS,OS_STATUS) ~ cox_cluster, data = mb_os)
print(fit_clst)
ggsurvplot( fit_clst, data = mb_os, pval = T, conf.int = FALSE, xlab="OS Time/Months")

# Conditional logistic model from TCGA using only BMAL1 and CLOCK
brca_clscore <- clock_score_1(metabric)
BC_Score <- ifelse(brca_clscore > median(brca_clscore), "high", "low")
fit_clst = survfit(Surv(OS_MONTHS,OS_STATUS) ~ BC_Score, data = mb_os)
print(fit_clst)
ggsurvplot(fit_clst, data = mb_os, pval = T)


subtypes_pam50 <- subtypes$pam50
subtypes_hist <- subtypes$Hist
fit_clst = survfit(Surv(OS_MONTHS,OS_STATUS) ~ subtypes_hist, data = mb_os)
print(fit_clst)
ggsurvplot(fit_clst, data = mb_os, pval = T)




## 
clock_scores <- data.frame(
  clock_score_1(metabric),
  clock_score_5(metabric))
clock_scores$hist <- subtypes[rownames(clock_scores),"Hist"]
clock_scores$pam50 <- subtypes[rownames(clock_scores),"pam50"]
colnames(clock_scores)[1:2] <- c("score_1","score_5")
clc_plots <- list(rep(0,2))
for (i in 1:2){
  over_median <- ifelse(clock_scores[,i]>= median(clock_scores[,i]),1,0)
  fit_clst = survfit(Surv(OS_MONTHS,OS_STATUS) ~ over_median, data = mb_os)
  print(fit_clst)
  clc_plots[i] <- ggsurvplot(fit_clst, data = mb_os, pval = T)
}
clc_plots

clock_scores_z <- data.frame(
  clock_score_1(meta_z),
  clock_score_5(meta_z))
colnames(clock_scores_z)[1:2] <- c("score_1","score_5")
clock_scores_z$hist <- clock_scores$hist
#clock_scores_z$pam50 <- clock_scores$pam50
clc_plots_z <- list(rep(0,2))
ggviolin(clock_scores_z, x = "hist", y = "score_1",
         color = "hist", pallette = "npg",
         add = c("jitter","mean_sd"), remove = "",
         label = "",
         xlab = "Subtype") +
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
ggviolin(clock_scores_all_sample, x = "case", y = "score_1",
         color = "case", pallette = "npg",
         add = c("jitter","mean_sd"), remove = "",
         label = "",
         xlab = "")

    ### Plot single genes by hist and pam50 subtype
meta_clc_plot <- data.frame(brca_clc,clock_scores$hist,clock_scores$pam50)
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
##
meta_coxscore %>%
  as.data.frame() %>%
  ggplot(aes(x=cox_score)) +
  geom_density()
