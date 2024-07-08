## Functions for tcga.R script

calc_clock_score <- function(mat, model, mode){
  
}

clock_score_1 <- function(mat){
  clock_score <- mat[,"ARNTL"]+mat[,"CLOCK"]-(mat[,"CRY1"]+mat[,"CRY2"]+mat[,"PER1"]+mat[,"PER2"]+mat[,"PER3"]+mat[,"NR1D1"]+mat[,"NR1D2"])
  return(clock_score) 
}
clock_score_2 <- function(mat){
  clock_score_2 <- 0.6466*0*mat[,"ARNTL"]+2.3007*mat[,"CLOCK"]+0.7781*0*mat[,"CRY1"]+1.5367*mat[,"CRY2"]-
  0.3514*0*mat[,"PER1"]-0*mat[,"PER2"]+1.2147*mat[,"PER3"]+0.6358*mat[,"NR1D1"]+1.5115*mat[,"NR1D2"]
  return(clock_score_2)
}

# Coefficients from logistic regression model on all core clock genes in TCGA
clock_score_3 <- function(mat){
  clock_score_3 <- -1.4099*mat[,"ARNTL"] + 1.9499 * mat[,"CLOCK"]+0.3856*0*mat[,"CRY1"]-3.6606*mat[,"CRY2"]-1.7047*mat[,"PER1"] -
  1.2683*0*mat[,"PER2"]+1.6071*0*mat[,"PER3"]+1.2580*0*mat[,"NR1D1"]-3.5315*mat[,"NR1D2"]
  return(clock_score_3)
}
# Coefficients from conditional logistic model ran on only BMAL1 and CLOCK
clock_score_4 <- function(mat){
  clock_score_4 <- 0.6591 * mat[,"CLOCK"] -0.8459 * mat[,"ARNTL"]
  return(clock_score_4)
}

#multi-variate cox from TCGA data
clock_score_5 <- function(mat){
  clock_score_5 <- 0.254713 * mat[,"CLOCK"] - 0.336401 * mat[,"PER2"] + 0.223933 * mat[,"NR1D1"]
  return(clock_score_5)
}

# LASSO penalized multi-variate cox when lambda = 0.01
clock_score_6 <- function(mat){
  clock_score <- 0.138 * mat[,"CLOCK"] - 0.164 * mat[,"PER2"] + 0.06 * mat[,"NR1D1"] 
  return(clock_score)
}
plot_tumor_vs_normal <- function(gene_name){
  tumor_vs_norm <- data.frame(brca_norm[,gene_name],brca_norm_paired_tumor[,gene_name])
  colnames(tumor_vs_norm) <- c("normal","tumor")
  ggscatter(tumor_vs_norm, x = "normal", y = "tumor",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "normal", ylab = "tumor", title = gene_name)
}


