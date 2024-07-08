 ## Created by Yuanzhong Pan on Feb 13, 2024. For analyzing STARR-seq data.

library("ggplot2")
library("plotly")
library("tidyverse")
library("uwot")
library("Rtsne")
library("magrittr")
library("ggpubr")
library('pheatmap')
library('paletteer')
library('corrplot')
library('RColorBrewer')
library('GenomicRanges')
library('Biostrings')
library('BSgenome.Hsapiens.UCSC.hg38')
library('universalmotif')
library("TFBSTools")
library('factoextra')

# Prepare TFBS motifs and naming vectors
motif_ebassoc <- file.path("./motifs/EBTFassoc.meme")
motif_jaspar <- file.path("./motifs/JASPAR_2022_matrix_clustering_vertebrates_UNVALIDATED_cluster_root_motifs.tf")
ebassoc <- read_meme(motif_ebassoc)
jaspar_nonred <- read_transfac(motif_jaspar)
  
ebtf_names <- c('EBOX', 'SP', 'E2F', 'NFY', 'ETS', 'IRF', 'SIX5', 'PBX3', 'YY1', 'RUNX', 'JUN')
treatments <- c("DMSO", "SHP1705", "MG132", "COMBO")
var_names <- c("dmso", "shp", "mg", "combo")


# Data preparation and initial TFBS counting 
for (i in 1:4){
  print(paste("processiong group",var_names[i],sep = " "))
  # read in files and get peak sequences
  peak_name <- file.path(paste(treatments[i],"peaks.xls",sep = "_"))
  starr_peaks <- read.csv(peak_name, sep = "\t", skip = 21, header = TRUE)
  starr_peaks <- GRanges(starr_peaks)
  seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38,starr_peaks)
  starr_peaks$seq <- seq
  assign(paste("peak",var_names[i],sep = "_"), starr_peaks)
    ## count occurrence of EBassocs 
  ebcounts <- tf_mat(seq,ebassoc, match.score = "80%")
  colnames(ebcounts) <- ebtf_names
  assign(paste("ebcounts",var_names[i], sep = "_"),ebcounts)
  count_sum <- colSums(ebcounts)
  assign(paste("ebcounts_sum",var_names[i], sep = "_"),count_sum)
    ## count all non-redundant TFs with threshold 80%
  tfcounts <- tf_mat(seq,jaspar_nonred, match.score = "80%")
  assign(paste("tfcounts",var_names[i], sep = "_"),tfcounts)
    ## save image every cycle to prevent data loss
  save.image('starr.RData')
}
################################
# Plotting heat map of counting all TFs. This part was not included in the loop because the pheatmap function 
# takes long time to calculate row clusters. Data is log2-transformed to prevent color inflation due to high-count motifs.
# pheatmap function also only takes a maximum rows of 65535, for shp and mg group row 1-65535 are maintained to plot heatmap.

tfcounts_log <- log2(tfcounts_mg+1)
tfcounts_log <- tfcounts_log[1:65535,]
tfheat <- pheatmap(tfcounts_log, cluster_cols = FALSE)
tfheat_mg <- tfheat
################################

# Utility line for checking if there are motifs that does not appear in any peaks, 
# thus causing an all-zero vector. This might cause problems in later analysis.
# count.empty.cluster(tfcounts)

# Summarize counts of EBassoc TFs in a matrix
ebcounts_sum <- rbind(ebcounts_sum_dmso,ebcounts_sum_shp,ebcounts_sum_mg,ebcounts_sum_combo)

################################
# Calculate normalized motif representation in 
ebsum_dmso <- colSums(ebcounts_dmso*peak_dmso$pileup)/dim(ebcounts_dmso)[1]/11.37
ebsum_shp <- colSums(ebcounts_shp*peak_shp$pileup)/dim(ebcounts_shp)[1]/34.99
ebsum_mg <- colSums(ebcounts_mg*peak_mg$pileup)/dim(ebcounts_mg)[1]/11.18
ebsum_combo <- colSums(ebcounts_combo*peak_combo$pileup)/dim(ebcounts_combo)[1]/12.04
ebsum <- rbind(ebsum_dmso,ebsum_shp,ebsum_mg,ebsum_combo)
ebsum_norm <- ebsum
for (i in 1:4) { ebsum_norm[i,] <- ebsum[i,]/ebsum[1,]}
################################


# EBASSOC co-occurrence analysis
for (i in 3:4){
  print(i)
  tfmat <- eval(as.symbol(paste("ebcounts", var_names[i], sep = "_")))
  tflist <- ebassoc
  tf.coocur <- co.occur(tfmat, tflist)
  assign(paste("ebco", var_names[i], sep = "_"), tf.coocur)
  print("assign")
  ## Bonferroni correction: for 11 EBTFs p should be 0.001818 (alpha 0.05) or 0.0003636(alpha 0.01)
  tf.co.sig <- matrix(as.numeric(tf.coocur < 0.001818), ncol = 11)
  assign(paste("tfco_sig", var_names[i], sep = "_"), tf.co.sig)
  colnames(tf.co.sig)<- colnames(tf.coocur)
  rownames(tf.co.sig)<- colnames(tf.coocur)
  print("plotting heatmap")
  co_heat <- pheatmap(tf.co.sig, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, legend = FALSE, main = treatments[i])
  assign(paste("ebco_heat", var_names[i], sep = "_"), co_heat)
}

pheatmap(1-ebco_shp, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, legend = FALSE,main = "SHP1705")


# All Jaspar TFBS co-occurrence analysis
for (i in 1:4){
  print(i)
  tfmat <- eval(as.symbol(paste("tfcounts", var_names[i], sep = "_")))
  tflist <- jaspar_nonred
  tf.coocur <- co.occur(tfmat, tflist)
  print("cooccur analysis done")
  assign(paste("tfco", var_names[i], sep = "_"), tf.coocur)
  ## Bonferroni correction: for 258-23=225 non-redundant TFs, p should be 1.98e-6 (significant level 0.05)
  tf.co.sig <- matrix(as.numeric(tf.coocur < 0.00000198), ncol = 258)
  assign(paste("tfco_sig", var_names[i], sep = "_"), tf.co.sig)
  colnames(tf.co.sig)<- colnames(tf.coocur)
  rownames(tf.co.sig)<- colnames(tf.coocur)
  co_heat <- pheatmap(tf.co.sig, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 20, legend = FALSE,show_rownames = FALSE,show_colnames = FALSE)
  assign(paste("tfco_heat", var_names[i], sep = "_"), co_heat)
}



## Analyze existence of motifs (binary: is or not)
for (i in var_names){
  print(i)
  motif_bi <- eval(as.symbol(paste("ebcounts", i, sep = "_")))
  motif_bi <- motif_bi>0
  motif_bi <- matrix(as.numeric(motif_bi), ncol = 11)
  assign(paste("binary", i, sep = "_"), motif_bi)
  motif_bi <- distinct(as.data.frame(motif_bi))
  ranks <- nrow(motif_bi)
  rownames(motif_bi) <- paste(rep(i, times = ranks), seq(ranks), sep = "_")
  colnames(motif_bi) <- ebtf_names
  assign(paste("uniq_binary", i, sep = "_"), motif_bi)
}
all_binary <- rbind(uniq_binary_dmso, uniq_binary_shp, uniq_binary_mg, uniq_binary_combo)
 # ebtypes lists all ebassoc types 
ebtypes <- all_binary[!duplicated(all_binary),]
 ### merge function will take the intersection of dmso and combo types
dmso_combo <- merge(uniq_binary_combo,uniq_binary_dmso)
 ### Construct a vector in which the positions corresponding to dmso types will indicate if the type is also observed in combo.
 ### If not, it is a potentially repressed type.
logi_vec <-  duplicated(rbind(dmso_combo, uniq_binary_dmso))
combo_suppressed <- uniq_binary_dmso[!logi_vec[584:1278],]

# Calculate normalized binary motif representation
ebsum_bi_dmso <- colSums(binary_dmso*peak_dmso$pileup)/dim(ebcounts_dmso)[1]/11.37
ebsum_bi_shp <- colSums(binary_shp*peak_shp$pileup)/dim(ebcounts_shp)[1]/34.99
ebsum_bi_mg <- colSums(binary_mg*peak_mg$pileup)/dim(ebcounts_mg)[1]/11.18
ebsum_bi_combo <- colSums(binary_combo*peak_combo$pileup)/dim(ebcounts_combo)[1]/12.04
ebsum_bi <- rbind(ebsum_bi_dmso,ebsum_bi_shp,ebsum_bi_mg,ebsum_bi_combo)
ebsum_bi_norm <- ebsum_bi
for (i in 1:4) { ebsum_bi_norm[i,] <- ebsum_bi[i,]/ebsum_bi[1,]}



ebtypes_dmso <- score_type(ebcounts_dmso, ebtypes,peak_dmso)
ebtypes_shp <- score_type(ebcounts_shp, ebtypes,peak_shp)
ebtypes_mg <- score_type(ebcounts_mg, ebtypes,peak_mg)
ebtypes_combo <- score_type(ebcounts_combo, ebtypes,peak_combo)

ebtypes_summary <- cbind(ebtypes_dmso,ebtypes_combo,ebtypes_shp,ebtypes_mg)
combo_supp_type <- ebtypes[ebtypes_summary[,1]/ebtypes_summary[,2] >2,]
combo_supp_type <- combo_supp_type[!is.na(combo_supp_type[,1]),]
pheatmap(combo_supp_type, cluster_cols = FALSE, color = c("#FFFFFF","#99FF99"), fontsize_col = 20, show_rownames = FALSE, legend = FALSE)

combo_sup_peaks <- match_type(ebcounts_dmso, combo_supp_type)
combo_sup_peaks_seq <- peak_dmso[combo_sup_peaks,]$seq
names(combo_sub_peaks_seq) <- peak_dmso[combo_sup_peaks,]$name
writeXStringSet(combo_sub_peaks_seq, "combo_sup.fa")
 # Unility for plotting eb-types of any group in a heatmap
pheatmap(ebtypes, cluster_cols = FALSE, color = c("#FFFFFF","#99FF99"), fontsize_col = 20, show_rownames = FALSE, legend = FALSE)

is_combo_supp <- match_type(uniq_binary_dmso,combo_supp_type)
supp_type_logit <- cbind(uniq_binary_dmso,is_combo_supp)
sup_logistic <- glm(is_combo_supp ~ ., family = binomial(link = 'logit'), data = uniq_binary_dmso)
summary(sup_logistic)

supp_co_type <- as.data.frame(cooccur_vec_from_binary(uniq_binary_dmso))
supp_co_logistic <- glm(is_combo_supp ~ ., family = binomial(link = 'logit'), data = supp_co_type)
summary(supp_co_logistic)
accuracy <- sum((predict(supp_co_logistic, supp_co_type) > -1.682962) & is_combo_supp)/188
write.table(sort(supp_co_logistic$coefficients, decreasing = TRUE), "supp_co_logistic_coefficients.txt", quote = FALSE, sep = "\t")

# Plot PCA
pc <- uniq_binary_dmso
pc <- pc[,which(apply(pc, 2, var)!=0)]

pc_pca <- prcomp(pc, scale. = TRUE)
pc_eigenvalues <- pc_pca$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)),
                         variance = pc_eigenvalues) %>%
                  mutate(pct = variance/sum(variance)*100) %>%
                  mutate(pct_cum = cumsum(pct))
pc_eigenvalues %>% 
  ggplot(aes (x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) +
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal Componant", y = "Fraction variance explained")

pc_pc_scores <- pc_pca$x %>%
  as_tibble(rownames = "sample")
pc_pc_scores %>%
  ggplot(aes(x = PC1, y = PC2, color = is_combo_supp, label = NA)) +
  geom_point() +
  geom_text()

pccluster_shp <- eclust(pc_pca$x, "kmeans", hc_metric = "euclidean", k=2)

fviz_pca_ind(pc_pca, col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
)

sub_pc <- matrix(as.numeric(sub_pc), ncol = 11)
motif_umap <- umap(pc)
plot(
  motif_umap,
  cex = 0.1,
  main = "R uwot::umap",
  xlab = "",
  ylab = ""
)
