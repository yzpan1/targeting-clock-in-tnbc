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
library('motifStack')
library('universalmotif')
library("TFBSTools")
library('factoextra')
library("dplyr")
library("tidyr")
library('EnhancedVolcano')
source('utils.r')

# Prepare TFBS motifs and naming vectors
motif_ebassoc <- file.path("./motif_database/EBTFassoc.meme")
motif_jaspar <- file.path("./motif_database/JASPAR_2022_matrix_clustering_vertebrates_UNVALIDATED_archive/interactive_trees/JASPAR_2022_matrix_clustering_vertebrates_UNVALIDATED_cluster_root_motifs.tf")
ebassoc <- read_meme(motif_ebassoc)
jaspar_nonred <- read_transfac(motif_jaspar)

  ## Used for testing different definition of PFMs
  #motif_ebassoc <- file.path("./motif_database/EBTFassoc_jas_core_E2F.txt")
  #ebassoc <- read_transfac(motif_ebassoc)

###########################
# Marking if the motif is palindromic. After some discussion with computational biologists I think this section is unnecessary, since
# DNA binding is direction-selective and palindromic motifs can be considered 2 motifs since it can bind in two directions and have double
# the chance to bind to DNA in a random environment, even though the piece of DNA that is bound is the same. 
  ## EBOX and E2F motifs are palindromic
#  ebassoc[[1]]['extrainfo'] <- "1"
#  ebassoc[[3]]['extrainfo'] <- "1"
  ## Assign if the motif is palindromic, the 1/0 list was hand curated according to the JASPAR 2022 matrix cluster results
#  is.palin <- read.csv("JASPAR_2022_matrix_clustering_vertebrates_UNVALIDATED_cluster_if_palindromic.txt",sep = " ", row.names = 1, header = FALSE)
#  for (i in 1:length(jaspar_nonred)){jaspar_nonred[[i]]['extrainfo'] <- as.character(is.palin[jaspar_nonred[[i]]['name'],2])}
############################
  
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
for (i in 1:4){
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
# Create a dataframe containing all possible combinations of 11 motifs. 
all_types <- expand.grid(replicate(11, c(0,1), simplify = FALSE))
rownames(all_types) <- paste0("type_", seq_len(nrow(all_types)))
colnames(all_types) <- colnames(ebtypes)

# Create quantification dataframe containing each peaks type and pile up
peak_dmso$type <- annotate_type(ebcounts_dmso,all_types)
quant_dmso <- mcols(peak_dmso)[,c("name","pileup","type")] %>% as.data.frame()
quant_dmso <- quant_dmso[order(quant_dmso$type),]

peak_shp$type <- annotate_type(ebcounts_shp,all_types)
quant_shp <- mcols(peak_shp)[,c("name","pileup","type")] %>% as.data.frame()
quant_shp <- quant_shp[order(quant_shp$type),]

peak_mg$type <- annotate_type(ebcounts_mg,all_types)
quant_mg <- mcols(peak_mg)[,c("name","pileup","type")] %>% as.data.frame()
quant_mg <- quant_mg[order(quant_mg$type),]

peak_combo$type <- annotate_type(ebcounts_combo,all_types)
quant_combo <- mcols(peak_combo)[,c("name","pileup","type")] %>% as.data.frame()
quant_combo <- quant_combo[order(quant_combo$type),]

quant_dmso$tpm <- quant_dmso$pileup/sum(quant_dmso$pileup) * 100000
quant_shp$tpm <- quant_shp$pileup/sum(quant_shp$pileup) * 100000
quant_mg$tpm <- quant_mg$pileup/sum(quant_mg$pileup) * 100000
quant_combo$tpm <- quant_combo$pileup/sum(quant_combo$pileup) * 100000

#
#
#
#




hist(quant_shp$pileup,)
ggscatter(quant_shp, x="type",y="pileup", 
          #ylim = c(0,2000), 
          #xlim = c(3, 5)
)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# check the total counts for each group
sum(quant_dmso$pileup)
sum(quant_shp$pileup)
sum(quant_mg$pileup)
sum(quant_combo$pileup)

#type_counts <- as.data.frame(table(quant_combo$type))
#colnames(type_counts) <- c("type","counts")
#ggbarplot(type_counts, x = "type", y = "counts", ylim = c(0,8000))+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot scatter to visualize the relationship between peak counts and read counts
peak_count_pileup_dmso <- aggregate(pileup ~ type, data = quant_dmso, sum)
peak_count_pileup_dmso$peak_count <- table(quant_dmso$type)
ggscatter(peak_count_pileup_dmso, x = "peak_count", y = "pileup")
ggscatter(peak_count_pileup_combo, x = "peak_count", y = "pileup",
          xlim = c(0,1000), 
          ylim = c(0,250000),
          #label = "type",
          #label.select = peak_count_pileup_dmso$type[peak_count_pileup_dmso$peak_count > 1400],
          repel = TRUE)+
  xlab("Peak Counts") +
  ylab("Read Counts")

peak_count_pileup_shp <- aggregate(pileup ~ type, data = quant_shp, sum)
peak_count_pileup_shp$peak_count <- table(quant_shp$type)
ggscatter(peak_count_pileup_shp, x = "peak_count", y = "pileup")
ggscatter(peak_count_pileup_shp, x = "peak_count", y = "pileup",
          xlim = c(0,7000), 
          ylim = c(0,3e6),
          label = "type",
          label.select = peak_count_pileup_dmso$type[peak_count_pileup_dmso$peak_count > 1400],
          repel = TRUE)+
  xlab("Peak Counts") +
  ylab("Read Counts")

peak_count_pileup_mg <- aggregate(pileup ~ type, data = quant_mg, sum)
peak_count_pileup_mg$peak_count <- table(quant_mg$type)
ggscatter(peak_count_pileup_mg, x = "peak_count", y = "pileup")
ggscatter(peak_count_pileup_mg, x = "peak_count", y = "pileup",
          xlim = c(0,7000), 
          ylim = c(0,3e6),
          label = "type",
          label.select = peak_count_pileup_dmso$type[peak_count_pileup_dmso$peak_count > 1400],
          repel = TRUE)+
  xlab("Peak Counts") +
  ylab("Read Counts")

peak_count_pileup_combo <- aggregate(pileup ~ type, data = quant_combo, sum)
peak_count_pileup_combo$peak_count <- table(quant_combo$type)
ggscatter(peak_count_pileup_combo, x = "peak_count", y = "pileup",
          xlim = c(0,7000), 
          ylim = c(0,3e6),
          label = "type",
          label.select = peak_count_pileup_combo$type[peak_count_pileup_combo$pileup > 500000],
          repel = TRUE,
          ) +
  xlab("Peak Counts") +
  ylab("Read Counts")




tpm_dmso <- aggregate(tpm ~ type, data = quant_dmso, mean)
tpm_dmso$var <- aggregate(tpm ~ type, data = quant_dmso, var)$tpm
tpm_dmso$rep <- peak_count_pileup_dmso$peak_count
tpm_dmso$logtpm <- log2(tpm_dmso$tpm)
tpm_dmso$logvar <- log2(tpm_dmso$var)

tpm_shp <- aggregate(tpm ~ type, data = quant_shp, mean)
tpm_shp$var <- aggregate(tpm ~ type, data = quant_shp, var)$tpm
tpm_shp$rep <- peak_count_pileup_shp$peak_count
tpm_shp$logtpm <- log2(tpm_shp$tpm)
tpm_shp$logvar <- log2(tpm_shp$var)

tpm_mg <- aggregate(tpm ~ type, data = quant_mg, mean)
tpm_mg$var <- aggregate(tpm ~ type, data = quant_mg, var)$tpm
tpm_mg$rep <- peak_count_pileup_mg$peak_count
tpm_mg$logtpm <- log2(tpm_mg$tpm)
tpm_mg$logvar <- log2(tpm_mg$var)

tpm_combo <- aggregate(tpm ~ type, data = quant_combo, mean)
tpm_combo$var <- aggregate(tpm ~ type, data = quant_combo, var)$tpm
tpm_combo$rep <- peak_count_pileup_combo$peak_count
tpm_combo$logtpm <- log2(tpm_combo$tpm)
tpm_combo$logvar <- log2(tpm_combo$var)
ggscatter(tpm_mg, x = "logtpm", y = "logvar",
          #xlim = c(0,9), 
          #ylim = c(0,5)
          )
hist(tpm_dmso$var)




type_summary <- full_join(tpm_dmso, tpm_combo, by = "type", suffix = c("_dmso", "_combo"))
type_summary[is.na(type_summary)] <- 0
type_summary <- type_summary[type_summary$rep_dmso > 1 & type_summary$rep_combo > 1,]

type_summary_shp <- full_join(tpm_dmso, tpm_shp, by = "type", suffix = c("_dmso", "_shp"))
type_summary_shp[is.na(type_summary_shp)] <- 0
type_summary_shp <- type_summary_shp[type_summary_shp$rep_dmso > 1 & type_summary_shp$rep_shp > 1,]

type_summary_mg <- full_join(tpm_dmso, tpm_mg, by = "type", suffix = c("_dmso", "_mg"))
type_summary_mg[is.na(type_summary_mg)] <- 0
type_summary_mg <- type_summary_mg[type_summary_mg$rep_dmso > 1 & type_summary_mg$rep_mg > 1,]

de_res <- data.frame(matrix(0, nrow = dim(type_summary)[1], ncol = 2))
rownames(de_res) <- type_summary$type
colnames(de_res) <- c("foldChange", "pvalue")
de_res$foldChange <- type_summary$tpm_combo/type_summary$tpm_dmso
de_res$pvalue <- mapply(
  welch_p,
  type_summary$logtpm_dmso, type_summary$logvar_dmso, type_summary$rep_dmso,
  type_summary$logtpm_combo, type_summary$logvar_combo, type_summary$rep_combo
)
de_res$logFoldChange <- log2(de_res$foldChange)

de_res_shp <- data.frame(matrix(0, nrow = dim(type_summary_shp)[1], ncol = 2))
rownames(de_res_shp) <- type_summary_shp$type
colnames(de_res_shp) <- c("foldChange", "pvalue")
de_res_shp$foldChange <- type_summary_shp$tpm_shp/type_summary_shp$tpm_dmso
de_res_shp$pvalue <- mapply(
  welch_p,
  type_summary_shp$logtpm_dmso, type_summary_shp$logvar_dmso, type_summary_shp$rep_dmso,
  type_summary_shp$logtpm_shp, type_summary_shp$logvar_shp, type_summary_shp$rep_shp
)
de_res_shp$logFoldChange <- log2(de_res_shp$foldChange)

de_res_mg <- data.frame(matrix(0, nrow = dim(type_summary_mg)[1], ncol = 2))
rownames(de_res_mg) <- type_summary_mg$type
colnames(de_res_mg) <- c("foldChange", "pvalue")
de_res_mg$foldChange <- type_summary_mg$tpm_mg/type_summary_mg$tpm_dmso
de_res_mg$pvalue <- mapply(
  welch_p,
  type_summary_mg$logtpm_dmso, type_summary_mg$logvar_dmso, type_summary_mg$rep_dmso,
  type_summary_mg$logtpm_mg, type_summary_mg$logvar_mg, type_summary_mg$rep_mg
)
de_res_mg$logFoldChange <- log2(de_res_mg$foldChange)

EnhancedVolcano(de_res,
                x = "logFoldChange",
                y = "pvalue",
                lab = rownames(de_res),
                FCcutoff = 0.5,
                pCutoff = 0.000001,
                labSize = 6,
                pointSize = 3,
                selectLab = c("type_578", "type_642", "type_970", "type_577", "type_513"),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE
                )

sig_de <- de_res[de_res$pvalue < 0.05 & !is.na(de_res$pvalue),]
sig_type <- all_types[rownames(all_types) %in% rownames(sig_de),]
rep_de <- sig_de[sig_de$foldChange<1,]
rep_type <- all_types[rownames(all_types) %in% rownames(rep_de),]
rep_rank <- rep_de[order(rep_de$foldChange),]
rep_rank <- rep_rank[rep_rank$pvalue<0.000001,]
rep_rank_type <- all_types[rownames(all_types) %in% rownames(rep_rank),]
e_type <- rep_type[rep_type$EBOX>0,]
#pheatmap(as.matrix(rep_rank_type), 
#         cluster_cols = FALSE, 
#         cluster_rows = FALSE,
#         color = c("#FFFFFF","#99FF99"), 
#         fontsize = 20,
#         show_rownames = TRUE, 
#         legend = FALSE,
#         labels_row = rownames(rep_rank_type))
library(ComplexHeatmap)
Heatmap(rep_rank_type,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        border = TRUE,
        col = c("#FFFFFF","#99FF99"),
        show_heatmap_legend = FALSE,
        column_names_gp = gpar(fontsize = 15),
        column_names_centered = TRUE,
        column_names_rot = 0,
        row_names_gp = gpar(col = ifelse(rownames(rep_rank_type) %in% c("type_578", "type_642", "type_970"), "red", "black"), fontsize = 15),
        row_order = rownames(rep_rank),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x, y, width, height, 
                    gp = gpar(col = "grey60", fill = NA, lwd = 0.5))
        }
        )

type_dmso <- all_types[rownames(all_types) %in% peak_count_pileup_dmso$type,]

is_combo_rep <- match_type(type_dmso,rep_rank_type)
supp_type_logit <- cbind(all_types,is_combo_rep)
sup_logistic <- glm(is_combo_rep ~ ., family = binomial(link = 'logit'), data = type_dmso)
summary(sup_logistic)

supp_co_type <- as.data.frame(cooccur_vec_from_binary(type_dmso))
supp_co_logistic <- glm(is_combo_rep ~ ., family = binomial(link = 'logit'), data = supp_co_type)
summary(supp_co_logistic)
accuracy <- sum((predict(supp_co_logistic, supp_co_type) > -1.682962) & is_combo_supp)/188
write.table(sort(supp_co_logistic$coefficients, decreasing = TRUE), "supp_co_logistic_coefficients.txt", quote = FALSE, sep = "\t")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("org.Hs.eg.db")
covplot(peak_dmso, chrs=c("chr17", "chr18"))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakHeatmap(peak_dmso, TxDb=txdb, upstream=3000, downstream=3000)
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
peakAnno <- annotatePeak(peak_dmso, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
annot_df <- as.data.frame(peakAnno)
ggscatter(annot_df, x = "type", y = "distanceToTSS")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

prom_peak <- annot_df[annot_df$annotation %in% c("Promoter (1-2kb)", "Promoter (<=1kb)"),]
quant_dmso <- quant_dmso[quant_dmso$name %in% prom_peak$name,]

peakAnno <- annotatePeak(peak_shp, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
annot_df <- as.data.frame(peakAnno)
prom_peak <- annot_df[annot_df$annotation %in% c("Promoter (1-2kb)", "Promoter (<=1kb)"),]
quant_shp <- quant_shp[quant_shp$name %in% prom_peak$name,]

peakAnno <- annotatePeak(peak_mg, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
annot_df <- as.data.frame(peakAnno)
prom_peak <- annot_df[annot_df$annotation %in% c("Promoter (1-2kb)", "Promoter (<=1kb)"),]
quant_mg <- quant_mg[quant_mg$name %in% prom_peak$name,]

peakAnno <- annotatePeak(peak_combo, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
annot_df <- as.data.frame(peakAnno)
prom_peak <- annot_df[annot_df$annotation %in% c("Promoter (1-2kb)", "Promoter (<=1kb)"),]
quant_combo <- quant_combo[quant_combo$name %in% prom_peak$name,]

rep_peaks <- peak_dmso[peak_dmso$type %in% rownames(rep_rank_type),]
rep_peak_anno <- annotatePeak(peak_combo, tssRegion=c(-2000, 2000),
                              TxDb=txdb, annoDb="org.Hs.eg.db")



plotPeakProf2(peak = rep_peaks, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, ignore_strand = F)
upsetplot(rep_peak_anno, vennpie=TRUE)
plotAnnoPie(rep_peak_anno)+
  theme(legend.text = element_text(size=25))
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(rep_peak_anno)$geneId)
head(pathway1, 2)
dotplot(pathway1)

library("plyranges")
mb231_atac <- read_bigwig("../GSM4557168_WT_atac.bw")
hits <- findOverlaps(rep_peaks, mb231_atac, type = "within")
rep_peaks_in_open <- rep_peaks[queryHits(hits)]
rep_peak_anno <- annotatePeak(peak_combo, tssRegion=c(-2000, 2000),
                              TxDb=txdb, annoDb="org.Hs.eg.db")



x <- peak_count_pileup_dmso
x$dot_color <- ifelse(peak_count_pileup_dmso$type %in% rownames(rep_rank), "red", "black")
ggscatter(x, x = "peak_count", y = "pileup",
          xlim = c(0,7000), 
          ylim = c(0,3e6),
          #label = "type",
          color = "dot_color",
          repel = TRUE)+
  xlab("Peak Counts") +
  ylab("Read Counts")
