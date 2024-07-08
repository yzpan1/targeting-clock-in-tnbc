library("GenomeInfoDb")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("EnhancedVolcano")
#################  feature counting using "GenomicAlignments"  package  #################
directory <- "~/Documents/projects/20220923_MB231_MG132/"
csvfile <- file.path(".", "sample_names_MB231.txt")
sampleTable <- read.csv(csvfile, row.names = 1)

# Generate summarizedExperiment
count_path <-file.path(".","MB231_simple_counts.txt")
simple_counts <- read.csv(count_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# Create DESeq class object
dds <- DESeqDataSetFromMatrix(countData = simple_counts, colData = sampleTable, design = ~ condition)
# Only save non-empty rows
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds$condition <- relevel(dds$condition, ref = "DMSO")

# Run Differential Expression
dds <- DESeq(dds)
shp <- results(dds, name = "condition_SHP_vs_DMSO")
mg <- results(dds, name = "condition_MG132_vs_DMSO")
combo <-results(dds, name = "condition_Combo_vs_DMSO")

shpShrunk <- lfcShrink(dds, coef = "condition_SHP_vs_DMSO")
mgShrunk <- lfcShrink(dds, coef = "condition_MG132_vs_DMSO")
comboShrunk <- lfcShrink(dds, coef = "condition_Combo_vs_DMSO")
shpSig <- subset(shpShrunk, padj < 0.05)
mgSig <- subset(mgShrunk, padj < 0.05)
comboSig <- subset(comboShrunk, padj < 0.05)

# Plot counts of single gene
geneCounts <- plotCounts(dds, gene = "DBP", intgroup = c("condition"), returnData = TRUE)
ggplot(geneCounts, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(width=.1,height=0), size=3)
# Plot MA-Plot, using shrunk log2FoldChange is recommended
plotMA(shpShrunk, ylim = c(-8,8))
plotMA(mgShrunk, ylim = c(-8,8))
plotMA(comboShrunk, ylim = c(-8,8))
#rlog transformation to make the data homoskedastic for PCA, may also use variance stabilizing transformation 
rld <- rlog(dds)
head(assay(rld),3)

# Plot Volcano with EnhancedVolcano" package
EnhancedVolcano(comboShrunk,
                lab = rownames(shpShrunk),
                x = 'log2FoldChange',
                y = 'pvalue')


# Plot sample distance heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Run )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
# Plot MDS
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=condition)) + geom_point(size=3)
# Plot PCA
plotPCA(rld, intgroup = "condition")

# Gene Clustering
library("genefilter")
topVarGenes <- head(order(rowVars(assay(shprld)),decreasing=TRUE),50)
coreClockGebes <- c('ARNTL','ARNTL2','CLOCK','NPAS2','CRY1','CRY2','PER1','PER2','PER3','NR1D1','NR1D2',
                    'RORA','RORB','RORC','DBP')
mat <- assay(rld)[ comboshp_names, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[c("Run","condition")])
pheatmap(mat, annotation_col=df)

comboshp_names <- scan("Combo+SHP_names.txt", what = "")
comboshp_names <- as.data.frame(comboshp)
write.csv(comboshp, file = "combo+shp_names.csv", row.names = F, quote = F)

comboshpname = c("MB231.Combo.1.sorted.bam","MB231.Combo.2.sorted.bam", "MB231.Combo.3.sorted.bam",
                 "MB231.SHP.1.sorted.bam", "MB231.SHP.2.sorted.bam", "MB231.SHP.3.sorted.bam")
shpname = c("MB231.DMSO.1.sorted.bam","MB231.DMSO.2.sorted.bam", "MB231.DMSO.3.sorted.bam",
            "MB231.SHP.1.sorted.bam", "MB231.SHP.2.sorted.bam", "MB231.SHP.3.sorted.bam")
simple_counts_comboshp <- subset(simple_counts, select = colnames(simple_counts) %in% comboshpname )
simple_count_shp <- subset(simple_counts, select = colnames(simple_counts) %in% shpname )
sample_table_comboshp <- sampleTable[-(4:9),]
sample_table_shp <- sampleTable[-(7:9),]
sample_table_shp <- sample_table_shp[-(1:3),]
rownames(sample_table_shp) <- shpname
dds_shp <- DESeqDataSetFromMatrix(countData = simple_count_shp, 
                                       colData = sample_table_shp, design = ~ condition)
dds_shp <- dds_shp[ rowSums(counts(dds_shp)) > 10, ]
shprld <- rlog(dds_shp)

mat <- assay(shprld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(shprld)[c("Run","condition")])
pheatmap(mat, annotation_col=df)


#################  Export CSV Files  #################
shpSigByP <- shpSig[order(shpSig$padj),]
shpSigByPDF <- as.data.frame(shpSigByP)
shpSigByLFC <- shpSig[order(shpSig$log2FoldChange),]
shpSigByLFC <- as.data.frame(shpSigByLFC)
# Column name of gene is not included. Add "Symbol" manually in the resulting csv file is needed to running gsea.r 
shpSigName <- as.data.frame(row.names(shpSigByLFC))
write.csv(shpSigName, file = "shpDESigName.csv", row.names = F, quote = F)
write.csv(shpSigByP, file = "shpDESigByP.csv")
write.csv(shpSigByLFC, file = "shpDESigbyLFC.csv")

mgSigByP <- mgSig[order(mgSig$padj),]
mgSigByPDF <- as.data.frame(mgSigByP)
mgSigByLFC <- mgSig[order(mgSig$log2FoldChange),]
mgSigByLFC <- as.data.frame(mgSigByLFC)
write.csv(mgSigByP, file = "mgDESigByP.csv")
write.csv(mgSigByLFC, file = "mgDESigByLFC.csv")
mgSigName <- as.data.frame(row.names(mgSig))
write.csv(mgSigName, file = "mgDESigName.csv", row.names = F, quote = F)

comboSigByP <- comboSig[order(comboSig$padj),]
comboSigByPDF <- as.data.frame(comboSigByP)
comboSigByLFC <- comboSig[order(comboSig$log2FoldChange),]
comboSigByLFC <- as.data.frame(comboSigByLFC)
write.csv(comboSigByP, file = "comboDESigByP.csv")
write.csv(comboSigByLFC, file = "comboDESigByLFC.csv")
comboSigName <- as.data.frame(row.names(comboSig))
write.csv(comboSigName, file = "comboDESigName.csv", row.names = F, quote = F)


comboOnly <- file.path(".", "ComboOnly_names.txt")
ComboOnly_names <- scan(comboOnly, what = "")
write.csv(comboOnly_names, file = "comboOnly_names.csv", row.names = F, quote = F)

save.image(file = "MB231_MG132.RData")
