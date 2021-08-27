if (! require(edgeR)) {
   source("https://bioconductor.org/biocLite.R")
   biocLite("edgeR")
   library(edgeR)
}

if (! require(DESeq2)) {
   source("https://bioconductor.org/biocLite.R")
   biocLite("DESeq2")
   library(DESeq2)
}

data = read.table("/scratch/pawsey0149/hhu/mainwork/bar/barley_pantranscriptome/EC_RNA/fastq/clean/sort/count/genes.counts.matrix", header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,5,6)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("ECI-2-120_C", 3), rep("ECI-2-120_T", 3))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","ECI-2-120_C","ECI-2-120_T")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "ECI-2-120_C"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "ECI-2-120_T"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="ECI-2-120_C", sampleB="ECI-2-120_T", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])
write.table(res, file='genes.counts.matrix.ECI-2-120_C_vs_ECI-2-120_T.DESeq2.DE_results', sep=' ', quote=FALSE)
write.table(rnaseqMatrix, file='genes.counts.matrix.ECI-2-120_C_vs_ECI-2-120_T.DESeq2.count_matrix', sep='      ', quote=FALSE)
source("/scratch/pawsey0149/hhu/tool/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf("genes.counts.matrix.ECI-2-120_C_vs_ECI-2-120_T.DESeq2.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
dev.off()
