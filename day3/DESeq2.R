setwd("~/day3")
# system("mkdir -p rsem")
# system("ln -s ~kjyi/day3/rsem/* ./rsem")
library(tidyverse)
library(DESeq2)
files <- dir("rsem", full.names = TRUE)
rsem.tximport <- tximport::tximport(files, type = "rsem")
sample_names <- basename(files) %>% stringr::str_replace(".rsem.genes.results", "")
condition <- stringr::str_extract(sample_names, "[A-Z]*$") %>% factor()
samples <- data.frame(sample_names, condition)
# change effective length 0 regions to 1 to avoid error
rsem.tximport$length[rsem.tximport$length == 0] <- 1
dds <- DESeqDataSetFromTximport(rsem.tximport, colData = samples, design = ~condition)

# Prefiltering low count genes < 10
keep <- rowSums(counts(dds)) >= 10; table(keep)
dds <- dds[keep, ]

# Differential expression analysis

dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "condition_TIL_vs_NTIL", type = "apeglm")
resLFC


res05 <- results(dds, alpha = 0.05)
summary(res05)



library(IHW)
resIHW <- results(dds, filterFun = ihw)
summary(resIHW)

# plotting
plotMA(res, ylim = c(-4,4))
plotMA(resLFC, ylim = c(-5,5))
plotCounts(dds, 
           gene = which.min(res$padj),
           intgroup = "condition")
res@rownames[which.min(res$padj)]
# TAGLN2, a smooth muscle marker

# Count data transformations
# rlog and vst
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
head(assay(vsd), 3)

library("pheatmap")
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("condition","type")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)
select <- dds %>%
  counts(normalized = TRUE) %>% 
  rowMeans() %>%
  order(decreasing = TRUE) %>% 
  .[1:20]

select <- order(res$padj[which(res$padj < 0.05)], decreasing = TRUE)

res %>% as.data.frame() %>%
  rownames_to_column("geneID") %>%
  select(geneID, log2FoldChange , padj) %>%
  na.omit() %>%
  filter(padj < 0.001) %>%
  arrange(desc(log2FoldChange)) %>%
  .$geneID -> selected

topVarGenes <- assay(vsd) %>% rowVars() %>% order(decreasing = TRUE) %>% head(100)

qval0.001_order_LFC <- as.data.frame(res) %>%
  rownames_to_column("geneID") %>%
  select(geneID, log2FoldChange , padj) %>%
  filter(padj < 0.001) %>%
  arrange(desc(log2FoldChange)) %>%
  .$geneID


df <- as.data.frame(colData(dds)[,c("condition")])

pheatmap(assay(vsd)[topVarGenes,], cluster_rows = TRUE, show_rownames = FALSE,
         cluster_cols = TRUE, annotation_col = df)
pheatmap(assay(vsd)[qval0.001_order_LFC,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = TRUE, annotation_col = df)
x <- assay(vsd)

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method = "min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
    df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

x_quantile_normalized <- quantile_normalisation(x)
pheatmap(x_quantile_normalized[qval0.001_order_LFC,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = TRUE, annotation_col = df)
pheatmap(scale(x)[qval0.001_order_LFC,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = TRUE, annotation_col = df)
# same code - pheatmap(... , scale=row)
pheatmap(t(scale(t(x)))[qval0.001_order_LFC,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = TRUE, annotation_col = df)
# same code - pheatmap(... , scale=column)

# limit log value to -3 ~ 3
# change distance function to Pearson correlation than Euclidean
# change color to blue-black-yelow
t(scale(t(x)))[qval0.001_order_LFC,] %>%
{.[. > 3] <- 3;.[. < -3] <- -3; .} %>%
pheatmap(cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = TRUE, annotation_col = df, 
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
t(scale(t(x)))[qval0.001_order_LFC,] %>%
{.[. > 3] <- 3;.[. < -3] <- -3; .} %>%
  pheatmap(cluster_rows = FALSE, show_rownames = FALSE,
           cluster_cols = TRUE, annotation_col = df, 
           clustering_distance_cols = "correlation",
           color = colorRampPalette(c("#29469E", "grey15", "yellow"))(50))
