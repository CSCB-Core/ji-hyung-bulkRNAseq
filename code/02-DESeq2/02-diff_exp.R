# I'd like to have a general transcriptomic analysis including
# 1)expression level of individual genes in each samples,
# 2) differentially expressed genes between groups, and
# 3) pathway analysis such as gene ontology test and gene set enrichment test.
#
# I have 16 samples from 4 different groups and each group has 4 samples.
#
#  I'd like to compare the gene expression between 1) Group1 vs Group2, 2) Group1 vs Group3, 3) Group 3 vs Group4, and 4) Group 2 vs Group 4.
#
# It would be perfect if I can get general graphic images from analysis such as volcano plot, PCA analysis, heatmaps and bar graph of affected pathways.
#
# Please let me know if you need more information. I am also happy to get your advice and comments on analysis if you have some.
### Bioconductor and CRAN libraries used
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
# I had to install the dependency, lasso2, from the archive
# https://cran.r-project.org/src/contrib/Archive/lasso2/
# https://stackoverflow.com/questions/24194409/how-do-i-install-a-package-that-has-been-archived-from-cran
library(DEGreport) # lasso2
library(tximport)
library(ggplot2)
library(ggrepel)
library(here)
library(org.Mm.eg.db)
library(data.table)
library(clusterProfiler)

setwd("/mnt/isilon/cscb/Projects/lexogen/lee_ji-hyung/SCTC-JHL-44/lee_bulkRNA")

here::i_am("code/02-DESeq2/02-diff_exp.R")

counts <- list()

## Pre-EDA ####
gene_counts <-
  fread(
    here(
      "processed/01-corall_processing/run4/star_salmon/salmon.merged.gene_counts.tsv"
    )
  )

gene_counts$means <- rowMeans(gene_counts[, 3:18])

gene_counts <-
  gene_counts[order(x = gene_counts$means, decreasing = TRUE),]

# create a metadata table with the column names of the gene_counts table as row names
metadata <-
  data.frame(row.names = colnames(gene_counts[, 3:18]))

# add the metadata
metadata$group <- c(
  "group1",
  "group1",
  "group1",
  "group1",
  "group2",
  "group2",
  "group2",
  "group2",
  "group3",
  "group3",
  "group3",
  "group3",
  "group4",
  "group4",
  "group4",
  "group4"
)

metadata$sample_names <- rownames(metadata)

# removing all of the SYMBOLS that did not have a corresponding Ensembl ID

pdf(here("processed/02-DESeq2/EDA_Plots.pdf"))

# generate ggplot plot of top 10 genes by mean expression, ignore NAs
ggplot(gene_counts[1:10,], aes(x = gene_name, y = means)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 10 genes by mean expression") +
  xlab("Gene") +
  ylab("Mean expression")

## EDA ####

ggplot(gene_counts) +
  geom_histogram(aes(x = means), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  labs(title = "Distribution of mean raw expression counts")

# show the same plot but with a log scale on the x-axis
ggplot(gene_counts) +
  geom_histogram(aes(x = means), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Log-transformed number of genes") +
  labs(title = "Distribution of mean raw expression counts") +
  scale_y_log10()

## Mean vs Variance ####
mean_counts <-
  apply(gene_counts[, 3:18], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns
variance_counts <- apply(gene_counts[, 3:17], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x = mean_counts, y = variance_counts)) +
  scale_y_log10(limits = c(1, 1e9)) +
  scale_x_log10(limits = c(1, 1e9)) +
  geom_abline(intercept = 0,
              slope = 1,
              color = "red") +
  labs(title = "Mean vs Variance") +
  xlab("Mean") +
  ylab("Variance")

graphics.off()

# https://lashlock.github.io/compbio/R_presentation.html
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/01b_DGE_setup_and_overview.html
# https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html
# sum(duplicated(gene_counts$gene_id))

# gene_counts$gene_id <-
#   make.names(gene_counts$gene_id, unique = TRUE)

column_sums <-
  apply(gene_counts[, 3:ncol(gene_counts)], 2, function(x)
    sum(x > 0))

gene_counts %>% select_if(is.numeric) %>% colSums

de_gene_counts <- gene_counts[, 3:18] %>% as.matrix()

rownames(de_gene_counts) <- gene_counts$gene_name

# 1) Group1 vs Group2, 2) Group1 vs Group3, 3) Group 3 vs Group4, and 4) Group 2 vs Group 4.
# One control, one knockout, one control with 2nd condition and knockout with 2nd condition.

metadata$genome <- "WT"

metadata$genome[metadata$group %in% c('group2', 'group4')] <- 'KO'

metadata$condition[metadata$group %in% c('group1', 'group2')] <-
  'condition1'
metadata$condition[metadata$group %in% c('group3', 'group4')] <-
  'condition2'

metadata$condition <- as.factor(metadata$condition)
metadata$genome <- as.factor(metadata$genome)

metadata$genome <- relevel(metadata$genome, ref = "WT")

metadata$group <- as.factor(metadata$group)

dds <- DESeqDataSetFromMatrix(
  countData = round(de_gene_counts),
  # negatives?
  colData = metadata,
  design = ~ 0 + genome:condition,
  tidy = FALSE
)

dds <- estimateSizeFactors(dds)

## Vis ####

### Transform counts for data visualization
rld <- rlog(dds, blind = TRUE)

### Plot PCA

# save the PCA plot as a pdf
pdf(here("processed/02-DESeq2/PCA_GroupxSample.pdf"))

plotPCA(rld, intgroup = "group")
plotPCA(rld, intgroup = "sample_names")

graphics.off()

### Extract the rlog matrix from the object
rld_mat <- assay(rld)

rld_out <- as.data.frame(rld_mat)

colnames(rld_out) <- colnames(dds)
rownames(rld_out) <- rownames(dds)

write.csv(rld_out,
          here(
            "processed/02-DESeq2/rlog2_counts_allSamples_JH_RNAseq.csv"
          ))

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

pdf(here("processed/02-DESeq2/SampleCorrelationHeatmap.pdf"))

### Plot heatmap using the correlation matrix and the metadata object
pheatmap(
  rld_cor,
  annotation = metadata,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Sample Correlation Heatmap"
)

graphics.off()

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds, quiet = FALSE)

plotDispEsts(dds)

# 1) Group1 vs Group2, 2) Group1 vs Group3, 3) Group 3 vs Group4, and 4) Group 2 vs Group 4.
# One control, one knockout, one control with 2nd condition and knockout with 2nd condition.

# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## Change this up a bit ####

# > resultsNames(dds)
# [1] "genomeWT.conditioncondition1" "genomeKO.conditioncondition1" "genomeWT.conditioncondition2" "genomeKO.conditioncondition2"

## Contrast condition 1 ####

res_WTvKO_cond1 <-
  results(dds,
          contrast = list(
            "genomeKO.conditioncondition1",
            "genomeWT.conditioncondition1"
          ))

res_WTvKO_cond1 <- res_WTvKO_cond1 %>% na.omit()

res_WTvKO_cond1 <- res_WTvKO_cond1[order(res_WTvKO_cond1$padj),]

### Condition 1 plots ####
# MA plot and volcano plot for each contrast
pdf(
  here("processed/02-DESeq2/WTvKO_cond1_DESeq2_plots.pdf"),
  height = 8,
  width = 11
)

plotMA(res_WTvKO_cond1)

p_data_cond1 <- res_WTvKO_cond1 %>% na.omit() %>% as_tibble()

p_data_cond1$sig <- "NS"

p_data_cond1[p_data_cond1$padj < 0.1, ]$sig <- "A"

## Volcano plot; color by significance (adjusted p-value)
ggplot(p_data_cond1, aes(
  x = log2FoldChange,
  y = -log10(padj),
  color = sig
)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted p-value") +
  ggtitle("Volcano plot of KO vs WT in Condition 1 interaction") +
  scale_color_manual(values = c("black", "blue"))

normalized_counts <- counts(dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  as_tibble() %>%
  filter(gene_name %in% rownames(res_WTvKO_cond1)[res_WTvKO_cond1$padj < 0.1])

heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(
  normalized_counts[, which(metadata$condition == "condition1") + 1],
  color = heat_colors,
  cluster_rows = T,
  show_rownames = F,
  annotation = metadata,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  height = 20,
  scale = "row",
  main = "Heatmap plot of Differentially-Expressed Genes with a < 0.1 Adjusted P-Value, KO vs WT in Condition 1 Interaction"
) # Row-scaled

plotCounts(dds,
           gene = which.min(res_WTvKO_cond1$padj),
           intgroup = "genome")

graphics.off()

write.csv(
  res_WTvKO_cond1,
  file = here("processed/02-DESeq2/genomeKOvsWT.conditioncondition1.csv")
)

## Contrast Condition 2 ####

res_WTvKO_cond2 <-
  results(dds,
          contrast = list(
            "genomeKO.conditioncondition2",
            "genomeWT.conditioncondition2"
          ))

res_WTvKO_cond2 <- res_WTvKO_cond2 %>% na.omit()

res_WTvKO_cond2 <- res_WTvKO_cond2[order(res_WTvKO_cond2$padj),]

### Condition 2 Plots ####
pdf(
  here("processed/02-DESeq2/WTvKO_cond2_DESeq2_plots.pdf"),
  height = 8,
  width = 11
)

plotMA(res_WTvKO_cond2)

p_data_cond2 <- res_WTvKO_cond2 %>% na.omit() %>% as_tibble()

p_data_cond2$sig <- "NS"

p_data_cond2[p_data_cond2$padj < 0.1, ]$sig <- "A"

## Volcano plot; color by significance (adjusted p-value)
ggplot(p_data_cond2, aes(
  x = log2FoldChange,
  y = -log10(padj),
  color = sig
)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted p-value") +
  ggtitle("Volcano plot of KO vs WT in Condition 2 interaction") +
  scale_color_manual(values = c("black", "blue"))

normalized_counts <- counts(dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  as_tibble() %>%
  filter(gene_name %in% rownames(res_WTvKO_cond2)[res_WTvKO_cond2$padj < 0.1])

heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(
  normalized_counts[, which(metadata$condition == "condition2") + 1],
  color = heat_colors,
  cluster_rows = T,
  show_rownames = F,
  annotation = metadata,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  height = 20,
  scale = "row",
  main = "Heatmap plot of Differentially-Expressed Genes with a < 0.1 Adjusted P-Value, KO vs WT in Condition 2 Interaction"
) # Row-scaled

plotCounts(dds,
           gene = which.min(res_WTvKO_cond2$padj),
           intgroup = "genome")

graphics.off()

write.csv(
  res_WTvKO_cond2,
  file = here("processed/02-DESeq2/genomeKOvsWT.conditioncondition2.csv")
)

## Custom Constrasts ####

# 1) Group1 vs Group2, 2) Group1 vs Group3, 3) Group 3 vs Group4, and 4) Group 2 vs Group 4.
# One control, one knockout, one control with 2nd condition and knockout with 2nd condition.

### WT2vsWT1 ####
res_WT2vsWT1 <-
  results(dds,
          contrast = list(
            "genomeWT.conditioncondition2",
            "genomeWT.conditioncondition1"
          ))

res_WT2vsWT1 <- res_WT2vsWT1 %>% na.omit()

res_WT2vsWT1 <- res_WT2vsWT1[order(res_WT2vsWT1$padj),]

pdf(
  here("processed/02-DESeq2/WT2vsWT1_DESeq2_plots.pdf"),
  height = 8,
  width = 11
)

plotMA(res_WT2vsWT1)

p_data_WT2vWT1 <- res_WT2vsWT1 %>% na.omit() %>% as_tibble()

p_data_WT2vWT1$sig <- "NS"

p_data_WT2vWT1[p_data_WT2vWT1$padj < 0.1, ]$sig <- "A"

## Volcano plot; color by significance (adjusted p-value)
ggplot(p_data_WT2vWT1, aes(
  x = log2FoldChange,
  y = -log10(padj),
  color = sig
)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted p-value") +
  ggtitle("Volcano plot of WT in Condition 2 vs WT in Condition 1") +
  scale_color_manual(values = c("black", "blue"))

normalized_counts <- counts(dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  as_tibble() %>%
  filter(gene_name %in% rownames(res_WT2vsWT1)[res_WT2vsWT1$padj < 0.1])

heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(
  normalized_counts[, which(metadata$genome == "WT") + 1],
  # WT 1 vs WT 2
  color = heat_colors,
  cluster_rows = T,
  show_rownames = F,
  annotation = metadata,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  height = 20,
  scale = "row",
  main = "Heatmap plot of Differentially-Expressed Genes with a < 0.1 Adjusted P-Value, WT 1 vs WT 2 Interaction"
) # Row-scaled

plotCounts(dds,
           gene = which.min(res_WT2vsWT1$padj),
           intgroup = "genome")

graphics.off()

write.csv(res_WT2vsWT1,
          file = here("processed/02-DESeq2/genomeWT.conditioncondition2vs1.csv"))

## KO2vsKO1 ####

res_KO2vsKO1 <-
  results(dds,
          contrast = list(
            "genomeKO.conditioncondition2",
            "genomeKO.conditioncondition1"
          ))

res_KO2vsKO1 <- res_KO2vsKO1 %>% na.omit()

res_KO2vsKO1 <- res_KO2vsKO1[order(res_KO2vsKO1$padj),]

pdf(
  here("processed/02-DESeq2/KO2vsKO1_DESeq2_plots.pdf"),
  height = 8,
  width = 11
)

plotMA(res_KO2vsKO1)

p_data_KO2vsKO1 <- res_KO2vsKO1 %>% na.omit() %>% as_tibble()

p_data_KO2vsKO1$sig <- "NS"

p_data_KO2vsKO1[p_data_KO2vsKO1$padj < 0.1, ]$sig <- "A"

## Volcano plot; color by significance (adjusted p-value)
ggplot(p_data_KO2vsKO1, aes(
  x = log2FoldChange,
  y = -log10(padj),
  color = sig
)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted p-value") +
  ggtitle("Volcano plot of KO vs WT in Condition 2 interaction") +
  scale_color_manual(values = c("black", "blue"))

normalized_counts <- counts(dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  as_tibble() %>%
  filter(gene_name %in% rownames(res_KO2vsKO1)[res_KO2vsKO1$padj < 0.1])

heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(
  normalized_counts[, which(metadata$genome == "KO") + 1],
  # WT 1 vs WT 2
  color = heat_colors,
  cluster_rows = T,
  show_rownames = F,
  annotation = metadata,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  height = 20,
  scale = "row",
  main = "Heatmap plot of Differentially-Expressed Genes with a < 0.1 Adjusted P-Value, KO 1 vs KO 2 Interaction"
) # Row-scaled

plotCounts(dds,
           gene = which.min(res_KO2vsKO1$padj),
           intgroup = "genome")

graphics.off()

write.csv(res_KO2vsKO1,
          file = here("processed/02-DESeq2/genomeKO.conditioncondition2vs1.csv"))

# save.image(here("preGO_GSEA.RData"))

#----------------------------------------------
# PATHWAY ANALYSIS AND GENE SET ENRICHMENT ####
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html#visualize-enriched-go-terms-as-a-directed-acyclic-graph


gene_WTvKO_cond1 <-
  rownames(res_WTvKO_cond1)[(res_WTvKO_cond1$padj < 0.1) &
                              (res_WTvKO_cond1$log2FoldChange > 0)]

gene_WTvKO_cond2 <-
  rownames(res_WTvKO_cond2)[(res_WTvKO_cond2$padj < 0.1) &
                              (res_WTvKO_cond2$log2FoldChange > 0)]

gene_KO2vsKO1 <-
  rownames(res_KO2vsKO1)[(res_KO2vsKO1$padj < 0.1) &
                           (res_KO2vsKO1$log2FoldChange > 0)]

gene_WT2vsWT1 <-
  rownames(res_WT2vsWT1)[(res_WT2vsWT1$padj < 0.1) &
                           (res_WT2vsWT1$log2FoldChange > 0)]

# gene.df <- bitr(
#   gene_WTvKO_cond2,
#   fromType = "SYMBOL",
#   toType = c("ENSEMBL"),
#   OrgDb = org.Mm.eg.db
# )

gene_WTvKO_cond1_sig_genes_tx <-
  gene_counts[gene_counts$gene_name %in% gene_WTvKO_cond1, ]$gene_id
gene_WTvKO_cond2_sig_genes_tx <-
  gene_counts[gene_counts$gene_name %in% gene_WTvKO_cond2, ]$gene_id
gene_KO2vsKO1_sig_genes_tx <-
  gene_counts[gene_counts$gene_name %in% gene_KO2vsKO1, ]$gene_id
gene_WT2vsWT1_sig_genes_tx <-
  gene_counts[gene_counts$gene_name %in% gene_WT2vsWT1, ]$gene_id

enrichGO_gene_WTvKO_cond2 <- enrichGO(
  gene = gene_WTvKO_cond2_sig_genes_tx,
  universe = gene_counts$gene_id,
  OrgDb = org.Mm.eg.db,
  keyType = 'ENSEMBL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.1,
  readable = TRUE
)
enrichGO_gene_WTvKO_cond1 <- enrichGO(
  gene = gene_WTvKO_cond1_sig_genes_tx,
  universe = gene_counts$gene_id,
  OrgDb = org.Mm.eg.db,
  keyType = 'ENSEMBL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.1,
  readable = TRUE
)
enrichGO_gene_KO2vsKO1 <- enrichGO(
  gene = gene_KO2vsKO1_sig_genes_tx,
  universe = gene_counts$gene_id,
  OrgDb = org.Mm.eg.db,
  keyType = 'ENSEMBL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.1,
  readable = TRUE
)
enrichGO_gene_WT2vsWT1 <- enrichGO(
  gene = gene_WT2vsWT1_sig_genes_tx,
  universe = gene_counts$gene_id,
  OrgDb = org.Mm.eg.db,
  keyType = 'ENSEMBL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.1,
  readable = TRUE
)

prepGeneList <- function(DESeq_res, gene_counts, rank = "padj") {
  geneList <-
    DESeq_res %>% as.data.frame() %>% rownames_to_column("gene_id")
  
  # geneList_res_WTvKO_cond2 <- geneList_res_WTvKO_cond2[geneList_res_WTvKO_cond2$gene_id %in% gene_counts$gene_name,]
  
  gene_count_sub <-
    gene_counts[gene_counts$gene_name %in% geneList$gene_id,][, c("gene_name", "gene_id")] %>% as.data.frame()
  
  geneList <-
    merge(gene_count_sub,
          geneList,
          by.x = "gene_name",
          by.y = "gene_id")
  
  if (rank == "padj") {
    geneList <-
      geneList[, c("gene_id", "padj")]
    geneList_input <-
      geneList$padj
  } else if (rank == "log2") {
    geneList <-
      geneList[, c("gene_id", "log2FoldChange")]
    geneList_input <-
      geneList$log2FoldChange
  }
  
  names(geneList_input) <-
    geneList$gene_id
  
  geneList_input <-
    geneList_input[order(geneList_input, decreasing = TRUE)]
  
  return(geneList_input)
}

geneList_res_WTvKO_cond1_input <-
  prepGeneList(res_WTvKO_cond1, gene_counts, rank = "log2")
geneList_res_WTvKO_cond2_input <-
  prepGeneList(res_WTvKO_cond2, gene_counts)
geneList_res_WT2vsWT1_input <-
  prepGeneList(res_WT2vsWT1, gene_counts)
geneList_res_KO2vsKO1_input <-
  prepGeneList(res_KO2vsKO1, gene_counts)


gsea_gene_WTvKO_cond1_sig_genes_tx <- gseGO(
  geneList     = geneList_res_WTvKO_cond1_input,
  OrgDb        = org.Mm.eg.db,
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  ont = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.10,
  verbose      = TRUE
)

gsea_gene_WTvKO_cond2_sig_genes_tx <- gseGO(
  geneList     = geneList_res_WTvKO_cond2_input,
  OrgDb        = org.Mm.eg.db,
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  ont = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.10,
  verbose      = TRUE
)

gsea_gene_WT2vsWT1_sig_genes_tx <- gseGO(
  geneList     = geneList_res_WT2vsWT1_input,
  OrgDb        = org.Mm.eg.db,
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  ont = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.10,
  verbose      = TRUE
)

gsea_gene_KO2vsKO1_sig_genes_tx <- gseGO(
  geneList     = geneList_res_KO2vsKO1_input,
  OrgDb        = org.Mm.eg.db,
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  ont = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.10,
  verbose      = TRUE
)


pdf(
  here(
    "processed/02-DESeq2/GO_Enrichment_Overrrepresentation.pdf"
  ),
  height = 16,
  width = 22
)

goplot(enrichGO_gene_WTvKO_cond1)
goplot(enrichGO_gene_WTvKO_cond2)
goplot(enrichGO_gene_KO2vsKO1)
goplot(enrichGO_gene_WT2vsWT1)

graphics.off()

pdf(
  here(
    "processed/02-DESeq2/GO_GSEA_Pathways.pdf"
  ),
  height = 16,
  width = 22
)

goplot(gsea_gene_WTvKO_cond1_sig_genes_tx)
goplot(gsea_gene_WTvKO_cond2_sig_genes_tx)
goplot(gsea_gene_WT2vsWT1_sig_genes_tx) # NA
goplot(gsea_gene_KO2vsKO1_sig_genes_tx) # NA

graphics.off()

## KEGG ####
