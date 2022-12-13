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

# > gene_counts[gene_counts$gene_name == "Tinf2",]
# gene_id gene_name HTS019_01_JH135.M_LU_S1_L001 HTS019_02_JH135.M_LB_S2_L001 HTS019_03_JH135.M_RU_S3_L001 HTS019_04_JH140.M_LB_S4_L001 HTS019_05_JH134.M_LU_S5_L001
# 1: ENSMUSG00000007589     Tinf2                           26                           17                            5                           11                           11
# HTS019_06_JH134.M_LB_S6_L001 HTS019_07_JH134.M_LR_S7_L001 HTS019_08_JH139.M_LU_S8_L001 HTS019_09_JH133.M_LU_S9_L001 HTS019_10_JH133.M_LB_S10_L001 HTS019_11_JH140.M_LU_S11_L001
# 1:                            3                           20                            7                       13.004                             6                         9.212
# HTS019_12_JH140.M_RU_S12_L001 HTS019_13_JH133.M_RU_S13_L001 HTS019_14_JH133.M_RB_S14_L001 HTS019_15_JH139.M_LB_S15_L001 HTS019_16_JH139.M_RU_S16_L001
# 1:                            15                             7                             3                             7                             3

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

## REMOVE SAMPLE 15
gene_counts <- gene_counts[,!"HTS019_15_JH139.M_LB_S15_L001"]

metadata <- metadata[-15,]

column_sums <-
  apply(gene_counts[, 3:ncol(gene_counts)], 2, function(x)
    sum(x > 0))

gene_counts %>% select_if(is.numeric) %>% colSums

gene_counts$means <- rowMeans(gene_counts[, 3:17])

gene_counts <-
  gene_counts[order(x = gene_counts$means, decreasing = TRUE),]

de_gene_counts <- gene_counts[, 3:17] %>% as.matrix()

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

# SCTC-JHL-49 - CONDITION 1 ONLY ####

metadata <- subset(metadata, select = -condition)

metadata <- metadata[-(9:15),]

de_gene_counts <- de_gene_counts[,-(9:15)]


# Prep DESeq2 ####
dds <- DESeqDataSetFromMatrix(
  countData = round(de_gene_counts),
  # negatives?
  colData = metadata,
  design = ~ 0 + genome,
  tidy = FALSE
)

dds <- estimateSizeFactors(dds)

## Vis ####

### Transform counts for data visualization
rld <- rlog(dds, blind = TRUE)

### Plot PCA

# save the PCA plot as a pdf
pdf(here("processed/02-DESeq2/PCA_GroupxSample_condition1_only.pdf"))

plotPCA(rld, intgroup = "group")
plotPCA(rld, intgroup = "sample_names")

graphics.off()

### Extract the rlog matrix from the object
rld_mat <- assay(rld)

rld_out <- as.data.frame(rld_mat)

colnames(rld_out) <- colnames(dds)
rownames(rld_out) <- rownames(dds)

#Tinf2 present
write.csv(rld_out,
          here(
            "processed/02-DESeq2/rlog2_counts_allSamples_JH_RNAseq_Condition1_only.csv"
          ))

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

pdf(
  here(
    "processed/02-DESeq2/CSCB_JHL_44_RNAseq_Analysis_sans15/SampleCorrelationHeatmap_condition1_only.pdf"
  ),
  height = 8,
  width = 11
)

### Plot heatmap using the correlation matrix and the metadata object
pheatmap(
  rld_cor,
  annotation = metadata,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Sample Correlation Heatmap"
)

pheatmap(
  rld_cor[1:8, 1:8],
  annotation = metadata,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Sample Correlation Heatmap, Condition 1"
)
# 
# pheatmap(
#   rld_cor[9:15, 9:15],
#   annotation = metadata,
#   cluster_rows = FALSE,
#   cluster_cols = FALSE,
#   main = "Sample Correlation Heatmap, Condition 2"
# )

graphics.off()

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# > counts(dds[rownames(dds) == "Tinf2",])
# HTS019_01_JH135.M_LU_S1_L001 HTS019_02_JH135.M_LB_S2_L001 HTS019_03_JH135.M_RU_S3_L001 HTS019_04_JH140.M_LB_S4_L001 HTS019_05_JH134.M_LU_S5_L001 HTS019_06_JH134.M_LB_S6_L001
# Tinf2                           26                           17                            5                           11                           11                            3
# HTS019_07_JH134.M_LR_S7_L001 HTS019_08_JH139.M_LU_S8_L001 HTS019_09_JH133.M_LU_S9_L001 HTS019_10_JH133.M_LB_S10_L001 HTS019_11_JH140.M_LU_S11_L001 HTS019_12_JH140.M_RU_S12_L001
# Tinf2                           20                            7                           13                             6                             9                            15
# HTS019_13_JH133.M_RU_S13_L001 HTS019_14_JH133.M_RB_S14_L001 HTS019_16_JH139.M_RU_S16_L001
# Tinf2                             7                             3                             3

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
            "genomeKO",
            "genomeWT"
          ))

write.csv(
  res_WTvKO_cond1,
  file = here(
    "processed/02-DESeq2/CSCB_JHL_44_RNAseq_Analysis_sans15/genomeKOvsWT_condition1_only.csv"
  )
)

res_WTvKO_cond1 <- res_WTvKO_cond1 %>% na.omit()

res_WTvKO_cond1 <- res_WTvKO_cond1[order(res_WTvKO_cond1$padj),]

### Condition 1 plots ####
# MA plot and volcano plot for each contrast
pdf(
  here(
    "processed/02-DESeq2/CSCB_JHL_44_RNAseq_Analysis_sans15/KOvsWT_cond1_DESeq2_plots.pdf"
  ),
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

# plotCounts(dds,
#            gene = which.min(res_WTvKO_cond1$padj),
#            intgroup = "genome")

graphics.off()

## Contrast Condition 2 ####
# 
# res_WTvKO_cond2 <-
#   results(dds,
#           contrast = list(
#             "genomeKO.conditioncondition2",
#             "genomeWT.conditioncondition2"
#           ))
# 
# write.csv(
#   res_WTvKO_cond2,
#   file = here(
#     "processed/02-DESeq2/CSCB_JHL_44_RNAseq_Analysis_sans15/genomeKOvsWT.conditioncondition2_sans15.csv"
#   )
# )
# 
# res_WTvKO_cond2 <- res_WTvKO_cond2 %>% na.omit()
# 
# res_WTvKO_cond2 <- res_WTvKO_cond2[order(res_WTvKO_cond2$padj),]

### Condition 2 Plots ####
# pdf(
#   here("processed/02-DESeq2/KOvsWT_cond2_DESeq2_plots.pdf"),
#   height = 8,
#   width = 11
# )
# 
# plotMA(res_WTvKO_cond2)
# 
# p_data_cond2 <- res_WTvKO_cond2 %>% na.omit() %>% as_tibble()
# 
# p_data_cond2$sig <- "NS"
# 
# p_data_cond2[p_data_cond2$padj < 0.1, ]$sig <- "A"
# 
# ## Volcano plot; color by significance (adjusted p-value)
# ggplot(p_data_cond2, aes(
#   x = log2FoldChange,
#   y = -log10(padj),
#   color = sig
# )) +
#   geom_point() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none") +
#   xlab("log2 Fold Change") +
#   ylab("-log10 Adjusted p-value") +
#   ggtitle("Volcano plot of KO vs WT in Condition 2 interaction") +
#   scale_color_manual(values = c("black", "blue"))
# 
# normalized_counts <- counts(dds, normalized = T) %>%
#   data.frame() %>%
#   rownames_to_column(var = "gene_name") %>%
#   as_tibble() %>%
#   filter(gene_name %in% rownames(res_WTvKO_cond2)[res_WTvKO_cond2$padj < 0.1])
# 
# heat_colors <- brewer.pal(6, "YlOrRd")
# 
# ### Run pheatmap using the metadata data frame for the annotation
# pheatmap(
#   normalized_counts[, which(metadata$condition == "condition2") + 1],
#   color = heat_colors,
#   cluster_rows = T,
#   show_rownames = F,
#   annotation = metadata,
#   border_color = NA,
#   fontsize = 10,
#   fontsize_row = 10,
#   height = 20,
#   scale = "row",
#   main = "Heatmap plot of Differentially-Expressed Genes with a < 0.1 Adjusted P-Value, KO vs WT in Condition 2 Interaction"
# ) # Row-scaled
# 
# # plotCounts(dds,
# #            gene = which.min(res_WTvKO_cond2$padj),
# #            intgroup = "genome")
# 
# graphics.off()

## Custom Constrasts ####

# 1) Group1 vs Group2, 2) Group1 vs Group3, 3) Group 3 vs Group4, and 4) Group 2 vs Group 4.
# One control, one knockout, one control with 2nd condition and knockout with 2nd condition.

### WT2vsWT1 ####
# res_WT2vsWT1 <-
#   results(dds,
#           contrast = list(
#             "genomeWT.conditioncondition2",
#             "genomeWT.conditioncondition1"
#           ))
# 
# write.csv(
#   res_WT2vsWT1,
#   file = here(
#     "processed/02-DESeq2/CSCB_JHL_44_RNAseq_Analysis_sans15/genomeWT.conditioncondition2vs1_sans15.csv"
#   )
# )
# 
# res_WT2vsWT1 <- res_WT2vsWT1 %>% na.omit()
# 
# res_WT2vsWT1 <- res_WT2vsWT1[order(res_WT2vsWT1$padj),]
# 
# pdf(
#   here("processed/02-DESeq2/WT2vsWT1_DESeq2_plots.pdf"),
#   height = 8,
#   width = 11
# )
# 
# plotMA(res_WT2vsWT1)
# 
# p_data_WT2vWT1 <- res_WT2vsWT1 %>% na.omit() %>% as_tibble()
# 
# p_data_WT2vWT1$sig <- "NS"
# 
# p_data_WT2vWT1[p_data_WT2vWT1$padj < 0.1, ]$sig <- "A"
# 
# ## Volcano plot; color by significance (adjusted p-value)
# ggplot(p_data_WT2vWT1, aes(
#   x = log2FoldChange,
#   y = -log10(padj),
#   color = sig
# )) +
#   geom_point() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none") +
#   xlab("log2 Fold Change") +
#   ylab("-log10 Adjusted p-value") +
#   ggtitle("Volcano plot of WT in Condition 2 vs WT in Condition 1") +
#   scale_color_manual(values = c("black", "blue"))
# 
# normalized_counts <- counts(dds, normalized = T) %>%
#   data.frame() %>%
#   rownames_to_column(var = "gene_name") %>%
#   as_tibble() %>%
#   filter(gene_name %in% rownames(res_WT2vsWT1)[res_WT2vsWT1$padj < 0.1])
# 
# heat_colors <- brewer.pal(6, "YlOrRd")
# 
# ### Run pheatmap using the metadata data frame for the annotation
# pheatmap(
#   normalized_counts[, which(metadata$genome == "WT") + 1],
#   # WT 1 vs WT 2
#   color = heat_colors,
#   cluster_rows = T,
#   show_rownames = F,
#   annotation = metadata,
#   border_color = NA,
#   fontsize = 10,
#   fontsize_row = 10,
#   height = 20,
#   scale = "row",
#   main = "Heatmap plot of Differentially-Expressed Genes with a < 0.1 Adjusted P-Value, WT 1 vs WT 2 Interaction"
# ) # Row-scaled
# 
# # plotCounts(dds,
# #            gene = which.min(res_WT2vsWT1$padj),
# #            intgroup = "genome")
# 
# graphics.off()

## KO2vsKO1 ####

# res_KO2vsKO1 <-
#   results(dds,
#           contrast = list(
#             "genomeKO.conditioncondition2",
#             "genomeKO.conditioncondition1"
#           ))
# 
# write.csv(
#   res_KO2vsKO1,
#   file = here(
#     "processed/02-DESeq2/CSCB_JHL_44_RNAseq_Analysis_sans15/genomeKO.conditioncondition2vs1_sans15.csv"
#   )
# )
# 
# res_KO2vsKO1 <- res_KO2vsKO1 %>% na.omit()
# 
# res_KO2vsKO1 <- res_KO2vsKO1[order(res_KO2vsKO1$padj),]
# 
# pdf(
#   here("processed/02-DESeq2/KO2vsKO1_DESeq2_plots.pdf"),
#   height = 8,
#   width = 11
# )
# 
# plotMA(res_KO2vsKO1)
# 
# p_data_KO2vsKO1 <- res_KO2vsKO1 %>% na.omit() %>% as_tibble()
# 
# p_data_KO2vsKO1$sig <- "NS"
# 
# p_data_KO2vsKO1[p_data_KO2vsKO1$padj < 0.1, ]$sig <- "A"
# 
# ## Volcano plot; color by significance (adjusted p-value)
# ggplot(p_data_KO2vsKO1, aes(
#   x = log2FoldChange,
#   y = -log10(padj),
#   color = sig
# )) +
#   geom_point() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none") +
#   xlab("log2 Fold Change") +
#   ylab("-log10 Adjusted p-value") +
#   ggtitle("Volcano plot of KO vs WT in Condition 2 interaction") +
#   scale_color_manual(values = c("black", "blue"))
# 
# normalized_counts <- counts(dds, normalized = T) %>%
#   data.frame() %>%
#   rownames_to_column(var = "gene_name") %>%
#   as_tibble() %>%
#   filter(gene_name %in% rownames(res_KO2vsKO1)[res_KO2vsKO1$padj < 0.1])
# 
# heat_colors <- brewer.pal(6, "YlOrRd")
# 
# ### Run pheatmap using the metadata data frame for the annotation
# pheatmap(
#   normalized_counts[, which(metadata$genome == "KO") + 1],
#   # WT 1 vs WT 2
#   color = heat_colors,
#   cluster_rows = T,
#   show_rownames = F,
#   annotation = metadata,
#   border_color = NA,
#   fontsize = 10,
#   fontsize_row = 10,
#   height = 20,
#   scale = "row",
#   main = "Heatmap plot of Differentially-Expressed Genes with a < 0.1 Adjusted P-Value, KO 1 vs KO 2 Interaction"
# ) # Row-scaled
# 
# # plotCounts(dds,
# #            gene = which.min(res_KO2vsKO1$padj),
# #            intgroup = "genome")
# 
# graphics.off()

# save.image(here("preGO_GSEA.RData"))
# load("preGO_GSEA.RData")

#----------------------------------------------
# PATHWAY ANALYSIS AND GENE SET ENRICHMENT ####
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html#visualize-enriched-go-terms-as-a-directed-acyclic-graph


gene_WTvKO_cond1 <-
  rownames(res_WTvKO_cond1)[(res_WTvKO_cond1$padj < 0.1) &
                              (res_WTvKO_cond1$log2FoldChange > 0)]

# gene_WTvKO_cond2 <-
#   rownames(res_WTvKO_cond2)[(res_WTvKO_cond2$padj < 0.1) &
#                               (res_WTvKO_cond2$log2FoldChange > 0)]
# 
# gene_KO2vsKO1 <-
#   rownames(res_KO2vsKO1)[(res_KO2vsKO1$padj < 0.1) &
#                            (res_KO2vsKO1$log2FoldChange > 0)]
# 
# gene_WT2vsWT1 <-
#   rownames(res_WT2vsWT1)[(res_WT2vsWT1$padj < 0.1) &
#                            (res_WT2vsWT1$log2FoldChange > 0)]

gene_WTvKO_cond1_sig_genes_tx <-
  gene_counts[gene_counts$gene_name %in% gene_WTvKO_cond1, ]$gene_id
# gene_WTvKO_cond2_sig_genes_tx <-
#   gene_counts[gene_counts$gene_name %in% gene_WTvKO_cond2, ]$gene_id
# gene_KO2vsKO1_sig_genes_tx <-
#   gene_counts[gene_counts$gene_name %in% gene_KO2vsKO1, ]$gene_id
# gene_WT2vsWT1_sig_genes_tx <-
#   gene_counts[gene_counts$gene_name %in% gene_WT2vsWT1, ]$gene_id

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

# enrichGO_gene_WTvKO_cond2 <- enrichGO(
#   gene = gene_WTvKO_cond2_sig_genes_tx,
#   universe = gene_counts$gene_id,
#   OrgDb = org.Mm.eg.db,
#   keyType = 'ENSEMBL',
#   ont = "BP",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.01,
#   qvalueCutoff = 0.1,
#   readable = TRUE
# )
# 
# enrichGO_gene_KO2vsKO1 <- enrichGO(
#   gene = gene_KO2vsKO1_sig_genes_tx,
#   universe = gene_counts$gene_id,
#   OrgDb = org.Mm.eg.db,
#   keyType = 'ENSEMBL',
#   ont = "BP",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.01,
#   qvalueCutoff = 0.1,
#   readable = TRUE
# )
# 
# enrichGO_gene_WT2vsWT1 <- enrichGO(
#   gene = gene_WT2vsWT1_sig_genes_tx,
#   universe = gene_counts$gene_id,
#   OrgDb = org.Mm.eg.db,
#   keyType = 'ENSEMBL',
#   ont = "BP",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.01,
#   qvalueCutoff = 0.1,
#   readable = TRUE
# )

prepGeneList <- function(DESeq_res, gene_counts, rank = "log2") {
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
  prepGeneList(res_WTvKO_cond1, gene_counts)

# geneList_res_WTvKO_cond2_input <-
#   prepGeneList(res_WTvKO_cond2, gene_counts)
# 
# geneList_res_WT2vsWT1_input <-
#   prepGeneList(res_WT2vsWT1, gene_counts)
# 
# geneList_res_KO2vsKO1_input <-
#   prepGeneList(res_KO2vsKO1, gene_counts)

# GSEA supports ENSEMBL IDs
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

# gsea_gene_WTvKO_cond2_sig_genes_tx <- gseGO(
#   geneList     = geneList_res_WTvKO_cond2_input,
#   OrgDb        = org.Mm.eg.db,
#   keyType = "ENSEMBL",
#   pAdjustMethod = "BH",
#   ont = "BP",
#   minGSSize    = 100,
#   maxGSSize    = 500,
#   pvalueCutoff = 0.10,
#   verbose      = TRUE
# )
# 
# gsea_gene_WT2vsWT1_sig_genes_tx <- gseGO(
#   geneList     = geneList_res_WT2vsWT1_input,
#   OrgDb        = org.Mm.eg.db,
#   keyType = "ENSEMBL",
#   pAdjustMethod = "BH",
#   ont = "BP",
#   minGSSize    = 100,
#   maxGSSize    = 500,
#   pvalueCutoff = 0.10,
#   verbose      = TRUE
# )
# 
# gsea_gene_KO2vsKO1_sig_genes_tx <- gseGO(
#   geneList     = geneList_res_KO2vsKO1_input,
#   OrgDb        = org.Mm.eg.db,
#   keyType = "ENSEMBL",
#   pAdjustMethod = "BH",
#   ont = "BP",
#   minGSSize    = 100,
#   maxGSSize    = 500,
#   pvalueCutoff = 0.10,
#   verbose      = TRUE
# )

# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

pdf(here("processed/02-DESeq2/GO_ORA-GSEA_Condition1_only.pdf"),
    height = 8,
    width = 11)


barplot(enrichGO_gene_WTvKO_cond1) + ggtitle(
  "Bar Plot of Top ORA Gene Ontologies in KO vs WT Condition 1 Differential Expression Analysis"
)
# barplot(enrichGO_gene_WTvKO_cond2) + ggtitle(
#   "Bar Plot of Top ORA Gene Ontologies in KO vs WT Condition 2 Differential Expression Analysis"
# )
# barplot(enrichGO_gene_WT2vsWT1) + ggtitle(
#   "Bar Plot of Top ORA Gene Ontologies in WT Condition 2 vs WT Condition 1 Differential Expression Analysis"
# )
# barplot(enrichGO_gene_KO2vsKO1) + ggtitle(
#   "Bar Plot of Top ORA Gene Ontologies in KO Condition 2 vs KO Condition 1 Differential Expression Analysis"
# )

dotplot(gsea_gene_WTvKO_cond1_sig_genes_tx) + ggtitle("Dot Plot of Top GO GSEA Gene Sets in KO vs WT Condition 1 Differential Expression Analysis")
# dotplot(gsea_gene_WTvKO_cond2_sig_genes_tx) + ggtitle("Dot Plot of Top GO GSEA Gene Sets in KO vs WT Condition 2 Differential Expression Analysis")
# dotplot(gsea_gene_WT2vsWT1_sig_genes_tx) + ggtitle(
#   "Dot Plot of Top GO GSEA Gene Sets in WT Condition 2 vs WT Condition 1 Differential Expression Analysis"
# )
# dotplot(gsea_gene_KO2vsKO1_sig_genes_tx) + ggtitle(
#   "Dot Plot of Top GO GSEA Gene Sets in KO Condition 2 vs WT Condition 1 Differential Expression Analysis"
# )

# emapplot(gsea_gene_WTvKO_cond1_sig_genes_tx, showCategory = 10)
#
# cnetplot(gsea_gene_WTvKO_cond1_sig_genes_tx, categorySize="pvalue", foldChange=geneList_res_WTvKO_cond1_input, showCategory = 5)


graphics.off()

## KEGG ####

KEGG_IDs_WTvKO1 <-
  bitr(
    names(geneList_res_WTvKO_cond1_input)[abs(geneList_res_WTvKO_cond1_input) > 1],
    # keep only logFC > 1
    fromType = 'ENSEMBL',
    toType = 'ENTREZID',
    OrgDb = org.Mm.eg.db
    # organism = "mmu"
  )

# KEGG_IDs_WTvKO2 <-
#   bitr(
#     names(geneList_res_WTvKO_cond2_input)[abs(geneList_res_WTvKO_cond2_input) > 1],
#     # keep only logFC > 1
#     fromType = 'ENSEMBL',
#     toType = 'ENTREZID',
#     OrgDb = org.Mm.eg.db
#     # organism = "mmu"
#   )
# 
# KEGG_IDs_WT2vWT1 <-
#   bitr(
#     names(geneList_res_WT2vsWT1_input)[abs(geneList_res_WT2vsWT1_input) > 1],
#     # keep only logFC > 1
#     fromType = 'ENSEMBL',
#     toType = 'ENTREZID',
#     OrgDb = org.Mm.eg.db
#     # organism = "mmu"
#   )
# 
# KEGG_IDs_KO2vsKO1 <-
#   bitr(
#     names(geneList_res_KO2vsKO1_input)[abs(geneList_res_KO2vsKO1_input) > 1],
#     # keep only logFC > 1
#     fromType = 'ENSEMBL',
#     toType = 'ENTREZID',
#     OrgDb = org.Mm.eg.db
#     # organism = "mmu"
#   )

# geneList_res_WTvKO_cond2 <- geneList_res_WTvKO_cond2[geneList_res_WTvKO_cond2$gene_id %in% gene_counts$gene_name,]

### WTvsKO1 ####
KEGG_IDs_WTvKO1 <-
  KEGG_IDs_WTvKO1[KEGG_IDs_WTvKO1$ENSEMBL %in% names(geneList_res_WTvKO_cond1_input),]

geneList_res_WTvKO_cond1_input_df <-
  geneList_res_WTvKO_cond1_input %>% as.data.frame() %>% rownames_to_column(var = "ENSEMBL")

geneList_res_WTvKO_cond1_input_df <-
  merge(KEGG_IDs_WTvKO1,
        geneList_res_WTvKO_cond1_input_df,
        by = "ENSEMBL")

colnames(geneList_res_WTvKO_cond1_input_df)[3] <-
  c("log2FoldChange")

geneList_res_WTvKO_cond1_input_df_input <-
  geneList_res_WTvKO_cond1_input_df$log2FoldChange

names(geneList_res_WTvKO_cond1_input_df_input) <-
  geneList_res_WTvKO_cond1_input_df$ENTREZID

geneList_res_WTvKO_cond1_input_df_input <-
  geneList_res_WTvKO_cond1_input_df_input[order(geneList_res_WTvKO_cond1_input_df_input, decreasing = TRUE)]

gseKEGG_WTvKO_cond1 <-
  gseKEGG(
    geneList     = geneList_res_WTvKO_cond1_input_df_input,
    organism     = 'mmu',
    # minGSSize    = 120,
    pvalueCutoff = 0.10,
    verbose      = TRUE
  )

### WTvsKO2 ####

# KEGG_IDs_WTvKO2 <-
#   KEGG_IDs_WTvKO2[KEGG_IDs_WTvKO2$ENSEMBL %in% names(geneList_res_WTvKO_cond2_input),]
# 
# geneList_res_WTvKO_cond2_input_df <-
#   geneList_res_WTvKO_cond2_input %>% as.data.frame() %>% rownames_to_column(var = "ENSEMBL")
# 
# geneList_res_WTvKO_cond2_input_df <-
#   merge(KEGG_IDs_WTvKO2,
#         geneList_res_WTvKO_cond2_input_df,
#         by = "ENSEMBL")
# 
# colnames(geneList_res_WTvKO_cond2_input_df)[3] <-
#   c("log2FoldChange")
# 
# geneList_res_WTvKO_cond2_input_df_input <-
#   geneList_res_WTvKO_cond2_input_df$log2FoldChange
# 
# names(geneList_res_WTvKO_cond2_input_df_input) <-
#   geneList_res_WTvKO_cond2_input_df$ENTREZID
# 
# geneList_res_WTvKO_cond2_input_df_input <-
#   geneList_res_WTvKO_cond2_input_df_input[order(geneList_res_WTvKO_cond2_input_df_input, decreasing = TRUE)]
# 
# gseKEGG_WTvKO_cond2 <-
#   gseKEGG(
#     geneList     = geneList_res_WTvKO_cond2_input_df_input,
#     organism     = 'mmu',
#     # minGSSize    = 120,
#     pvalueCutoff = 0.10,
#     verbose      = TRUE
#   )

### WT1vsWT2 ####

# KEGG_IDs_WT2vWT1 <-
#   KEGG_IDs_WT2vWT1[KEGG_IDs_WT2vWT1$ENSEMBL %in% names(geneList_res_WT2vsWT1_input),]
# 
# geneList_res_WT2vsWT1_input_df <-
#   geneList_res_WT2vsWT1_input %>% as.data.frame() %>% rownames_to_column(var = "ENSEMBL")
# 
# geneList_res_WT2vsWT1_input_df <-
#   merge(KEGG_IDs_WT2vWT1,
#         geneList_res_WT2vsWT1_input_df,
#         by = "ENSEMBL")
# 
# colnames(geneList_res_WT2vsWT1_input_df)[3] <- c("log2FoldChange")
# 
# geneList_res_WT2vsWT1_input_df_input <-
#   geneList_res_WT2vsWT1_input_df$log2FoldChange
# 
# names(geneList_res_WT2vsWT1_input_df_input) <-
#   geneList_res_WT2vsWT1_input_df$ENTREZID
# 
# geneList_res_WT2vsWT1_input_df_input <-
#   geneList_res_WT2vsWT1_input_df_input[order(geneList_res_WT2vsWT1_input_df_input, decreasing = TRUE)]
# 
# gseKEGG_WT2vsWT1 <-
#   gseKEGG(
#     geneList     = geneList_res_WT2vsWT1_input_df_input,
#     organism     = 'mmu',
#     # minGSSize    = 120,
#     pvalueCutoff = 0.10,
#     verbose      = TRUE
#   )

### KO1vsKO2 ####

# KEGG_IDs_KO2vsKO1 <-
#   KEGG_IDs_KO2vsKO1[KEGG_IDs_KO2vsKO1$ENSEMBL %in% names(geneList_res_KO2vsKO1_input),]
# 
# geneList_res_KO2vsKO1_input_df <-
#   geneList_res_KO2vsKO1_input %>% as.data.frame() %>% rownames_to_column(var = "ENSEMBL")
# 
# geneList_res_KO2vsKO1_input_df <-
#   merge(KEGG_IDs_KO2vsKO1,
#         geneList_res_KO2vsKO1_input_df,
#         by = "ENSEMBL")
# 
# colnames(geneList_res_KO2vsKO1_input_df)[3] <- c("log2FoldChange")
# 
# geneList_res_KO2vsKO1_input_df_input <-
#   geneList_res_KO2vsKO1_input_df$log2FoldChange
# 
# names(geneList_res_KO2vsKO1_input_df_input) <-
#   geneList_res_KO2vsKO1_input_df$ENTREZID
# 
# geneList_res_KO2vsKO1_input_df_input <-
#   geneList_res_KO2vsKO1_input_df_input[order(geneList_res_KO2vsKO1_input_df_input, decreasing = TRUE)]
# 
# gseKEGG_KO2vsKO1 <-
#   gseKEGG(
#     geneList     = geneList_res_KO2vsKO1_input_df_input,
#     organism     = 'mmu',
#     # minGSSize    = 120,
#     pvalueCutoff = 0.10,
#     verbose      = TRUE
#   )

### KEGG Plots ####

pdf(here("processed/02-DESeq2/KEGG-GSEA_Condition1_only.pdf"),
    height = 8,
    width = 11)

dotplot(
  gseKEGG_WTvKO_cond1,
  showCategory = 10,
  title = "Enriched Pathways" ,
  split = ".sign"
) + facet_grid(. ~ .sign) + ggtitle(
  "Dot Plot of Top KEGG Enriched Pathways in KO vs WT Condition 1 Differential Expression Analysis"
)

# dotplot(
#   gseKEGG_WTvKO_cond2,
#   showCategory = 10,
#   title = "Enriched Pathways" ,
#   split = ".sign"
# ) + facet_grid(. ~ .sign) + ggtitle(
#   "Dot Plot of Top KEGG Enriched Pathways in KO vs WT Condition 2 Differential Expression Analysis"
# ) # NO Enrichment Terms Found
# 
# dotplot(
#   gseKEGG_WT2vsWT1,
#   showCategory = 10,
#   title = "Enriched Pathways" ,
#   split = ".sign"
# ) + facet_grid(. ~ .sign) + ggtitle("Dot Plot of Top KEGG Enriched Pathways in WT 2 vs WT 1 Differential Expression Analysis")
# 
# dotplot(
#   gseKEGG_KO2vsKO1,
#   showCategory = 10,
#   title = "Enriched Pathways" ,
#   split = ".sign"
# ) + facet_grid(. ~ .sign) + ggtitle("Dot Plot of Top KEGG Enriched Pathways in KO 2 vs KO 1 Differential Expression Analysis")

graphics.off()

## Write CSVs for PCA and all enrichment objects

# GO
# write.csv(
#   enrichGO_gene_WT2vsWT1@result,
#   file = here("processed/02-DESeq2/enrichGO_gene_WT2vsWT1_results.csv")
# )
# 
# write.csv(
#   enrichGO_gene_KO2vsKO1@result,
#   file = here("processed/02-DESeq2/enrichGO_gene_KO2vsKO1_results.csv")
# )

write.csv(
  enrichGO_gene_WTvKO_cond1@result,
  file = here(
    "processed/02-DESeq2/enrichGO_gene_WTvKO_cond1_results.csv"
  )
)

# write.csv(
#   enrichGO_gene_WTvKO_cond2@result,
#   file = here(
#     "processed/02-DESeq2/enrichGO_gene_WTvKO_cond2_results.csv"
#   )
# )

# GSEA

# write.csv(
#   gsea_gene_WT2vsWT1_sig_genes_tx@result,
#   file = here("processed/02-DESeq2/gsea_gene_WT2vsWT1_sig_genes_tx.csv")
# )
# 
# write.csv(
#   gsea_gene_KO2vsKO1_sig_genes_tx@result,
#   file = here("processed/02-DESeq2/gsea_gene_KO2vsKO1_sig_genes_tx.csv")
# )

write.csv(
  gsea_gene_WTvKO_cond1_sig_genes_tx@result,
  file = here(
    "processed/02-DESeq2/gsea_gene_WTvKO_cond1_sig_genes_tx.csv"
  )
)

# write.csv(
#   gsea_gene_WTvKO_cond2_sig_genes_tx@result,
#   file = here(
#     "processed/02-DESeq2/gsea_gene_WTvKO_cond2_sig_genes_tx.csv"
#   )
# )

# KEGG_GSEA

# write.csv(gseKEGG_KO2vsKO1@result,
#           file = here("processed/02-DESeq2/gseKEGG_KO2vsKO1.csv"))
# 
# write.csv(gseKEGG_WT2vsWT1@result,
#           file = here("processed/02-DESeq2/gseKEGG_WT2vsWT1.csv"))

write.csv(
  gseKEGG_WTvKO_cond1@result,
  file = here("processed/02-DESeq2/gseKEGG_WTvKO_cond1.csv")
)
# 
# write.csv(
#   gseKEGG_WTvKO_cond2@result,
#   file = here("processed/02-DESeq2/gseKEGG_WTvKO_cond2.csv")
# )

# PCA

pca_table <- plotPCA(rld, intgroup = "group")$data

pca_table$group.1 <- NULL

write.csv(pca_table, file = "pca_table_condition1_only.csv")
