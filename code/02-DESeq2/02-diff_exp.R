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

bams <- Sys.glob(
  here(
    "processed/01-corall_processing/run3/star_rsem/HTS019_*_S*_L001.markdup.sorted.bam"
  )
)

for (i in range(1:length(bams))) {
  counts[[i]] <-
    featureCounts(
      files = bams[i],
      isPairedEnd = TRUE,
      requireBothEndsMapped = TRUE,
      nthreads = 4
    )
}
## Pre-EDA ####
gene_counts <-
  fread(
    here(
      "processed/01-corall_processing/run4/star_salmon/salmon.merged.gene_counts.tsv"
    )
  )

# tx_counts <-
#   fread(
#     here(
#       "processed/01-corall_processing/run3/star_rsem/rsem.merged.transcript_counts.tsv"
#     )
#   )


# dropping sample 10
# gene_counts$`HTS019_10_JH133-M_LB_S10_L001` <- NULL

# gene_counts <- as.data.frame(gene_counts)

gene_counts$means <- rowMeans(gene_counts[, 3:18])

gene_counts <-
  gene_counts[order(x = gene_counts$means, decreasing = TRUE), ]

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
ggplot(gene_counts[1:10, ], aes(x = gene_name, y = means)) +
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
  design =  ~ genome * condition,
  tidy = FALSE
)

dds <- estimateSizeFactors(dds)

sum(apply(counts(dds), 1, function(x)
  sum(x) < 10))

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

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

pdf(here("processed/02-DESeq2/SampleCorrelationHeatmap.pdf"))

### Plot heatmap using the correlation matrix and the metadata object
pheatmap(rld_cor, annotation = metadata)

graphics.off()

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds, quiet = FALSE)

plotDispEsts(dds)

# 1) Group1 vs Group2, 2) Group1 vs Group3, 3) Group 3 vs Group4, and 4) Group 2 vs Group 4.
# One control, one knockout, one control with 2nd condition and knockout with 2nd condition.

# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## Save the unshrunken results to compare

resultsNames(dds)

res <-
  results(dds, name = c("condition_condition2_vs_condition1")) # need to modify the base levels in the third arguments here

# Apply fold change shrinkage


plotMA(res)

# > resultsNames(dds)
# [1] "Intercept"                          "genome_KO_vs_WT"                    "condition_condition2_vs_condition1" "genomeKO.conditioncondition2"

res_int <- results(dds, name = c("genomeKO.conditioncondition2"))

sum(res_int$padj < 0.1, na.rm = TRUE)

png(here(
  "processed/02-DESeq2/genomeKO.conditioncondition2_plotMA.png"
))

plotMA(res_int)

dev.off()

res_int_lfcShrink <-
  lfcShrink(dds, coef = "genomeKO.conditioncondition2", type = "apeglm")

png(here(
  "processed/02-DESeq2/genomeKO.conditioncondition2_LFCShrink_plotMA.png"
))

plotMA(res_int_lfcShrink)

dev.off()

res_int_lfcShrink_noNA <- na.omit(res_int_lfcShrink)

res_int_lfcShrink_noNA <-
  res_int_lfcShrink_noNA[order(res_int_lfcShrink_noNA$padj),]

sum(res_int_lfcShrink_noNA$padj < 0.1)

write.csv(res_int_lfcShrink_noNA, file = "processed/02-DESeq2/genomeKO.conditioncondition2_Interaction_naOmit.csv")

## Plots ####
p_data <- res_int %>% na.omit() %>% as_tibble()

p_data$sig <- "NS"

p_data[p_data$padj < 0.1, ]$sig <- "A"

# TODO make red into black; keep sig as blue
# TODO update the title and stuff
# TODO Save plot
## Volcano plot
ggplot(p_data) +
  geom_point(aes(
    x = log2FoldChange,
    y = -log10(padj),
    colour = sig
  )) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  )

# TODO do heatmap stuff later
normalized_counts <- counts(dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  as_tibble()

norm_OEsig <- normalized_counts %>%
  filter(gene %in% rownames(group1_vs_group2))

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(
  normalized_counts,
  color = heat_colors,
  cluster_rows = T,
  show_rownames = F,
  annotation = metadata,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  height = 20
)

group1_vs_group3 <- group1_vs_group3 %>% as_tibble() %>%
  mutate(group1_vs_group3 = padj < 0.05 &
           abs(log2FoldChange) >= 0.58)

## Change this up a bit ####
# model.matrix(~ 1 + genome:condition, metadata)

dds2 <- DESeqDataSetFromMatrix(
  countData = round(de_gene_counts),
  # negatives?
  colData = metadata,
  design = ~ 0 + genome:condition,
  tidy = FALSE
)

keep2 <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]

dds2 <- DESeq(dds2, quiet = FALSE)

resultsNames(dds2)

# > resultsNames(dds2)
# [1] "genomeWT.conditioncondition1" "genomeKO.conditioncondition1" "genomeWT.conditioncondition2" "genomeKO.conditioncondition2"

## Contrast condition 1 ####

res_WTvKO <-
  results(dds2,
          contrast = list(
            "genomeKO.conditioncondition1",
            "genomeWT.conditioncondition1"
          ))

res_WTvKO <- res_WTvKO %>% na.omit()

res_WTvKO <- res_WTvKO[order(res_WTvKO$padj), ]

# TODO Do heatmap, MA plot and volcano plot for each contrast

colData(dds2) <- metadata

plotCounts(dds2, gene = which.min(res_WTvKO$padj), intgroup = "genome")

write.csv(res_WTvKO,
          file = here("processed/02-DESeq2/genomeKOvsWT.conditioncondition1.csv"))

## Contrast Condition 2 ####

res_WTvKO_cond2 <-
  results(dds2,
          contrast = list(
            "genomeKO.conditioncondition2",
            "genomeWT.conditioncondition2"
          ))

res_WTvKO_cond2 <- res_WTvKO_cond2 %>% na.omit()

res_WTvKO_cond2 <- res_WTvKO_cond2[order(res_WTvKO_cond2$padj), ]

# TODO Do heatmap, MA plot and volcano plot for each contrast

plotCounts(dds2,
           gene = which.min(res_WTvKO_cond2$padj),
           intgroup = "genome")

write.csv(
  res_WTvKO_cond2,
  file = here("processed/02-DESeq2/genomeKOvsWT.conditioncondition2.csv")
)

# PATHWAY ANALYSIS AND GENE SET ENRICHMENT ####

gene_WTvKO_cond2 <-
  rownames(res_WTvKO_cond2)[(res_WTvKO_cond2$padj < 0.1) &
                              (res_WTvKO_cond2$log2FoldChange > 0)]

rownames(res_WTvKO_cond2)

# gene.df <- bitr(
#   gene_WTvKO_cond2,
#   fromType = "SYMBOL",
#   toType = c("ENSEMBL"),
#   OrgDb = org.Mm.eg.db
# )

ego <- enrichGO(
  gene = gene.df$ENSEMBL,
  universe = rownames(dds2),
  OrgDb = org.Mm.eg.db,
  keyType = 'ENSEMBL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.1,
  readable = TRUE
)
