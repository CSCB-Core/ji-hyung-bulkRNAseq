# I'd like to have a general transcriptomic analysis including 1)expression level of individual genes in each samples, 2) differentially expressed genes between groups, and 3) pathway analysis such as gene ontology test and gene set enrichment test.
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
library(Rsubread)

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
      "processed/01-corall_processing/run3/star_rsem/rsem.merged.gene_counts.tsv"
    )
  )

# tx_counts <-
#   fread(
#     here(
#       "processed/01-corall_processing/run3/star_rsem/rsem.merged.transcript_counts.tsv"
#     )
#   )


# dropping sample 10
gene_counts$`HTS019_10_JH133-M_LB_S10_L001` <- NULL

# gene_counts <- as.data.frame(gene_counts)

m_data <- gene_counts$gene_id

annots <- select(
  org.Mm.eg.db,
  keys = m_data,
  columns = c("SYMBOL", "ENSEMBLTRANS"),
  keytype = "ENSEMBL"
)

gene_counts_symbols <-
  merge(gene_counts, annots, by.x = "gene_id", by.y = "ENSEMBL")

gene_counts_symbols$means <- rowMeans(gene_counts[, 3:17])

gene_counts_symbols <-
  gene_counts_symbols[order(x = gene_counts_symbols$means, decreasing = TRUE), ]

# create a metadata table with the column names of the gene_counts table as row names
metadata <-
  data.frame(row.names = colnames(gene_counts_symbols[, 3:17]))

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
  "group4",
  "group4",
  "group4",
  "group4"
)

metadata$sample_names <- rownames(metadata)

# removing all of the SYMBOLS that did not have a corresponding Ensembl ID
gene_counts_symbols_noNA <- na.omit(gene_counts_symbols)

pdf(here("processed/02-DESeq2/EDA_Plots.pdf"))

# generate ggplot plot of top 10 genes by mean expression, ignore NAs
ggplot(gene_counts_symbols[1:10, ], aes(x = SYMBOL, y = means)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 10 genes by mean expression") +
  xlab("Gene") +
  ylab("Mean expression")


ggplot(gene_counts_symbols_noNA[1:10, ], aes(x = SYMBOL, y = means)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 10 genes by mean expression") +
  xlab("Gene") +
  ylab("Mean expression")

## EDA ####

ggplot(gene_counts_symbols) +
  geom_histogram(aes(x = means), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  labs(title = "Distribution of mean raw expression counts")

# show the same plot but with a log scale on the x-axis
ggplot(gene_counts_symbols) +
  geom_histogram(aes(x = means), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Log-transformed number of genes") +
  labs(title = "Distribution of mean raw expression counts") +
  scale_y_log10()

# gene_counts[,3:18] <- round(gene_counts[,3:18])

## Mean vs Variance ####
mean_counts <-
  apply(gene_counts[, 3:17], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns
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

de_gene_counts <- gene_counts[, 3:17] %>% as.matrix()

rownames(de_gene_counts) <- gene_counts$gene_id

dds <- DESeqDataSetFromMatrix(
  countData = round(de_gene_counts),
  # negatives?
  colData = metadata,
  design =  ~ group,
  tidy = FALSE
)

dds <- estimateSizeFactors(dds)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

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

dds <- DESeq(dds, quiet = FALSE)

# 1) Group1 vs Group2, 2) Group1 vs Group3, 3) Group 3 vs Group4, and 4) Group 2 vs Group 4.
# One control, one knockout, one control with 2nd condition and knockout with 2nd condition.
group1_vs_group2 <-
  results(dds,
          contrast = c("group", "group2", "group1"),
          alpha = 0.05) # need to modify the base levels in the third arguments here
group1_vs_group3 <-
  results(dds,
          contrast = c("group", "group3", "group1"),
          alpha = 0.05)
group3_vs_group4 <-
  results(dds,
          contrast = c("group", "group3", "group4"),
          alpha = 0.05)
group2_vs_group4 <-
  results(dds,
          contrast = c("group", "group2", "group4"),
          alpha = 0.05)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## Save the unshrunken results to compare
group1_vs_group2_unshrunken <- group1_vs_group2

resultsNames(dds)
# [1] "Intercept"              "group_group2_vs_group1" "group_group3_vs_group1" "group_group4_vs_group1"
# Apply fold change shrinkage
group1_vs_group2 <-
  lfcShrink(dds, coef = "group_group2_vs_group1", type = "apeglm")

## Plots ####
# normalized_counts <- counts(dds, normalized=T) %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>%
#   as_tibble() %>%
#   left_join(annots, by=c("gene" = "ENSEMBL"))
#
# norm_OEsig <- normalized_counts[,c(1:4,5:8)] %>%
#   filter(gene %in% rownames(group1_vs_group2))
# ### Set a color palette
# heat_colors <- brewer.pal(6, "YlOrRd")
#
# ### Run pheatmap using the metadata data frame for the annotation
# pheatmap(norm_OEsig[,2:ncol(norm_OEsig)],
#          color = heat_colors,
#          cluster_rows = T,
#          show_rownames = F,
#          annotation = meta2,
#          border_color = NA,
#          fontsize = 10,
#          fontsize_row = 10,
#          height = 20)

group1_vs_group3 <- group1_vs_group3 %>% as_tibble() %>%
  mutate(group1_vs_group3 = padj < 0.05 &
           abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(group1_vs_group3) +
  geom_point(aes(
    x = log2FoldChange,
    y = -log10(padj),
    colour = 0.05
  )) +
  ggtitle("Group 3 vs Group 1") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  )
