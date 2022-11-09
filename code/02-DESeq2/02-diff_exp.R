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

setwd("/mnt/isilon/cscb/Projects/lexogen/lee_ji-hyung/SCTC-JHL-44/lee_bulkRNA")

here::i_am("code/02-DESeq2/02-diff_exp.R")

gene_counts <-
  fread(
    here(
      "processed/01-corall_processing/run3/star_rsem/rsem.merged.gene_counts.tsv"
    )
  )

tx_counts <-
  fread(
    here(
      "processed/01-corall_processing/run3/star_rsem/rsem.merged.transcript_counts.tsv"
    )
  )

gene_counts <- as.data.frame(gene_counts)

m_data <- gene_counts$gene_id

annots <- select(
  org.Mm.eg.db,
  keys = m_data,
  columns = c("SYMBOL", "ENSEMBLTRANS"),
  keytype = "ENSEMBL"
)

gene_counts <-
  merge(gene_counts, annots, by.x = "gene_id", by.y = "ENSEMBL")

gene_counts$means <- rowMeans(gene_counts[, 3:18])

gene_counts <-
  gene_counts[order(x = gene_counts$means, decreasing = TRUE), ]

# create a metadata table with the column names of the gene_counts table as row names
metadata <- data.frame(row.names = colnames(gene_counts[, 3:18]))

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

# for (i in 1:nrow(gene_counts)) {
#   if (is.na(gene_counts$SYMBOL[i])) {
#     transcript_ids <- strsplit(gene_counts$`transcript_id(s)`[i], ",")[[1]]
#     for (j in 1:length(transcript_ids)) {
#       annots <- select(org.Mm.eg.db,
#                        keys = transcript_ids,
#                        columns = "SYMBOL",
#                        keytype = "ENSEMBL")
#       if (!is.na(annots$SYMBOL)) {
#         gene_counts$SYMBOL[i] <- annots$SYMBOL
#         break
#       }
#     } # end for j
#   } # end if
# } # end for i


# removing all of the SYMBOLS that did not have a corresponding Ensembl ID
gene_counts_noNA <- na.omit(gene_counts)

#find the 10 most highly expressed genes per sample (gene_counts[, 3:18])
top10 <-
  apply(gene_counts[, 3:18], 2, function(x)
    head(order(x, decreasing = TRUE), 10))

# generate ggplot plot of top 10 genes by mean expression, ignore NAs
ggplot(gene_counts[1:10, ], aes(x = SYMBOL, y = means)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 10 genes by mean expression") +
  xlab("Gene") +
  ylab("Mean expression")


ggplot(gene_counts_noNA[1:10, ], aes(x = SYMBOL, y = means)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 10 genes by mean expression") +
  xlab("Gene") +
  ylab("Mean expression")


gene_counts$gene_id <-
  make.names(gene_counts$gene_id, unique = TRUE)

de_gene_counts <- gene_counts[, 3:18]

rownames(de_gene_counts) <- gene_counts$gene_id

# set all columns of de_gene_counts to integers
de_gene_counts <- as.matrix(de_gene_counts)

# https://lashlock.github.io/compbio/R_presentation.html
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/01b_DGE_setup_and_overview.html
dds <- DESeqDataSetFromMatrix(
  countData = round(de_gene_counts), # negatives?
  colData = metadata,
  design =  ~ group,
  tidy = FALSE
)
