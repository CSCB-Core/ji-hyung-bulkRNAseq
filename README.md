# lee_bulkRNA CSCB Bulk Cell RNA-Seq and Differential Expression Analysis

- [lee_bulkRNA CSCB Bulk Cell RNA-Seq and Differential Expression Analysis](#lee_bulkrna-cscb-bulk-cell-rna-seq-and-differential-expression-analysis)
    - [MultiQC](#multiqc)
    - [DESeq2](#deseq2)
    - [GO/KEGG enrichment and GSEA](#gokegg-enrichment-and-gsea)

The following is a brief report that covers the attached plots and datasets produced for Ji-Hyung's bulk RNA-seq project, internally coded as `SCTC-JHL-44`. 

The input to this pipeline consisted of 16 samples coming from 4 groups, overall representing 2 genotypes and 2 conditions:

 - "HTS019_01_JH135.M_LU_S1_L001"
 - "HTS019_02_JH135.M_LB_S2_L001"
 - "HTS019_03_JH135.M_RU_S3_L001"
 - "HTS019_04_JH140.M_LB_S4_L001"
 - "HTS019_05_JH134.M_LU_S5_L001"
 - "HTS019_06_JH134.M_LB_S6_L001"
 - "HTS019_07_JH134.M_LR_S7_L001"
 - "HTS019_08_JH139.M_LU_S8_L001"
 - "HTS019_09_JH133.M_LU_S9_L001"
 - "HTS019_10_JH133.M_LB_S10_L001" 
 - "HTS019_11_JH140.M_LU_S11_L001" 
 - "HTS019_12_JH140.M_RU_S12_L001"
 - "HTS019_13_JH133.M_RU_S13_L001"
 - "HTS019_14_JH133.M_RB_S14_L001"
 - "HTS019_15_JH139.M_LB_S15_L001"
 - "HTS019_16_JH139.M_RU_S16_L001"

HTS019_15_JH139.M_LB_S15_L001 was excluded from the differential expression analysis due to exceptionally poor read quality.

More details about the sample metadata can be found in `JHL_44_RNAseq_sample_information.xlsx`

### MultiQC

The MultiQC report was generated from the standardized and pre-packaged RNA-seq QC and analysis pipeline, [nf-core/rnaseq][def]. STAR and Salmon were used for alignment and quantification against Ensembl's GRCm38 reference index. Those raw gene-by-sample counts were then used for our custom downstream analyses. Extracting `run4_mouse_star_salmon.zip` and opening `multiqc_report.html` give more details about the quality of the sequenced samples and percent aligned reads. 

### DESeq2

The counts generated from nf-core were brought into R for the purposes of differential expression and GO enrichment and GSEA.

The [R-log2-transformed counts][def4] were saved for the collaborator's own edification here:
* rlog2_counts_allSamples_JH_RNAseq.csv

DESeq2 contrasted the differential expression as such:
* Knockout condition 2 / knockout condition 1
* Knockout condition 1 / wildtype condition 1
* Knockout condition 2 / wildtype condition 2
* Wildtype condition 2 / wildtype condition 1

CSVs were generated following the differential expression calculation by [DESeq2][def2]:
* genomeKO.conditioncondition2vs1.csv
* genomeKOvsWT.conditioncondition1.csv
* genomeKOvsWT.conditioncondition2.csv
* genomeWT.conditioncondition2vs1.csv

PCA and Sample correlation heatmap:
* PCA_GroupxSample.pdf
* SampleCorrelationHeatmap.pdf

[MA plots][def3], volcano plots, and heatmaps were generated according to each contrasts:
* KO2vsKO1_DESeq2_plots.pdf
* KOvsWT_cond1_DESeq2_plots.pdf
* KOvsWT_cond2_DESeq2_plots.pdf
* WT2vsWT1_DESeq2_plots.pdf

### GO/KEGG enrichment and GSEA
The R package `clusterProfiler` was used to find gene ontology enrichment, GSEA, and KEGG GSEA.
* GO_ORA-GSEA.pdf
* KEGG-GSEA.pdf


[def]: https://nf-co.re/rnaseq/3.9
[def2]: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
[def3]: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot
[def4]: https://rdrr.io/bioc/DESeq2/man/rlog.html
