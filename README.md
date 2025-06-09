# PHGDH Deletion in Tumors: Single-Cell Transcriptomic Analysis of Coding and Non-Coding RNAs

## Background

**Phgdh (phosphoglycerate dehydrogenase)** is a key enzyme in the serine synthesis pathway and is implicated in tumor growth, metabolic reprogramming, and stress resistance. This study investigates the consequences of **Phgdh deletion** in the tumor context, focusing not only on changes in protein-coding gene expression but also on the emerging regulatory layer provided by **non-coding RNAs**.

**Non-coding RNAs (ncRNAs)**—including long non-coding RNAs (lncRNAs), microRNAs (miRNAs), and small nucleolar RNAs (snoRNAs)—do not encode proteins but play critical roles in gene regulation, chromatin structure, RNA stability, and cell fate decisions. While most single-cell transcriptomic studies focus on protein-coding transcripts, this analysis integrates ncRNA expression to capture a more comprehensive landscape of transcriptional regulation in tumors.

-   By combining both coding and non-coding transcripts, this study aims to:

-   Reveal tumor cell populations that are driven or distinguished by ncRNA activity.

-   Compare clustering patterns and marker gene signatures derived from coding vs. non-coding expression profiles.

-   Identify specific non-coding RNAs associated with Phgdh deletion and their potential role in tumor reprogramming.

## Methods

### Data Processing and Quality Control

-   Seurat was used to load and preprocess scRNA-seq data from selected GSMO and PHGDH samples.

-   Mitochondrial content and gene/UMI counts were used to filter low-quality cells.

-   Normalization, identification of highly variable genes, and dimensionality reduction (PCA, UMAP, t-SNE) were performed.

-   Louvain clustering was applied for unsupervised cluster identification.

### Non-Coding Gene Integration

-   GTF annotation (GRCm39.109) was parsed to classify genes into **coding** and **non-coding** categories (lncRNA, miRNA, snoRNA).

-   Gene biotype labels were added to Seurat objects using Seurat and SingleCellExperiment interfaces.

-   Per-cell average expression scores were calculated for coding and non-coding genes, enabling comparative UMAP visualizations and classification of cells as "Coding-High" or "NonCoding-High".

### Comparative Transcriptomic Analysis

-   Seurat objects were subset into **coding-only** and **non-coding-only** transcriptomes.

-   Independent clustering and UMAP analyses were conducted on both subsets to explore distinct structure and diversity.

-   Feature plots and violin plots were used to visualize key gene expression patterns.

-   Differential gene expression analysis was performed separately for coding and non-coding subsets, and top marker genes per cluster were identified.

## Results

**Cell Clustering and Annotation**:

-   Clusters were annotated to reflect distinct tumor states (e.g., quiescent, mitotic, differentiating), as well as stromal and immune populations.

**Non-Coding Transcript Usage**:

-   UMAP overlays revealed tumor subpopulations with higher non-coding gene expression.

-   A subset of tumor cells classified as “NonCoding-High,” suggesting potential functional roles for lncRNAs or miRNAs in the Phgdh-deleted context.

**Independent Analysis of Gene Types**:

-   Clustering using non-coding genes alone revealed distinct structures from the coding transcriptome, demonstrating the potential of non-coding profiles in capturing unique biological variance.

-   Several non-coding RNAs emerged as top differentially expressed features within tumor subpopulations

**Coding vs. Non-Coding Influence**:

-   Parallel analysis highlighted key differences in gene regulation and cell-type distribution driven by gene biotype.

-   The comparison supports the hypothesis that **Phgdh loss impacts both metabolic and non-coding transcriptional networks**.

    ![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/PC_elbow_plot.png?raw=true){width="600"}

    ![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/scatter_rnd_QC.png?raw=true){width="1200"}

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/tsne_louvain_res_p3.png?raw=true){width="500"}

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/umap_louvain_res_p3.png?raw=true){width="500"}

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/seurat_obj_mouse_with_gene_type_03_violin_exp2.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/seurat_obj_mouse_with_gene_type_03_feature_exp.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/umap_coding_vs_noncoding_subset.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/umap_coding_vs_noncoding_overlay.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/feature_plot_coding_vs_noncoding_02.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/var_features.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/var_features_noncoding.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/seurat_obj_mouse_with_gene_type_03_dotplot_top5.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/dotplot_top5_coding.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh/dotplot_top5_noncoding.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh_02/umap_coding_vs_noncoding_subset.png?raw=true)

![](https://github.com/chingyaousf/PHGDH-Deletion-in-Tumors-Single-Cell-Transcriptomic-Analysis-of-Coding-and-Non-Coding-RNAs/blob/main/plots/phgdh_02/cluster_correspondence_heatmap.png?raw=true)
