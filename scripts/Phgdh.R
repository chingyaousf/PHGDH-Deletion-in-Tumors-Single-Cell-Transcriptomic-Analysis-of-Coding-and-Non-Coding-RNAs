# Load necessary libraries
library(Seurat)

# Define the file path
seurat_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_GSMO_PHGDH_genotypes.RDS"

# Load the Seurat object
seurat_object <- readRDS(seurat_path)

# Check object structure
seurat_object

# Check metadata
head(seurat_object@meta.data)

# Check available assays
Assays(seurat_object)

# Check gene names
head(rownames(seurat_object))


# Check the total number of unique samples
table(seurat_object$sample)

# Check the distribution of RNA counts and features
summary(seurat_object$nFeature_RNA)
summary(seurat_object$nCount_RNA)

# Check if any metadata columns indicate treatment groups
table(seurat_object$treatment)  # If this doesn't exist, you'll define it next


# subset the Seurat object to keep only GSMO and PHGDH samples (excluding PSD samples)
mouse <- subset(seurat_object, 
                  sample %in% c("GSMO_1", "GSMO_2", "GSMO_3", 
                                "PHGDH_1", "PHGDH_2", "PHGDH_4"))

# Define treatment groups based on sample names
mouse$treatment <- substr(mouse$sample, 1, 2)
mouse$treatment[mouse$treatment == "GS"] <- "GSMO"
mouse$treatment[mouse$treatment == "PH"] <- "PHGDH"
mouse$treatment[mouse$treatment == "PS"] <- "PHGDH_SHMT1"

# Set treatment as identity
Idents(mouse) <- "treatment"

# Confirm groups
table(mouse$treatment)



## QUALITY CONTROL, NORMALISATION AND SCALING

# Define figure saving path
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"

# Function to save figures as PNG only
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  
  # Save as PNG
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

## QC Plots before filtering and scaling

# Compute mitochondrial percentage
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^mt-")

# Cell quality control (by sample)
plot <- VlnPlot(mouse, pt.size = 0, alpha = 0.99,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, group.by = 'sample')
SaveFigure(plot, "vln_QC_sample", width = 16, height = 6)

# Cell quality control (by treatment)
plot <- VlnPlot(mouse, pt.size = 0,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, group.by = 'treatment')
SaveFigure(plot, "vln_QC_treatment", width = 16, height = 6)

# Scatter plots for RNA count relationships
plot1 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "percent.mt", 
                        shuffle = TRUE, pt.size = 0.1)
plot2 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        shuffle = TRUE, pt.size = 0.1)
SaveFigure((plot1 + plot2), "scatter_rnd_QC", width = 12, height = 6, res = 200)

# Perform the filtering
macfeat = median(mouse$nFeature_RNA) + sd(mouse$nFeature_RNA) * 3
mouse <- subset(mouse, subset = nFeature_RNA > 200 & nFeature_RNA < macfeat & percent.mt < 5)

# Confirm number of cells after filtering
dim(mouse)

# Normalizing the data
mouse <- NormalizeData(mouse, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mouse), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouse)
plot2 <- LabelPoints(plot = plot1, points = top10, xnudge = 0, ynudge = 0, repel = TRUE)
SaveFigure((plot1 + plot2), "var_features", width = 12, height = 6)

# Scaling the data
mouse <- ScaleData(mouse)

# Save the Seurat object after scaling
seurat_save_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_GSMO_PHGDH_genotypes_02.RDS"
saveRDS(mouse, file = seurat_save_path)

# Load the Seurat object
mouse <- readRDS("/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_GSMO_PHGDH_genotypes_02.RDS")


## Identify Non-Coding Genes
# Use an Annotation File to Extract Non-Coding Genes
#Since Andrey used Mus_musculus.GRCm39.109, download the annotation file and extract non-coding genes.


# Load necessary libraries
library(data.table)
library(openxlsx)  # For saving Excel files

# Define paths
output_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"

# Define the GTF file path
gtf_file <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/Mus_musculus.GRCm39.109.gtf.gz"


# Read the GTF file (without comment.char argument)
gene_anno <- fread(gtf_file, header = FALSE, sep = "\t")

# Remove comment lines manually (lines starting with "#")
gene_anno <- gene_anno[!grepl("^#", gene_anno$V1), ]

# Print first few rows to check
head(gene_anno)



# Load necessary libraries
library(data.table)
library(openxlsx)  # For saving Excel files

# Define paths
output_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"

# Define the GTF file path
gtf_file <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/Mus_musculus.GRCm39.109.gtf.gz"


# Read the GTF file (without comment.char argument)
gene_anno <- fread(gtf_file, header = FALSE, sep = "\t")

# Remove comment lines manually (lines starting with "#")
gene_anno <- gene_anno[!grepl("^#", gene_anno$V1), ]

# Print first few rows to check
head(gene_anno)



# Extract only gene entries (avoid exons, transcripts, etc.)
gene_info <- gene_anno[V3 == "gene", c(1,9)]
colnames(gene_info) <- c("chromosome", "info")

# Extract gene names safely
gene_info$gene_name <- ifelse(grepl('gene_name "', gene_info$info),
                              sub('.*gene_name "(.*?)";.*', "\\1", gene_info$info),
                              NA)  # Assign NA if gene_name is missing

# Extract gene biotypes safely
gene_info$gene_biotype <- ifelse(grepl('gene_biotype "', gene_info$info),
                                 sub('.*gene_biotype "(.*?)";.*', "\\1", gene_info$info),
                                 NA)  # Assign NA if gene_biotype is missing

# Remove rows where gene_name is NA
gene_info <- gene_info[!is.na(gene_info$gene_name), ]

# Remove duplicate entries
gene_info <- unique(gene_info)

# Display first few rows of extracted gene data
cat("Head of gene_info:\n")
print(head(gene_info))

# Show total number of genes extracted
cat("\nTotal number of genes extracted:", nrow(gene_info), "\n")

# Save `gene_info` as an Excel file (overwrite existing file)
write.xlsx(gene_info, file = paste0(output_path, "gene_info.xlsx"))
cat("Saved gene_info.xlsx to:", output_path, "\n")



## Filter for non-coding genes (lncRNA, miRNA, snoRNA)
non_coding_genes <- unique(gene_info$gene_name[gene_info$gene_biotype %in% c("lncRNA", "miRNA", "snoRNA")])

# Convert `non_coding_genes` to a data frame before saving
non_coding_genes_df <- data.frame(non_coding_genes)

# Display first few non-coding genes
cat("\nHead of non-coding genes:\n")
print(head(non_coding_genes_df))

# Show total number of non-coding genes extracted
cat("\nTotal number of non-coding genes:", nrow(non_coding_genes_df), "\n")

# Save `non_coding_genes` as an Excel file
write.xlsx(non_coding_genes_df, file = paste0(output_path, "non_coding_genes.xlsx"))
cat("Saved non_coding_genes.xlsx to:", output_path, "\n")



###  Classifies genes as either "coding" or "non_coding" and stores this in metadata
# Define path to non-coding gene list
non_coding_genes_file <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/non_coding_genes.xlsx"

# Read non-coding genes from Excel file
non_coding_genes_df <- read.xlsx(non_coding_genes_file, sheet = 1)

# Convert to a vector
non_coding_genes <- as.vector(non_coding_genes_df$non_coding_genes)  # Ensure column name matches saved file

# Confirm loading
cat("Total non-coding genes loaded:", length(non_coding_genes), "\n")


## Classifies genes as either "coding" or "non_coding" and stores this in metadata--(method 1)
# Load required library
library(SingleCellExperiment)

# Load the Seurat object
mouse <- readRDS("/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_GSMO_PHGDH_genotypes_02.RDS")

# Create a named vector indicating whether each gene is coding or non-coding
gene_type_vector <- ifelse(rownames(mouse) %in% non_coding_genes, "non_coding", "coding")
names(gene_type_vector) <- rownames(mouse)  # Ensure names match gene names


# Convert Seurat object to SingleCellExperiment (SCE)
sce <- as.SingleCellExperiment(mouse)

# Add gene metadata
rowData(sce)$gene_type <- gene_type_vector[rownames(sce)]  # Assign to rowData

# Convert back to Seurat object
mouse <- as.Seurat(sce)

# Save the updated Seurat object
saveRDS(mouse, "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_with_gene_type.RDS")

# Load the Seurat object
mouse <- readRDS("/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_with_gene_type.RDS")
str(mouse)
head(mouse@meta.data)

table(mouse@assays$RNA@meta.features$gene_type)
head(mouse@assays$RNA@meta.features)

## Verify that row names match features
all(rownames(mouse@assays$RNA@meta.features) %in% rownames(mouse))



## Classifies genes as either "coding" or "non_coding" and stores this in metadata--(method 2)
# Load the Seurat object
mouse <- readRDS("/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_GSMO_PHGDH_genotypes_02.RDS")

# Define coding and non-coding gene lists
non_coding_genes_file <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/non_coding_genes.xlsx"
non_coding_genes_df <- read.xlsx(non_coding_genes_file, sheet = 1)
non_coding_genes <- as.vector(non_coding_genes_df$non_coding_genes)

# Create a named vector to classify genes
gene_type_vector <- ifelse(rownames(mouse) %in% non_coding_genes, "non_coding", "coding")

# Convert to a data frame with correct row names
gene_metadata <- data.frame(gene_type = gene_type_vector, row.names = rownames(mouse))

# Assign to Seurat v5 meta.data inside the RNA assay
mouse@assays$RNA@meta.data <- gene_metadata  # Correct placement in Seurat v5

# Save the updated Seurat object
saveRDS(mouse, "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_with_gene_type_02.RDS")

# Check first few entries in meta.data inside RNA assay
head(mouse@assays$RNA@meta.data)

# Confirm total counts of coding vs non-coding genes
table(mouse@assays$RNA@meta.data$gene_type)

# Verify that row names match features
all(rownames(mouse@assays$RNA@meta.data) %in% rownames(mouse))

all(rownames(mouse@assays$RNA@meta.data) %in% rownames(mouse@assays$RNA))



### run PCA, cluster analysis, Differential Gene Expression Analysis
# Load necessary libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(stringr)
library(future)
library(data.table)
library(openxlsx)

# Define paths
seurat_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_with_gene_type_02.RDS"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/"

# Define SaveFigure function
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  
  # Save as PNG
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

## Load Seurat Object & Check Metadata
# Load the Seurat object
mouse <- readRDS(seurat_path)

# Check object structure
str(mouse)

# Check metadata
head(mouse@meta.data)

# Check available assays
Assays(mouse)

# Verify gene names
head(rownames(mouse))

# Check total unique samples
table(mouse$sample)

# Check RNA counts and features distribution
summary(mouse$nFeature_RNA)
summary(mouse$nCount_RNA)

# Check treatment groups
table(mouse$treatment)

# Identify highly variable features
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)

# Confirm variable features
length(VariableFeatures(mouse))  # Should return 2000


## Perform Dimensional Reduction (PCA)
# Run PCA
mouse <- RunPCA(mouse)

# Elbow plot to determine optimal PCs
plot <- ElbowPlot(mouse, ndims = 50)
SaveFigure(plot, "PC_elbow_plot", width = 8, height = 10)


## Clustering Analysis
# Find neighbors and clusters
mouse <- FindNeighbors(mouse, dims = 1:30)
mouse <- FindClusters(mouse, resolution = 1)


## Non-Linear Dimensional Reduction (UMAP & tSNE)
# Run UMAP
mouse <- RunUMAP(mouse, dims = 1:30)

# UMAP Plots
plot <- DimPlot(mouse, reduction = "umap", label = TRUE, shuffle = TRUE) + NoAxes()
SaveFigure(plot, "umap_louvain_res_p3", width = 9, height = 8)

# UMAP by sample
plot_sample <- DimPlot(mouse, reduction = "umap", group.by = "sample", shuffle = TRUE) + NoAxes()
SaveFigure(plot_sample, "umap_by_sample", width = 9, height = 8)

# UMAP by treatment
plot_treatment <- DimPlot(mouse, reduction = "umap", group.by = "treatment", shuffle = TRUE) + NoAxes()
SaveFigure(plot_treatment, "umap_by_treatment", width = 9, height = 8)

# Run tSNE
mouse <- RunTSNE(mouse, dims = 1:30)

# tSNE Plots
plot <- DimPlot(mouse, reduction = "tsne", label = TRUE, shuffle = TRUE) + NoAxes()
SaveFigure(plot, "tsne_louvain_res_p3", width = 9, height = 8)

# tSNE by sample
plot_sample <- DimPlot(mouse, reduction = "tsne", group.by = "sample", shuffle = TRUE) + NoAxes()
SaveFigure(plot_sample, "tsne_by_sample", width = 9, height = 8)

# tSNE by treatment
plot_treatment <- DimPlot(mouse, reduction = "tsne", group.by = "treatment", shuffle = TRUE) + NoAxes()
SaveFigure(plot_treatment, "tsne_by_treatment", width = 9, height = 8)

table(mouse$seurat_clusters)

# Save the processed Seurat object
saveRDS(mouse, file = paste0(data_path, "seurat_obj_mouse_with_gene_type_03.RDS"))

##  Subset Clusters
# Remove Cluster 7 and retain selected clusters
mouse.subset <- subset(mouse, subset = seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24))

table(mouse.subset$seurat_clusters)

# Perform PCA on subset
mouse.subset <- RunPCA(mouse.subset)

# Find neighbors and clusters
mouse.subset <- FindNeighbors(mouse.subset, dims = 1:30)
mouse.subset <- FindClusters(mouse.subset, resolution = 1)

table(mouse.subset$seurat_clusters)

# Run UMAP
mouse.subset <- RunUMAP(mouse.subset, dims = 1:30)

# UMAP Plots
plot <- DimPlot(mouse.subset, reduction = "umap", label = TRUE, shuffle = TRUE) + NoAxes()
SaveFigure(plot, "umap_louvain_res_p3_subset", width = 9, height = 8)

# UMAP by sample
plot_sample <- DimPlot(mouse.subset, reduction = "umap", group.by = "sample", shuffle = TRUE) + NoAxes()
SaveFigure(plot_sample, "umap_by_sample_subset", width = 9, height = 8)

# UMAP by treatment
plot_treatment <- DimPlot(mouse.subset, reduction = "umap", group.by = "treatment", shuffle = TRUE) + NoAxes()
SaveFigure(plot_treatment, "umap_by_treatment_subset", width = 9, height = 8)

# Run tSNE on subset
mouse.subset <- RunTSNE(mouse.subset, dims = 1:30)

# tSNE Plots
plot <- DimPlot(mouse.subset, reduction = "tsne", label = TRUE, shuffle = TRUE) + NoAxes()
SaveFigure(plot, "tsne_louvain_res_p3_subset", width = 9, height = 8)

# tSNE by sample
plot_sample <- DimPlot(mouse.subset, reduction = "tsne", group.by = "sample", shuffle = TRUE) + NoAxes()
SaveFigure(plot_sample, "tsne_by_sample_subset", width = 9, height = 8)

# tSNE by treatment
plot_treatment <- DimPlot(mouse.subset, reduction = "tsne", group.by = "treatment", shuffle = TRUE) + NoAxes()
SaveFigure(plot_treatment, "tsne_by_treatment_subset", width = 9, height = 8)

# Redefine object to exclude cluster 7
mouse <- mouse.subset


## Differential Gene Expression Analysis
# Find differentially expressed genes in each cluster
mouse_markers <- FindAllMarkers(mouse, only.pos = TRUE)

# View top markers for each cluster
mouse_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Save DEGs
write.csv(mouse_markers, "mouse_cluster_markers_no7.csv", sep="\t", quote=F, col.names=T, row.names = FALSE)

# Violin plot for key genes
plot <- VlnPlot(mouse, features = c("Mki67", "Meg3", "Aqp4", "Dcn", "Aif1", "Sox10"), group.by = "seurat_clusters")
SaveFigure(plot, "vln_exp1", width = 16, height = 8)

# Violin plot with raw counts
plot <- VlnPlot(mouse, features = c("Mki67", "Meg3", "Aqp4", "Dcn", "Aif1", "Sox10"), layer = "counts", log = TRUE, group.by = "seurat_clusters")
SaveFigure(plot, "vln_exp2", width = 16, height = 8)

# Save and export mouse_markers
saveRDS(mouse_markers, file = paste0(data_path, "seurat_mouse_markers.RDS"))
write.csv(mouse_markers, file = paste0(data_path, "mouse_cluster_markers_no7.csv"), 
          sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# DotPlot for top 5 genes per cluster
top5 <- mouse_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)
plot <- DotPlot(mouse, features = to_plot, group.by = "seurat_clusters") + coord_flip()
SaveFigure(plot, "dplot_top5", width = 9, height = 20)

# Save the processed Seurat object
saveRDS(mouse, file = paste0(data_path, "seurat_obj_mouse_with_gene_type_04.RDS"))



### annotation seurat_obj_mouse_with_gene_type_03 (27 clusters)
# Load the Seurat object
mouse <- readRDS("/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/seurat_obj_mouse_with_gene_type_03.RDS")

# Annotate clusters
mouse.annotated <- mouse
mouse.annotated <- RenameIdents(mouse.annotated, 
                              "0" = "quiescent tumor",
                              "1" = "S+G2 tumor",
                              "2" = "quiescent tumor",
                              "3" = "quiescent tumor",
                              "4" = "mitotic tumor",
                              "5" = "early differentiation",
                              "6" = "mitotic tumor",
                              "8" = "S+G2 tumor",
                              "9" = "quiescent tumor",
                              "10" = "late differentiation",
                              "11" = "early differentiation",
                              "12" = "mitotic tumor",
                              "13" = "astrocytes",
                              "14" = "late differentiation",
                              "15" = "endothelial cells",
                              "16" = "S+G2 tumor",
                              "17" = "oligodendrocyte progenitors",
                              "18" = "fibroblasts",
                              "19" = "microglia",
                              "20" = "microglia Mrc+",
                              "21" = "pericytes",
                              "22" = "oligodendrocytes myelinated",
                              "23" = "mature neurons",
                              "24" = "ependymal")

# Save annotated cluster plots
plot <- DimPlot(mouse.annotated, reduction = "umap", shuffle = TRUE, label = TRUE) + NoAxes()
SaveFigure(plot, "umap_louvain_res_p3_subset_annotated", width = 12, height = 8)

plot <- DimPlot(mouse.annotated, reduction = "tsne", shuffle = TRUE, label = TRUE) + NoAxes()
SaveFigure(plot, "tsne_louvain_res_p3_subset_annotated", width = 12, height = 8)

# Save annotated Seurat object
saveRDS(mouse.annotated, file = paste0(data_path, "seurat_obj_mouse_annotated.RDS"))


### Plot UMAP with coding genes and non_coding genes expression
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Define paths
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"

# Define SaveFigure function
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

# Load the annotated Seurat object
mouse.annotated <- readRDS(paste0(data_path, "seurat_obj_mouse_annotated.RDS"))

table(mouse.annotated$seurat_clusters)

# a table showing the number of cells per annotated cluster
table(Idents(mouse.annotated))

head(mouse.annotated@meta.data)

str(mouse.annotated)

# Check gene type distribution
table(mouse.annotated@assays$RNA@meta.data$gene_type)


# Extract gene names directly from meta.data
all_genes <- rownames(gene_type_info)

head(all_genes)

# Get coding and non-coding genes
coding_genes <- rownames(gene_type_info)[gene_type_info$gene_type == "coding"]
noncoding_genes <- rownames(gene_type_info)[gene_type_info$gene_type == "non_coding"]

# Check counts
cat("Matched Coding Genes:", length(coding_genes), "\n")
cat("Matched Non-Coding Genes:", length(noncoding_genes), "\n")


# Extract expression data properly
expr_data <- GetAssayData(mouse.annotated, assay = "RNA", slot = "data")

# Check if row names exist
rownames(expr_data)[1:10]  # Display first 10 gene names

# Filter coding and non-coding genes from available gene names
all_genes <- rownames(expr_data)
coding_genes <- intersect(all_genes, rownames(mouse.annotated@assays$RNA@meta.data)[mouse.annotated@assays$RNA@meta.data$gene_type == "coding"])
noncoding_genes <- intersect(all_genes, rownames(mouse.annotated@assays$RNA@meta.data)[mouse.annotated@assays$RNA@meta.data$gene_type == "non_coding"])

# Compute per-cell average expression
coding_score <- Matrix::colMeans(expr_data[coding_genes, , drop = FALSE])
noncoding_score <- Matrix::colMeans(expr_data[noncoding_genes, , drop = FALSE])

# Store in Seurat metadata
mouse.annotated$Coding_Score <- coding_score
mouse.annotated$NonCoding_Score <- noncoding_score

# Verify metadata update
head(mouse.annotated@meta.data)


# Classify cells based on whether they express more coding or non-coding genes
mouse.annotated$Gene_Type <- ifelse(mouse.annotated$Coding_Score > mouse.annotated$NonCoding_Score, "Coding-High", "NonCoding-High")

# Check distribution
table(mouse.annotated$Gene_Type)

# Generate UMAP plot
plot <- DimPlot(mouse.annotated, group.by = "Gene_Type", cols = c("blue", "red")) + NoAxes()
SaveFigure(plot, "umap_coding_vs_noncoding", width = 12, height = 8)

# Generate the plot
plot <- FeaturePlot(mouse.annotated, features = c("Coding_Score", "NonCoding_Score"), blend = TRUE)

# Now save it using SaveFigure
SaveFigure(plot, "feature_plot_coding_vs_noncoding_02", width = 12, height = 8)

# Plot UMAP with coding gene expression
plot_coding <- FeaturePlot(mouse.annotated, features = "Coding_Score", cols = c("lightgray", "blue"))
SaveFigure(plot_coding, "umap_coding_score", width = 12, height = 8)

# Plot UMAP with non-coding gene expression
plot_noncoding <- FeaturePlot(mouse.annotated, features = "NonCoding_Score", cols = c("lightgray", "red"))
SaveFigure(plot_noncoding, "umap_noncoding_score", width = 12, height = 8)

# Load required library
library(patchwork)

# Plot UMAP with coding gene expression
plot_coding <- FeaturePlot(mouse.annotated, features = "Coding_Score", cols = c("lightgray", "blue")) + ggtitle("Coding Gene Expression")

# Plot UMAP with non-coding gene expression
plot_noncoding <- FeaturePlot(mouse.annotated, features = "NonCoding_Score", cols = c("lightgray", "red")) + ggtitle("Non-Coding Gene Expression")

# Combine the plots side by side
combined_plot <- plot_coding | plot_noncoding

# Save the combined plot
SaveFigure(combined_plot, "umap_coding_vs_noncoding_combined", width = 18, height = 8)




# ===========================================
# Approach 1: Subset Seurat Object for Coding and Non-Coding Genes
# ===========================================
# This approach subsets the Seurat object separately for coding and non-coding genes.
# It runs PCA, clustering, and UMAP separately on each subset to compare their differences.

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Define paths
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"

# Define SaveFigure function
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}


# Load the annotated Seurat object
mouse.annotated <- readRDS(paste0(data_path, "seurat_obj_mouse_annotated.RDS"))

# Extract coding and non-coding gene names
gene_type_info <- mouse.annotated@assays$RNA@meta.data
coding_genes <- rownames(gene_type_info)[gene_type_info$gene_type == "coding"]
noncoding_genes <- rownames(gene_type_info)[gene_type_info$gene_type == "non_coding"]

# Ensure genes exist in the dataset
coding_genes <- intersect(coding_genes, rownames(mouse.annotated))
noncoding_genes <- intersect(noncoding_genes, rownames(mouse.annotated))

# Subset Seurat object for coding genes only
mouse_coding <- subset(mouse.annotated, features = coding_genes)
mouse_coding <- NormalizeData(mouse_coding)  # Ensure normalization before PCA
mouse_coding <- ScaleData(mouse_coding)
mouse_coding <- RunPCA(mouse_coding, verbose = FALSE)
mouse_coding <- FindNeighbors(mouse_coding, dims = 1:30)
mouse_coding <- FindClusters(mouse_coding, resolution = 1)
mouse_coding <- RunUMAP(mouse_coding, dims = 1:30)

# Subset Seurat object for non-coding genes only
mouse_noncoding <- subset(mouse.annotated, features = noncoding_genes)
mouse_noncoding <- NormalizeData(mouse_noncoding)
mouse_noncoding <- ScaleData(mouse_noncoding)
mouse_noncoding <- RunPCA(mouse_noncoding, verbose = FALSE)
mouse_noncoding <- FindNeighbors(mouse_noncoding, dims = 1:10)
mouse_noncoding <- FindClusters(mouse_noncoding, resolution = 1)
mouse_noncoding <- RunUMAP(mouse_noncoding, dims = 1:10)

# Generate UMAP plots for both coding and non-coding gene clustering
plot_coding <- DimPlot(mouse_coding, reduction = "umap", label = TRUE) + ggtitle("UMAP - Coding Genes")
plot_noncoding <- DimPlot(mouse_noncoding, reduction = "umap", label = TRUE) + ggtitle("UMAP - Non-Coding Genes")

# Combine UMAP plots side by side and save the figure
library(patchwork)
combined_plot <- plot_coding | plot_noncoding
SaveFigure(combined_plot, "umap_coding_vs_noncoding_subset", width = 18, height = 8)

print(combined_plot)

# ===========================================
# Approach 2: Overlay Coding vs. Non-Coding Gene Expression on Annotated Seurat Object
# ===========================================
# This approach retains existing clustering and cell type annotations.
# It overlays coding and non-coding expression levels on the existing UMAP.

# Compute mean expression per cell for coding and non-coding genes
mouse.annotated$Coding_Score <- colMeans(as.matrix(GetAssayData(mouse.annotated, slot = "data")[coding_genes, , drop = FALSE]))
mouse.annotated$NonCoding_Score <- colMeans(as.matrix(GetAssayData(mouse.annotated, slot = "data")[noncoding_genes, , drop = FALSE]))

# Generate UMAP plots colored by coding vs. non-coding gene expression
plot_coding <- FeaturePlot(mouse.annotated, features = "Coding_Score", cols = c("lightgray", "blue")) + ggtitle("Coding Gene Expression")
plot_noncoding <- FeaturePlot(mouse.annotated, features = "NonCoding_Score", cols = c("lightgray", "red")) + ggtitle("Non-Coding Gene Expression")

# Combine the two plots and save
combined_plot <- plot_coding | plot_noncoding
SaveFigure(combined_plot, "umap_coding_vs_noncoding_overlay", width = 18, height = 8)
print(combined_plot)




### coding genes and non_coding genes top 10 most variable genes, 
# cluster umap, Plot UMAP with coding genes and non_coding genes expression,
# Differential Gene Expression Analysis for Coding and Non-Coding Genes 
# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Define paths
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"

# Define SaveFigure function
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

# ===========================================
# Load the Annotated Seurat Object
# ===========================================
mouse.annotated <- readRDS(paste0(data_path, "seurat_obj_mouse_annotated.RDS"))

# Extract coding and non-coding gene names
gene_type_info <- mouse.annotated@assays$RNA@meta.data
coding_genes <- rownames(gene_type_info)[gene_type_info$gene_type == "coding"]
noncoding_genes <- rownames(gene_type_info)[gene_type_info$gene_type == "non_coding"]

# Ensure genes exist in the dataset
coding_genes <- intersect(coding_genes, rownames(mouse.annotated))
noncoding_genes <- intersect(noncoding_genes, rownames(mouse.annotated))

# ===========================================
# Approach 1: Subset Seurat Object for Coding and Non-Coding Genes
# ===========================================
# This approach runs PCA, clustering, and UMAP separately on each subset to compare their differences.

# ---------- Subset and Process Coding Genes ----------
mouse_coding <- subset(mouse.annotated, features = coding_genes)
mouse_coding <- NormalizeData(mouse_coding)
mouse_coding <- ScaleData(mouse_coding)
mouse_coding <- FindVariableFeatures(mouse_coding, selection.method = "vst", nfeatures = 2000)

# Identify top 10 most variable genes
top10_coding <- head(VariableFeatures(mouse_coding), 10)

# Plot variable features for coding genes
plot1 <- VariableFeaturePlot(mouse_coding)
plot2 <- LabelPoints(plot = plot1, points = top10_coding, xnudge = 0, ynudge = 0, repel = TRUE)
SaveFigure((plot1 + plot2), "var_features_coding", width = 12, height = 6)

# Perform dimensional reduction and clustering
mouse_coding <- RunPCA(mouse_coding, verbose = FALSE)
mouse_coding <- FindNeighbors(mouse_coding, dims = 1:30)
mouse_coding <- FindClusters(mouse_coding, resolution = 1)
mouse_coding <- RunUMAP(mouse_coding, dims = 1:30)

# ---------- Subset and Process Non-Coding Genes ----------
mouse_noncoding <- subset(mouse.annotated, features = noncoding_genes)
mouse_noncoding <- NormalizeData(mouse_noncoding)
mouse_noncoding <- ScaleData(mouse_noncoding)
mouse_noncoding <- FindVariableFeatures(mouse_noncoding, selection.method = "vst", nfeatures = 2000)

# Identify top 10 most variable genes
top10_noncoding <- head(VariableFeatures(mouse_noncoding), 10)

# Plot variable features for non-coding genes
plot1 <- VariableFeaturePlot(mouse_noncoding)
plot2 <- LabelPoints(plot = plot1, points = top10_noncoding, xnudge = 0, ynudge = 0, repel = TRUE)
SaveFigure((plot1 + plot2), "var_features_noncoding", width = 12, height = 6)

# Perform dimensional reduction and clustering
mouse_noncoding <- RunPCA(mouse_noncoding, verbose = FALSE)
mouse_noncoding <- FindNeighbors(mouse_noncoding, dims = 1:10)
mouse_noncoding <- FindClusters(mouse_noncoding, resolution = 1)
mouse_noncoding <- RunUMAP(mouse_noncoding, dims = 1:10)

# ---------- Generate UMAP Plots for Both Subsets ----------
plot_coding <- DimPlot(mouse_coding, reduction = "umap", label = TRUE) + ggtitle("UMAP - Coding Genes")
plot_noncoding <- DimPlot(mouse_noncoding, reduction = "umap", label = TRUE) + ggtitle("UMAP - Non-Coding Genes")

# Combine UMAP plots side by side and save the figure
combined_plot <- plot_coding | plot_noncoding
SaveFigure(combined_plot, "umap_coding_vs_noncoding_subset", width = 18, height = 8)

# ===========================================
# Approach 2: Overlay Coding vs. Non-Coding Gene Expression on Annotated Seurat Object
# ===========================================
# This approach overlays coding and non-coding expression levels on the existing UMAP.


# Compute mean expression per cell for coding and non-coding genes
mouse.annotated$Coding_Score <- colMeans(as.matrix(GetAssayData(mouse.annotated, slot = "data")[coding_genes, , drop = FALSE]))
mouse.annotated$NonCoding_Score <- colMeans(as.matrix(GetAssayData(mouse.annotated, slot = "data")[noncoding_genes, , drop = FALSE]))

# Generate UMAP plots colored by coding vs. non-coding gene expression
plot_coding <- FeaturePlot(mouse.annotated, features = "Coding_Score", cols = c("lightgray", "blue")) + ggtitle("Coding Gene Expression")
plot_noncoding <- FeaturePlot(mouse.annotated, features = "NonCoding_Score", cols = c("lightgray", "red")) + ggtitle("Non-Coding Gene Expression")

# Combine the two plots and save
combined_plot <- plot_coding | plot_noncoding
SaveFigure(combined_plot, "umap_coding_vs_noncoding_overlay", width = 18, height = 8)
print(combined_plot)


# ===========================================
# Differential Gene Expression Analysis for Coding and Non-Coding Genes
# ===========================================
# Find differentially expressed genes in each cluster for coding and non-coding genes separately

# ---------- DE Analysis for Coding Subset ----------
coding_markers <- FindAllMarkers(mouse_coding, only.pos = TRUE)

# Save DEGs for coding genes
saveRDS(coding_markers, file = paste0(data_path, "seurat_mouse_coding_markers.RDS"))
write.csv(coding_markers, file = paste0(data_path, "mouse_cluster_markers_coding.csv"), 
          sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# DotPlot for top 5 genes per cluster
top5_coding <- coding_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot_coding <- unique(top5_coding$gene)
plot <- DotPlot(mouse_coding, features = to_plot_coding, group.by = "seurat_clusters") + coord_flip()
SaveFigure(plot, "dotplot_top5_coding", width = 9, height = 20)

# ---------- DE Analysis for Non-Coding Subset ----------
noncoding_markers <- FindAllMarkers(mouse_noncoding, only.pos = TRUE)

# Save DEGs for non-coding genes
saveRDS(noncoding_markers, file = paste0(data_path, "seurat_mouse_noncoding_markers.RDS"))
write.csv(noncoding_markers, file = paste0(data_path, "mouse_cluster_markers_noncoding.csv"), 
          sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# DotPlot for top 5 genes per cluster
top5_noncoding <- noncoding_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot_noncoding <- unique(top5_noncoding$gene)
plot <- DotPlot(mouse_noncoding, features = to_plot_noncoding, group.by = "seurat_clusters") + coord_flip()
SaveFigure(plot, "dotplot_top5_noncoding", width = 9, height = 20)

# ===========================================
# Save Processed Seurat Objects
# ===========================================
saveRDS(mouse_coding, file = paste0(data_path, "seurat_obj_mouse_coding.RDS"))
saveRDS(mouse_noncoding, file = paste0(data_path, "seurat_obj_mouse_noncoding.RDS"))




# ===========================================
# Seurat Differential Gene Expression Analysis, "seurat_obj_mouse_with_gene_type_03"
# ===========================================
# This script performs differential gene expression analysis using a Seurat object.
# It identifies differentially expressed genes, generates violin and dot plots,
# and saves the processed results.

# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ===========================================
# Define Paths and Helper Functions
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/PHGDH/phgdh/"

# Function to Save Figures
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

# ===========================================
# Load the Seurat Object
# ===========================================
seurat_filename <- "seurat_obj_mouse_with_gene_type_03"
mouse <- readRDS(paste0(data_path, seurat_filename, ".RDS"))

# ===========================================
# Differential Gene Expression Analysis
# ===========================================
# Identify differentially expressed genes in each cluster
mouse_markers <- FindAllMarkers(mouse, only.pos = TRUE)

# View top markers for each cluster
mouse_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Save DEGs to file
saveRDS(mouse_markers, file = paste0(data_path, seurat_filename, "_markers.RDS"))
write.csv(mouse_markers, file = paste0(data_path, seurat_filename, "_cluster_markers.csv"), 
          sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# ===========================================
# Visualization of Key Marker Genes
# ===========================================
# Selected marker genes for visualization
selected_features <- c("Aqp4", "Pecam1", "Sox10", "Mrc1", "Dcn", "Meg3", "Ptprc", "Mki67", "Ccnd1", "Nhlh1", "Gabra2")

# ---------- Violin Plot for Key Genes ----------
plot <- VlnPlot(mouse, features = selected_features, group.by = "seurat_clusters", pt.size = 0, ncol = 3)
SaveFigure(plot, paste0(seurat_filename, "_violin_exp1"), width = 20, height = 8)

# Violin plot with raw counts (log scale)
plot <- VlnPlot(mouse, features = selected_features, layer = "counts", log = TRUE, group.by = "seurat_clusters", pt.size = 0, ncol = 3)
SaveFigure(plot, paste0(seurat_filename, "_violin_exp2"), width = 20, height = 8)


# ---------- Feature Plots for Selected Genes ----------
feature_plot <- FeaturePlot(mouse, features = selected_features, ncol = 4)
SaveFigure(feature_plot, paste0(seurat_filename, "_feature_exp"), width = 20, height =12 )


# Count the number of cells in each cluster
table(mouse$seurat_clusters)



# ---------- DotPlot for Top 5 Genes per Cluster ----------
# Identify top 5 differentially expressed genes per cluster
top5 <- mouse_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)

# Generate and save dot plot
plot <- DotPlot(mouse, features = to_plot, group.by = "seurat_clusters") + coord_flip()
SaveFigure(plot, paste0(seurat_filename, "_dotplot_top5"), width = 9, height = 20)


