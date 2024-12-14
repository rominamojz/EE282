# Load necessary libraries
library(Seurat)  
library(Matrix) 
library(ggplot2)

# Step 1: Read in data
counts = t(readMM("/share/crsp/lab/seyedam/share/igvf_pipeline/IGVF_analysis/cellbender_tissues/annotated/Romina_counts.mtx")) 
meta = read.csv("/share/crsp/lab/seyedam/share/igvf_pipeline/IGVF_analysis/cellbender_tissues/annotated/Romina_cell_metadata.csv")     
genes = read.csv("/share/crsp/lab/seyedam/share/igvf_pipeline/IGVF_analysis/cellbender_tissues/annotated/Romina_genes.csv")        

# Step 2: Assign row and column names
rownames(counts) = genes$gene_name  
colnames(counts) = meta$cellID     
rownames(meta) = meta$cellID      

# Step 3: Create a Seurat object
obj = CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0, meta.data = meta)
# Creates a Seurat object, combining the count matrix and metadata, with no filtering applied

obj  # View a summary of the Seurat object

# Step 4: Add mitochondrial gene percentage to metadata
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
# Calculates the percentage of mitochondrial gene expression and adds it to metadata

# Step 5: Set cell identities based on tissue type
Idents(obj) = obj$Tissue  # Assigns cell identities to the "Tissue" column in metadata

# Step 6: Visualize quality control metrics with violin plots
plot <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "../../output/std_average.jpg", height = 7, width = 12, plot = plot, quality = 50)
plot  # Shows the violin plot

# Step 7: Scatter plots for QC metric relationships
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")  # RNA count vs mitochondrial percentage
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  # RNA count vs number of features
ggsave(filename = "../../output/rna_percent.jpg", height = 7, width = 12, plot = plot1, quality = 50, bg="white")
ggsave(filename = "../../output/rna_feature.jpg", height = 7, width = 12, plot = plot2, quality = 50, bg="white")
plot1  # Show scatter plot 1
plot2  # Show scatter plot 2

# Step 8: Normalize the data
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Normalizes the expression data by log transformation and scales to 10,000 counts per cell

# Step 9: Identify highly variable features
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
# Identifies 2000 most variable genes using variance-stabilizing transformation (VST)

# Step 10: Visualize the most variable features
top10 <- head(VariableFeatures(obj), 10)  # Get the top 10 highly variable genes
plot1 <- VariableFeaturePlot(obj)  # Plot all variable features
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)  # Annotate top 10 genes
ggsave(filename = "../../output/std_average.jpg", height = 7, width = 12, plot = plot2, quality = 50, bg="white")
plot2  # Show the annotated variable feature plot

# Step 11: Scale the data
all.genes <- rownames(obj)  # Get all genes in the dataset
obj <- ScaleData(obj, features = all.genes)  # Scales and centers gene expression data

# Step 12: Perform PCA
obj <- RunPCA(obj, features = VariableFeatures(object = obj))  # Run PCA on variable genes
print(obj[["pca"]], dims = 1:5, nfeatures = 5)  # Print top 5 genes contributing to first 5 PCs

# Visualize PCA results
vizPlot <- VizDimLoadings(obj, dims = 1:2, reduction = "pca")  # Visualize PCA loadings for dims 1 and 2
ggsave(filename = "../../output/PC1_PC2.jpg", height = 7, width = 12, plot = vizPlot, quality = 50, bg="white")
vizPlot  # Show PCA loadings plot

dimViz <- DimPlot(obj, reduction = "pca") + NoLegend()  # PCA scatter plot
ggsave(filename = "../../output/dim_pc1_pc2.jpg", height = 7, width = 12, plot = dimViz, quality = 50, bg="white")
dimViz  # Show PCA scatter plot

pc_std <- ElbowPlot(obj)  # Elbow plot to determine optimal number of PCs
ggsave(filename = "../../output/pc_std.jpg", height = 7, width = 12, plot = pc_std, quality = 50, bg="white")
pc_std  # Show elbow plot

# Step 13: Clustering
obj <- FindNeighbors(obj, dims = 1:10)  # Find nearest neighbors using first 10 PCs
obj <- FindClusters(obj, resolution = 0.5)  # Cluster cells into groups with resolution 0.5

# Step 14: UMAP dimensionality reduction
obj <- RunUMAP(obj, dims = 1:10)  # Run UMAP on first 10 PCs

# Step 15: Visualize UMAP results
plot <- DimPlot(obj, reduction = "umap")
ggsave(filename = "../../output/umap_ind_cluster.jpg", height = 7, width = 12, plot = plot, quality = 50, bg="white")
plot  # Show UMAP plot

# Additional visualizations grouped by metadata attributes
plot <- DimPlot(obj, reduction = "umap", group.by = 'celltype', label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "../../output/umap_celltype.jpg", height = 7, width = 12, plot = plot, quality = 50, bg="white")
plot

plot <- DimPlot(obj, reduction = "umap", group.by = 'Genotype', pt.size = 0.5)
ggsave(filename = "../../output/umap_genotype.jpg", height = 7, width = 12, plot = plot, quality = 50, bg="white")
plot

plot <- DimPlot(obj, reduction = "umap", group.by = 'Sex', pt.size = 0.5)
ggsave(filename = "../../output/umap_sex.jpg", height = 7, width = 12, plot = plot, quality = 50, bg="white")
plot

# Step 16: Differential expression analysis
Idents(obj) = "Sex"  # Set identities for differential expression analysis
monocyte.de.markers <- FindMarkers(obj, ident.1 = "Male", ident.2 = "Female")
head(monocyte.de.markers)  # View results of differential expression

# Additional heatmaps and violin plots
DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE)  # Heatmap for PCA dimension 1
DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE)  # Heatmap for first 15 dimensions

monocyte.de.markers['Xist',]  # View results for a specific gene

plot <- VlnPlot(obj, features = "Myh4", group.by = "celltype")
ggsave(filename = "../../output/myh4.jpg", height = 7, width = 12, plot = plot, quality = 50, bg="white")
plot

plot <- VlnPlot(obj, features = "Actn3", group.by = "celltype")
ggsave(filename = "../../output/actn3.jpg", height = 7, width = 12, plot = plot, quality = 50, bg="white")
plot

# Save the Seurat object for future use
# saveRDS(obj, file= "data.rds")

# Rename cluster identities for better labeling
new.cluster.ids <- c("Myh4", "Myh2")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

# UMAP visualization with updated labels
plot <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "../../output/umap.jpg", height = 7, width = 12, plot = plot, quality = 50, bg="white")
plot
