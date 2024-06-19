# Load necessary libraries
.libPaths(c("/data/Common_Folder/R/x86_64-pc-linux-gnu-library/4.1/", .libPaths()))
library(lifecycle)
library(dplyr)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(irlba)
library(ggplot2)
library(data.table)
library(Seurat)
library(future)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load the RNA and ATAC data
counts <- Read10X("/data/Common_Folder/Multiome_aim1c_data/filtered_featues_bc_matrices/Sample_1_filtered_feature_bc_matrix/filtered_feature_bc_matrix/")

# Extract RNA and ATAC data
rna_counts <- counts$`Gene Expression`
atac_counts <- counts$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
genome(annotations) <- "hg38"
frag.file <- "/data/Common_Folder/Multiome_aim1c_data/atac_fragments/Sample_14_atac_fragments/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

# Pre-QC plot
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

# QC
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# Post-QC plot
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Metadata assignment
pbmc@meta.data$Sample <- "Sample14"
pbmc@meta.data$Cond <- "GBA"
pbmc@meta.data$Clone <- "PL1C7"

# Save preprocessed object
saveRDS(pbmc, file = "/data/dehestani/multiome_aim1c/preprocessed_RNA_ATAC_objects_before_merge/s14.rds")

# Read individual samples
s1 <- readRDS("s1.rds")
s2 <- readRDS("s2.rds")
s3 <- readRDS("s3.rds")
s4 <- readRDS("s4.rds")
s5 <- readRDS("s5.rds")
s6 <- readRDS("s6.rds")
s7 <- readRDS("s7.rds")
s8 <- readRDS("s8.rds")
s9 <- readRDS("s9.rds")
s10 <- readRDS("s10.rds")
s11 <- readRDS("s11.rds")
s12 <- readRDS("s12.rds")
s13 <- readRDS("s13.rds")
s14 <- readRDS("s14.rds")

# Extract the RNA assay counts
rna_counts <- s14[["RNA"]]

# Create a new Seurat object with only the RNA counts
rna_seurat_obj <- CreateSeuratObject(counts = rna_counts)

# Transfer metadata if needed
rna_seurat_obj@meta.data <- s14@meta.data

saveRDS(rna_seurat_obj, file = "/data/dehestani/multiome_aim1c/Only_RNA_analysis/s14.rds")

# Merge samples
merged <- merge(s1, y = c(s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14), merge.data = TRUE, project = "All")
saveRDS(merged, file = "/data/dehestani/multiome_aim1c/Only_RNA_analysis/objects_singular_and_merged/merged_RNA_preprocessed.rds")

# Load merged object
pbmc <- readRDS("/data/dehestani/multiome_aim1c/Only_RNA_analysis/objects_singular_and_merged/merged_RNA_preprocessed.rds")

# Harmony for integration
DefaultAssay(pbmc) <- "RNA"
pbmc <- SplitObject(pbmc, split.by = "Sample")
for (i in 1:length(pbmc)) {
  pbmc[[i]] <- SCTransform(pbmc[[i]], verbose = T, vars.to.regress = "percent.mt")
}

features <- SelectIntegrationFeatures(object.list = pbmc, nfeatures = 1000)
options(future.globals.maxSize = 100000 * 1024^2)
pbmc <- PrepSCTIntegration(object.list = pbmc, anchor.features = features, verbose = FALSE)

pbmc14 <- pbmc$Sample14
pbmc$pbmc14 <- NULL
pbmc <- merge(pbmc14, y = c(pbmc), project = "pd", merge.data = TRUE)
pbmc <- RunPCA(object = pbmc, npcs = 20, assay = "SCT", features = features)
DimHeatmap(pbmc, dims = 1:15, cells = 1000, balanced = TRUE)

pbmc <- RunHarmony(
  object = pbmc,
  group.by.vars = 'Sample',
  reduction = 'pca',
  assay.use = 'SCT',
  project.dim = FALSE,
  reduction.save = "harmony_gex"
)

pbmc <- RunUMAP(object = pbmc, assay = "SCT", reduction = "harmony_gex", dims = 1:20)
pbmc <- FindNeighbors(object = pbmc, assay = "SCT", reduction = "harmony_gex", graph.name = "test", dims = 1:20)
pbmc <- FindClusters(object = pbmc, graph.name = "test", resolution = 0.2)
pbmc <- FindSubCluster(
  pbmc,
  '3',
  graph.name = "test",
  subcluster.name = "sub.cluster",
  resolution = 0.1,
  algorithm = 1
)

DimPlot(pbmc)
DimPlot(pbmc, reduction = "umap", group.by = "sub.cluster")
saveRDS(pbmc, file = "/data/dehestani/multiome_aim1c/Only_RNA_analysis/processed_objects/Aim1c_full_processed_res01_harmony_lognorm.rds")

# Cluster identity
DefaultAssay(pbmc) <- "SCT"
feature1 = c("RFX4", "HES1", "SLIT2") # Early-neuron Progenitor
feature2 = c("DLK1", "LGALS1", "VCAN") # Late-neuron Progenitor
feature3 = c("TH", "ATP1A3", "ZCCHC12", "MAP2", "SYT1") # Dopaminergic Neurons
feature4 = c("TPH1", "SLC18A1", "SLC18A2", "SNAP25") # Immature Dopaminergic Neurons
feature5 = c("HMGB2", "TOP2A", "MKI67") # Proliferating Floor Plate Progenitors
feature6 = c("KRT19", "KRT8", "COL17A1") # Neuroepithelial-like Cells
feature7 = c("MLF1", "STOML3", "FOXJ1") # Ependymal-like Cells

FeaturePlot(pbmc, features = feature1, reduction = "umap")
DotPlot(pbmc, features = feature1) + RotatedAxis()

# Set the active identity to 'sub.cluster'
Idents(pbmc) <- "sub.cluster"

# Rename the identities
pbmc <- RenameIdents(pbmc, `0` = "Early-neuron Progenitor", `1` = "Late-neuron Progenitor", `2` = "Dopaminergic Neurons", `3` = "Immature Dopaminergic Neurons", `4` = "Proliferating Floor Plate Progenitors", `5` = "Neuroepithelial-like Cells", `6` = "Ependymal-like Cells")
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# Save the final object
saveRDS(pbmc, file = "/data/dehestani/multiome_aim1c/Only_RNA_analysis/processed_objects/Aim1c_full_processed_final.rds")
