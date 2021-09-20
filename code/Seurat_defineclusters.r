install.packages("remotes")
install.packages("R.utils")
remotes::install_github("satijalab/seurat-wrappers")
install.packages("matrixStats")
remotes::install_github("HenrikBengtsson/matrixStats", ref="develop")
BiocManager::install("batchelor")
BiocManager::install("SingleCellExperiment")
library(Seurat)
library(sctransform)
library(ggplot2)
library(SeuratWrappers)

synapser::synLogin("lheath", "Q54!A!&9iCfl")


#raw counts matrix, not normalized, with labeled rownames and colnames:
p1 <- synapser::synGet('syn25871862')
counts <- readRDS(p1$path)
counts <- as(counts, "dgCMatrix")

p <- synapser::synGet('syn25871851')
mathys_metadata <- read.csv(p$path)
rownames(mathys_metadata) <- mathys_metadata$TAG
mathys_metadata$X<-NULL


#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
mathys_metadata <- mathys_metadata[order(match(mathys_metadata$TAG, colnames_counts$columnnames)),]

mathys2 <- CreateSeuratObject(counts = counts, project = "Clustering", min.cells = 3, min.features = 200)

head(mathys2@meta.data, 20)

mathys2 <- AddMetaData(mathys2, mathys_metadata)
head(x=mathys2[[]])


mathys2 <- NormalizeData(mathys2)
mathys2 <- FindVariableFeatures(mathys2)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

mathys2 <- RunFastMNN(object.list = SplitObject(mathys2, split.by = "ros_ids"))
mathys2 <- RunUMAP(mathys2, reduction="mnn", dims=1:30)
mathys2 <- FindNeighbors(mathys2, reduction="mnn", dims=1:30)
mathys2 <- FindClusters(mathys2)
DimPlot(mathys2, group.by="predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathys2, reduction="umap", label = TRUE, repel=TRUE)


mathys3 <- CreateSeuratObject(counts = counts, project = "Clustering")

mathys3 <- SCTransform(mathys3, do.scale = TRUE, verbose = TRUE)
#need to add the metadata to the seurat object:
mathys3 <- AddMetaData(mathys3, mathys_metadata)
head(x=mathys3[[]])

mathys3 <- RunPCA(mathys3, npcs = 30)
mathys3 <- RunUMAP(mathys3, dims = 1:30)

mathys3 <- FindNeighbors(mathys3, reduction = "pca", dims = 1:30)
mathys3 <- FindClusters(mathys3, verbose = TRUE)
DimPlot(mathys3, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathys3, reduction="umap", label = TRUE, repel=TRUE)


#subset the main microglia cluster (#16)
micro <- subset(mathys3, idents=16)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
dim(micro)
#relabel all cells as Microglia
micro$refined.subtype_label <- 'Micro'
saveRDS(micro, file='~/scRNAseq-subtype-mapping/Seurat_microcounts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_microcounts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
#export for trajectory analysis
Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_micrometa.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_micrometa.RDS', parentId='syn25871777')
file <- synapser::synStore(file)


#get astrocytes
DimPlot(mathys3, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathys3, reduction="umap", label = TRUE, repel=TRUE)
#subset the main astrocyte cluster (#8)
micro <- subset(mathys3, idents=8)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
dim(micro)
#relabel all cells as Microglia
micro$refined.subtype_label <- 'Micro'
saveRDS(micro, file='~/scRNAseq-subtype-mapping/Seurat_microcounts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_microcounts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
#export for trajectory analysis
Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_micrometa.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_micrometa.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
