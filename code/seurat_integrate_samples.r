install.packages("remotes")
install.packages("R.utils")
remotes::install_github("satijalab/seurat-wrappers")
#install.packages("matrixStats")
remotes::install_github("HenrikBengtsson/matrixStats", ref="develop")
BiocManager::install("batchelor")
BiocManager::install("SingleCellExperiment")
library(Seurat)
library(sctransform)
library(ggplot2)
library(SeuratWrappers)
library(patchwork)

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
mathys2 <- AddMetaData(mathys2, mathys_metadata)
head(x=mathys2[[]])

#split into male and female data sets
mathysF <- subset(x=mathys2, subset=sex=='female')



p4 <- synapser::synGet('syn26136476')
mathys_DEG <- read.csv(p4$path, header=TRUE)
mathys_DEG2 <- unique(mathys_DEG$DEGgenes)
names(mathys_DEG2)[names(mathys_DEG2) == "mathys_DEG2"] <- "genes"



# split the dataset into a list of seurat objects (batches)
mathys.list <- SplitObject(mathysF, split.by = "batch")

# normalize and identify variable features for each batch independently
mathys.list <- lapply(X = mathys.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mathys.list, nfeatures=10000)

mathys.anchors <- FindIntegrationAnchors(object.list = mathys.list, anchor.features = features)

# this command creates an 'integrated' data assay
mathys.combined <- IntegrateData(anchorset = mathys.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(mathys.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
mathys.combined <- ScaleData(mathys.combined)
mathys.combined <- RunPCA(mathys.combined, npcs = 30, verbose = FALSE)
mathys.combined <- RunUMAP(mathys.combined, reduction = "pca", dims = 1:30)
mathys.combined <- FindNeighbors(mathys.combined, reduction = "pca", dims = 1:30)
mathys.combined <- FindClusters(mathys.combined, resolution = 0.5)
DimPlot(mathys.combined, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathys.combined, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)

DimPlot(mathys.combined, reduction="umap", label = TRUE, repel=TRUE)




integrated <- GetAssayData(object = mathys.combined, assay = "integrated", slot = "data")



#subset the main microglia cluster (#16)
micro <- subset(mathys2, idents=16)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
dim(micro)
#relabel all cells as Microglia
micro$refined.subtype_label <- 'Micro'
#export counts & metadata for trajectory analysis
microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(microcounts, file='~/scRNAseq-subtype-mapping/Seurat_microcounts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_microcounts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

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
micro$refined.subtype_label <- 'Astro'
#export for trajectory analysis
microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(microcounts, file='~/scRNAseq-subtype-mapping/Seurat_astrocounts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_astrocounts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_astrometa.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_astrometa.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

#get oligodendrocytes by cluster and overall
DimPlot(mathysF, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysF, reduction="umap", label = TRUE, repel=TRUE)
#subset the main oligodendrocyte clusters (#0, 2, 6)
micro <- subset(mathys3, idents=c(0, 2, 6))
dim(micro)
table(mathys3$predicted.subclass_label)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)

#relabel all cells as Oligo
micro$refined.subtype_label <- 'Oligo'
microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(microcounts, file='~/scRNAseq-subtype-mapping/Seurat_oligocounts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligocounts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
#export metadata for trajectory analysis
Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_oligometa.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligometa.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

#Cluster 0
#subset the main oligodendrocyte clusters (#0)
micro <- subset(mathysF, idents=0)
dim(micro)
table(mathysF$predicted.subclass_label)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)

#relabel all cells as Oligo1
micro$refined.subtype_label <- 'Oligo1'
microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(microcounts, file='~/scRNAseq-subtype-mapping/Seurat_oligocounts1.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligocounts1.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
#export metadata for trajectory analysis
Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_oligometa1.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligometa1.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

#relabel all cells as Oligo2
#subset the main oligodendrocyte clusters (#0)
micro <- subset(mathys3, idents=2)
dim(micro)

DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
micro$refined.subtype_label <- 'Oligo1'
microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(microcounts, file='~/scRNAseq-subtype-mapping/Seurat_oligocounts2.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligocounts2.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
#export metadata for trajectory analysis
Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_oligometa2.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligometa2.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

#relabel all cells as Oligo3
#subset the main oligodendrocyte clusters (#0)
micro <- subset(mathys3, idents=6)
dim(micro)

DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
micro$refined.subtype_label <- 'Oligo3'
microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(microcounts, file='~/scRNAseq-subtype-mapping/Seurat_oligocounts3.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligocounts3.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
#export metadata for trajectory analysis
Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_oligometa3.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_oligometa3.RDS', parentId='syn25871777')
file <- synapser::synStore(file)


#get OPCs
DimPlot(mathys3, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathys3, reduction="umap", label = TRUE, repel=TRUE)
#subset the main OPC clusters (#12)
micro <- subset(mathys3, idents=12)
dim(micro)
table(mathys3$predicted.subclass_label)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)

#relabel all cells as Oligo
micro$refined.subtype_label <- 'OPC'
microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(microcounts, file='~/scRNAseq-subtype-mapping/Seurat_OPCcounts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_OPCcounts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
#export metadata for trajectory analysis
Seurat_micrometa <- micro@meta.data
saveRDS(Seurat_micrometa, file='~/scRNAseq-subtype-mapping/Seurat_OPCmeta.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Seurat_OPCmeta.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
