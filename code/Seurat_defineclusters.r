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

#split into male and female data sets
mathysF <- subset(x=mathys2, subset=sex=='female')

mathysF <- SCTransform(mathysF, do.scale = TRUE, verbose = TRUE)
#need to add the metadata to the seurat object:
mathysF <- AddMetaData(mathysF, mathys_metadata)
head(x=mathysF[[]])

mathysF <- RunPCA(mathysF, npcs = 30)
mathysF <- RunUMAP(mathysF, dims = 1:30)

mathysF <- FindNeighbors(mathysF, reduction = "pca", dims = 1:30)
mathysF <- FindClusters(mathysF, verbose = TRUE)
DimPlot(mathysF, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysF, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysF, reduction="umap", label = TRUE, repel=TRUE)








#subset the main microglia cluster (#15)
micro <- subset(mathysF, idents=15)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
dim(micro)

#use Seurat integration to correct for batch effects in this cluster
table(micro$batch)

# normalize and identify variable features for each batch independently
mathys.list <- SplitObject(micro, split.by = "batch")
mathys.list <- lapply(X = mathys.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 15000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mathys.list, nfeatures=15000)

mathys.anchors <- FindIntegrationAnchors(object.list = mathys.list, anchor.features = features)

# this command creates an 'integrated' data assay
mathys.combined <- IntegrateData(anchorset = mathys.anchors, k.weight = 34)

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
DimPlot(mathys.combined, group.by = "ros_ids", reduction="umap", label = TRUE, repel=TRUE)

DimPlot(mathys.combined, reduction="umap", label = TRUE, repel=TRUE)

integrated <- GetAssayData(object = mathys.combined, assay = "integrated", slot = "data")

#export counts & metadata for trajectory analysis
# microcounts <- GetAssayData(object = micro, slot = "counts")
saveRDS(integrated, file='~/scRNAseq-subtype-mapping/Seurat_microcounts.RDS')
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
DimPlot(mathysF, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysF, reduction="umap", label = TRUE, repel=TRUE)
#subset the main oligodendrocyte clusters (#1, 2, 5)
micro <- subset(mathysF, idents=c(1, 2, 5))
dim(micro)
table(mathysF$predicted.subclass_label)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, reduction="umap", label = TRUE, repel=TRUE)
#use Seurat integration to correct for batch effects in this cluster
table(micro$batch)

# normalize and identify variable features for each batch independently
mathys.list <- SplitObject(micro, split.by = "batch")
mathys.list <- lapply(X = mathys.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 15000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mathys.list, nfeatures=15000)

mathys.anchors <- FindIntegrationAnchors(object.list = mathys.list, anchor.features = features)

# this command creates an 'integrated' data assay
mathys.combined <- IntegrateData(anchorset = mathys.anchors, k.weight = 400)

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
DimPlot(mathys.combined, group.by = "ros_ids", reduction="umap", label = TRUE, repel=TRUE)

DimPlot(mathys.combined, reduction="umap", label = TRUE, repel=TRUE)

integrated <- GetAssayData(object = mathys.combined, assay = "integrated", slot = "data")






#try fastMNN
DimPlot(mathysF, reduction="umap", label = TRUE, repel=TRUE)
#subset the main oligodendrocyte clusters (#1, 2, 5)
micro <- subset(mathysF, idents=c(1, 2, 5))

micro <- NormalizeData(micro)
micro <- FindVariableFeatures(micro, nfeatures = 3000)
micro <- RunFastMNN(object.list = SplitObject(micro, split.by = "batch"), features=15000)
micro <- ScaleData(micro, verbose = FALSE)
micro <- RunUMAP(micro, reduction = "mnn", dims = 1:30)
micro <- FindNeighbors(micro, reduction = "mnn", dims = 1:30)
micro <- FindClusters(micro)

DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, group.by = "ros_ids", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, reduction="umap", label = TRUE, repel=TRUE)

mnn <- GetAssayData(object = micro, assay = "mnn.reconstructed", slot = "data")
dim(mnn)

p4 <- synapser::synGet('syn26136476')
mathys_DEG <- read.csv(p4$path, header=TRUE)
mathys_DEG2 <- unique(mathys_DEG$DEGgenes)
names(mathys_DEG2)[names(mathys_DEG2) == "mathys_DEG2"] <- "genes"

genes2 <- c()
for (gene in unique(c(as.vector(mathys_DEG$DEGgenes)))){
  
  
  if (gene %in% rownames(mnn)){
    genes2 <- c(genes2,which(rownames(mnn)==gene))
  }
}
mnn2 <- mnn[genes2,]
dim(mnn2)





###try SCT
DimPlot(mathysF, reduction="umap", label = TRUE, repel=TRUE)
#subset the main oligodendrocyte clusters (#1, 2, 5)
micro <- subset(mathysF, idents=c(1, 2, 5))
dim(micro)
micro.list <- SplitObject(micro, split.by="batch")
micro.list <- lapply(X=micro.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list=micro.list, nfeatures=3000)
micro.list <- PrepSCTIntegration(object.list=micro.list, anchor.features=features)

micro.anchors <- FindIntegrationAnchors(object.list=micro.list, normalization.method="SCT", anchor.features=features)
micro.combined.sct <- IntegrateData(anchorset=micro.anchors, normalization.method="SCT")

micro.combined.sct <- RunPCA(micro.combined.sct)
micro.combined.sct <- RunUMAP(micro.combined.sct, reduction="pca", dims=1:30)


DimPlot(micro.combined.sct, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)

sct <- GetAssayData(object = micro.combined.sct, assay = "SCT", slot = "counts")
sct_counts <- GetAssayData(object = micro.combined.sct, assay = "SCT", slot = "data")
dim(sct)
dim(sct_counts)
sct3 <- GetAssayData(object = micro.combined.sct, assay = "integrated", slot = "data")
dim(sct3)
p4 <- synapser::synGet('syn26136476')
mathys_DEG <- read.csv(p4$path, header=TRUE)
mathys_DEG2 <- unique(mathys_DEG$DEGgenes)
names(mathys_DEG2)[names(mathys_DEG2) == "mathys_DEG2"] <- "genes"

genes2 <- c()
for (gene in unique(c(as.vector(mathys_DEG$DEGgenes)))){
  
  
  if (gene %in% rownames(mnn)){
    genes2 <- c(genes2,which(rownames(mnn)==gene))
  }
}
mnn2 <- mnn[genes2,]
dim(mnn2)







#recluster oligos and run SCT, regress out batch
#need to add the metadata to the seurat object:
DimPlot(mathysF, reduction="umap", label = TRUE, repel=TRUE)
#subset the main oligodendrocyte clusters (#1, 2, 5)
micro2 <- subset(mathysF, idents=c(1, 2, 5))
DimPlot(micro2, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)

micro <- SCTransform(micro, do.scale = TRUE, verbose = TRUE, vars.to.regress="ros_ids")

micro <- RunPCA(micro, npcs = 30)
micro <- RunUMAP(micro, dims = 1:30)

micro <- FindNeighbors(micro, reduction = "pca", dims = 1:30)
micro <- FindClusters(micro, verbose = TRUE)
DimPlot(micro, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, reduction="umap", label = TRUE, repel=TRUE)
DimPlot(micro, group.by = "ros_ids", reduction="umap", label = TRUE, repel=TRUE)

sct_counts <- GetAssayData(object = micro, assay = "SCT", slot = "counts")
sct_data <- GetAssayData(object = micro, assay = "SCT", slot = "counts")

genes2 <- c()
for (gene in unique(c(as.vector(mathys_DEG$DEGgenes)))){
  
  
  if (gene %in% rownames(mnn)){
    genes2 <- c(genes2,which(rownames(mnn)==gene))
  }
}
mnn2 <- mnn[genes2,]
dim(mnn2)









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
micro <- subset(mathys3, idents=0)
dim(micro)
table(mathys3$predicted.subclass_label)
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
