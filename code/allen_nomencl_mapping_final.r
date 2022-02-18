
#upload the gene expression matrix from allen institute's M1 reference data:

p1 <- synapser::synGet('syn25871859')
M1ref <- readRDS(p1$path)
#turn the csv file into a sparse matrix and transpose, so cells are columns and genes are rows
#rownames(M1ref) <- M1ref$sample_name
#M1ref$sample_name<-NULL
#M1ref <- as.matrix(M1ref)
#M1ref <- Matrix(M1ref, sparse=TRUE)
#M1ref <- t(M1ref)
saveRDS(M1ref, file="~/scRNAseq-subtype-mapping/data/M1ref_counts.rds")
M1ref <- readRDS(file="~/scRNAseq-subtype-mapping/data/M1ref_counts.rds")

#read in metadata for cells:
p2 <- synapser::synGet('syn25871900')
metadata <- read.csv(p2$path)
rownames(metadata) <- metadata$sample_name
metadata$X<-NULL

#counts <- M1ref@assays$RNA@counts

dat <- CreateSeuratObject(counts = M1ref, project = "Mapping", min.cells = 3, min.features = 200)
head(dat@meta.data, 5)
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#MT genes already removed, no need to worry about those

FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dat <- SCTransform(dat, do.scale = TRUE, verbose = TRUE)

#need to add the metadata to the seurat object:
dat <- AddMetaData(dat, metadata)
head(x=dat[[]])

dat <- RunPCA(dat, npcs = 30, return.model = TRUE, verbose = TRUE)
dat <- RunUMAP(dat, dims = 1:30, return.model = TRUE, verbose = TRUE)

dat <- FindNeighbors(dat, reduction = "pca", dims = 1:30, cache.index = TRUE, verbose = TRUE)
dat <- FindClusters(dat, verbose = TRUE)
DimPlot(dat, group.by = "subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(dat, group.by = "class_label", reduction="umap", label = TRUE, repel=TRUE)


# read in mathys cds
#mathys <- readRDS(file="~/scAD_analysis/Mathys_scrannorm_cds.rds")
# p <- synapser::synGet('syn22363128')
# mathys <- readRDS(p$path)
# dim(mathys)
# class(mathys)

#get the metadata:
#mathys_meta <- pData(mathys)

# non-normalized count matrix for mathys:
p <- synapser::synGet('syn18686381')
counts <- readMM(p$path)
counts <- as(counts, "dgCMatrix")
#get all mathys metadata
#mathys_meta <- readRDS(file="~/scRNAseq-subtype-mapping/data/mathys_metadata.rds")
#p2 <- synapser::synGet('syn24171852')
p2 <- synapser::synGet('syn26530396')
mathys_meta <- read.csv(p2$path)
#get short gene names list and make them into rownames on counts file: filtered_gene_row_names.txt
p3 <- synapser::synGet('syn18686382')

colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
mathys_meta <- mathys_meta[order(match(mathys_meta$TAG, colnames_counts$columnnames)),]

rownames(counts) <- readLines(p3$path)
colnames(counts) <- mathys_meta[,1]

#save the data structures:
#saveRDS(counts, file="~/scRNAseq-subtype-mapping/data/mathys_scrannorm_counts.rds")
#saveRDS(mathys_meta, file="~/scRNAseq-subtype-mapping/data/mathys_metadata.rds")
#saveRDS(dat, file="~/scRNAseq-subtype-mapping/data/M1reference_seurat.RDS")


mathys2 <- CreateSeuratObject(counts = counts, project = "Mapping", min.cells = 3, min.features = 200)
head(mathys2@meta.data, 20)

mathys_meta <- as.data.frame(mathys_meta)
rownames(mathys_meta)<-colnames(mathys2)
mathys2 <- AddMetaData(mathys2, mathys_meta)
head(x=mathys2[[]])

mathys2 <- SCTransform(mathys2, do.scale = TRUE, verbose = TRUE)

### map each individual separately (for greater accuracy in mapping):
mathys2.batches <- SplitObject(mathys2, split.by = "ros_ids")
mathys2.batches <- lapply(X = mathys2.batches, FUN = SCTransform, verbose = TRUE)

anchors <- list()
for (i in 1:length(mathys2.batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = dat,
    query = mathys2.batches[[i]],
    normalization.method = "SCT",
    reference.reduction = "pca", 
    dims = 1:30
  )
}


for (i in 1:length(mathys2.batches)) {
  mathys2.batches[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = mathys2.batches[[i]],
    reference = dat, 
    refdata = list(
      class_label = "class_label",
      subclass_label = "subclass_label",
      cluster_label = "cluster_label"),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}
head(x=mathys2.batches[[1]])


p1 <- DimPlot(mathys2.batches[[1]], reduction = 'ref.umap', group.by = 'predicted.subclass_label', label.size = 3)
p2 <- DimPlot(mathys2.batches[[2]], reduction = 'ref.umap', group.by = 'predicted.subclass_label', label.size = 3)
p1 + p2 + patchwork::plot_layout(guides = "collect")
mathys3 <- merge(mathys2.batches[[1]], mathys2.batches[2:length(mathys2.batches)], merge.dr = "ref.umap")
DimPlot(mathys3, reduction = "ref.umap", group.by =  "predicted.subclass_label", label = FALSE, repel = TRUE, label.size = 1)
head(x=mathys3[[]])

hist(mathys3$predicted.subclass_label.score)
hist(mathys3$predicted.subclass_label.score)



#pull the metadata with the new predictions from the Seurat object:
metadata_update <- mathys3@meta.data


##need to save in synapse
write.csv(metadata_update, file='~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/mathys_new_celltypes_seurat4.csv', row.names=FALSE)
file <- synapser::File(path='~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/mathys_new_celltypes_seurat4.csv', parentId='syn25871777')
file <- synapser::synStore(file)


saveRDS(mathys3, file='~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/mathys_mapped_seuratObj.rds')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/mathys_mapped_seuratObj.rds', parentId='syn25871777')
file <- synapser::synStore(file)

