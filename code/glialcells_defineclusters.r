#install.packages("remotes")
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

#upload metadata with allen nomenclature mapping attached:
p <- synapser::synGet('syn25871851')
mathys_metadata <- read.csv(p$path)
rownames(mathys_metadata) <- mathys_metadata$TAG
mathys_metadata$X<-NULL


#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
mathys_metadata <- mathys_metadata[order(match(mathys_metadata$TAG, colnames_counts$columnnames)),]


mathys3 <- CreateSeuratObject(counts = counts, project = "Clustering")

seurObj <- CreateSeuratObject(counts = counts)
head(x=seurObj[[]])
#get the UMI counts:
UMIcounts <- seurObj@meta.data

mathys2 <- AddMetaData(mathys2, mathys_metadata)
head(x=mathys3[[]])

#split into male and female data sets
mathysF <- subset(x=mathys2, subset=sex=='female')
dim(mathysF)

mathysM <- subset(x=mathys2, subset=sex=='male')
dim(mathysM)


mathysF <- SCTransform(mathysF, do.scale = TRUE, verbose = TRUE)
mathysM <- SCTransform(mathysM, do.scale = TRUE, verbose = TRUE)

mathysF <- RunPCA(mathysF, npcs = 30)
mathysF <- RunUMAP(mathysF, dims = 1:30)
mathysM <- RunPCA(mathysM, npcs = 30)
mathysM <- RunUMAP(mathysM, dims = 1:30)


#for cleaner clusters, decrease resolution in FindNeighbors
mathysF <- FindNeighbors(mathysF, reduction = "umap", dims = 1:2)
mathysF <- FindClusters(mathysF, resolution = 0.1, verbose = TRUE)
DimPlot(mathysF, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysF, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysF, reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysF, group.by = "predicted.class_label", reduction="umap", label = TRUE, repel=TRUE)

mathysM <- FindNeighbors(mathysM, reduction = "umap", dims = 1:2)
mathysM <- FindClusters(mathysM, resolution = 0.05, verbose = TRUE)
DimPlot(mathysM, group.by = "predicted.subclass_label", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysM, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysM, reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysM, group.by = "predicted.class_label", reduction="umap", label = TRUE, repel=TRUE)


#metadata <- glialcells@meta.data
mathysF$clusters <- mathysF$seurat_clusters
mathysF$reclassify <- ifelse(mathysF$clusters==12, 'Micro-PVM',
                             ifelse(mathysF$clusters==6, 'Astro',
                                    ifelse(mathysF$clusters==0|mathysF$clusters==2, 'Oligo',
                                           ifelse(mathysF$clusters==11, 'OPC',
                                                  ifelse(mathysF$clusters==10|mathysF$clusters==16|mathysF$clusters==15|mathysF$clusters==9, 'In',
                                                      ifelse(mathysF$clusters==19|mathysF$clusters==20|mathysF$clusters==21, 'other', 'Ex'))))))
table(mathysF$reclassify)
table(mathysF$predicted.subclass_label)
table(mathysF$predicted.class_label)



mathysM$clusters <- mathysM$seurat_clusters
mathysM$reclassify <- ifelse(mathysM$clusters==11, 'Micro-PVM',
                             ifelse(mathysM$clusters==6, 'Astro', 
                                    ifelse(mathysM$clusters==0, 'Oligo',
                                           ifelse(mathysM$clusters==7, 'OPC',
                                                  ifelse(mathysM$clusters==12|mathysM$clusters==9|mathysM$clusters==4|mathysM$clusters==5, 'In',
                                                         ifelse(mathysM$clusters==15|mathysM$clusters==16|mathysM$clusters==17, 'other', 'Ex'))))))
table(mathysM$reclassify)
table(mathysM$predicted.subclass_label)
table(mathysM$predicted.class_label)

DimPlot(mathysF, group.by = "reclassify", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(mathysM, group.by = "reclassify", reduction="umap", label = TRUE, repel=TRUE)

#save seurat objects with new classifications, and matching metadata
metadataF <- mathysF@meta.data
saveRDS(mathysF, file="~/scRNAseq-subtype-mapping/data_objects/allcellsF_Seurat.RDS")
saveRDS(metadataF, file="~/scRNAseq-subtype-mapping/data_objects/allcellsF_metadata.RDS")

metadataM <- mathysM@meta.data
saveRDS(mathysM, file="~/scRNAseq-subtype-mapping/data_objects/allcellsM_Seurat.RDS")
saveRDS(metadataM, file="~/scRNAseq-subtype-mapping/data_objects/allcellsM_metadata.RDS")









############# based on the broad cell types in reclassify, find DE genes (non-conservative simple
#model to allow for more potential DE genes rather than a restrictive list)

### oligodendrocyte DE genes (control vs all pathology, control vs early, early vs late)

#reset levels in overall data set:
Idents(object=mathysF) <- mathysF@meta.data$Diagnosis
levels(mathysF)
table(mathysF$Diagnosis)


oligo <- subset(mathysF, reclassify=='Oligo')
dim(oligo)
table(oligo$Diagnosis)
CvsAll <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = NULL, logfc.threshold = 0.1)
head(CvsAll)
oligo_genes1 <- subset(CvsAll, p_val_adj<0.05)
oligo_genes1$DEgenes <- rownames(oligo_genes1)

CvsE <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Early", logfc.threshold = 0.1)
head(CvsE)
oligo_genes2 <- subset(CvsE, p_val_adj<0.05)
oligo_genes2$DEgenes <- rownames(oligo_genes2)

EvsL <- FindMarkers(oligo, ident.1 = "Early", ident.2 = "Late", logfc.threshold = 0.1)
head(EvsL)
oligo_genes3 <- subset(EvsL, p_val_adj<0.05)
oligo_genes3$DEgenes <- rownames(oligo_genes3)

CvsL <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Late", logfc.threshold = 0.1)
head(CvsL)
oligo_genes4 <- subset(CvsL, p_val_adj<0.05)
oligo_genes4$DEgenes <- rownames(oligo_genes4)

oligo_genes <- rbind(oligo_genes1, oligo_genes2, oligo_genes3, oligo_genes4)
oligo_DE <- unique(oligo_genes$DEgenes)

### Microglia DE genes (control vs all pathology, control vs early, early vs late)
oligo <- subset(mathysF, reclassify=='Micro-PVM')
dim(oligo)

CvsAll <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = NULL, logfc.threshold = 0.1)
head(CvsAll)
oligo_genes1 <- subset(CvsAll, p_val_adj<0.05)
oligo_genes1$DEgenes <- rownames(oligo_genes1)

CvsE <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Early", logfc.threshold = 0.1)
head(CvsE)
oligo_genes2 <- subset(CvsE, p_val_adj<0.05)
oligo_genes2$DEgenes <- rownames(oligo_genes2)

EvsL <- FindMarkers(oligo, ident.1 = "Early", ident.2 = "Late", logfc.threshold = 0.1)
head(EvsL)
oligo_genes3 <- subset(EvsL, p_val_adj<0.05)
oligo_genes3$DEgenes <- rownames(oligo_genes3)

CvsL <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Late", logfc.threshold = 0.1)
head(CvsL)
oligo_genes4 <- subset(CvsL, p_val_adj<0.05)
oligo_genes4$DEgenes <- rownames(oligo_genes4)

micro_genes <- rbind(oligo_genes1, oligo_genes2, oligo_genes3, oligo_genes4)
micro_DE <- unique(micro_genes$DEgenes)


### Astrocyte DE genes (control vs all pathology, control vs early, early vs late)



oligo <- subset(mathysF, reclassify=='Astro')
dim(oligo)
CvsAll <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = NULL, logfc.threshold = 0.1)
head(CvsAll)
oligo_genes1 <- subset(CvsAll, p_val_adj<0.05)
oligo_genes1$DEgenes <- rownames(oligo_genes1)

CvsE <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Early", logfc.threshold = 0.1)
head(CvsE)
oligo_genes2 <- subset(CvsE, p_val_adj<0.05)
oligo_genes2$DEgenes <- rownames(oligo_genes2)

EvsL <- FindMarkers(oligo, ident.1 = "Early", ident.2 = "Late", logfc.threshold = 0.1)
head(EvsL)
oligo_genes3 <- subset(EvsL, p_val_adj<0.05)
oligo_genes3$DEgenes <- rownames(oligo_genes3)

CvsL <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Late", logfc.threshold = 0.1)
head(CvsL)
oligo_genes4 <- subset(CvsL, p_val_adj<0.05)
oligo_genes4$DEgenes <- rownames(oligo_genes4)

astro_genes <- rbind(oligo_genes1, oligo_genes2, oligo_genes3, oligo_genes4)
astro_DE <- unique(astro_genes$DEgenes)

### OPC DE genes (control vs all pathology, control vs early, early vs late)
oligo <- subset(mathysF, reclassify=='OPC')
dim(oligo)
CvsAll <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = NULL, logfc.threshold = 0.1)
head(CvsAll)
oligo_genes1 <- subset(CvsAll, p_val_adj<0.05)
oligo_genes1$DEgenes <- rownames(oligo_genes1)

CvsE <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Early", logfc.threshold = 0.1)
head(CvsE)
oligo_genes2 <- subset(CvsE, p_val_adj<0.05)
oligo_genes2$DEgenes <- rownames(oligo_genes2)

EvsL <- FindMarkers(oligo, ident.1 = "Early", ident.2 = "Late", logfc.threshold = 0.1)
head(EvsL)
oligo_genes3 <- subset(EvsL, p_val_adj<0.05)
oligo_genes3$DEgenes <- rownames(oligo_genes3)

CvsL <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Late", logfc.threshold = 0.1)
head(CvsL)
oligo_genes4 <- subset(CvsL, p_val_adj<0.05)
oligo_genes4$DEgenes <- rownames(oligo_genes4)

opc_genes <- rbind(oligo_genes1, oligo_genes2, oligo_genes3, oligo_genes4)
opc_DE <- unique(opc_genes$DEgenes)


### Excitatory neurons DE genes (control vs all pathology, control vs early, early vs late)
oligo <- subset(mathysF, reclassify=='Ex')
dim(oligo)
Idents(object = oligo) <- oligo@meta.data$Diagnosis
levels(oligo)
CvsAll <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = NULL, logfc.threshold = 0.1)
head(CvsAll)
oligo_genes1 <- subset(CvsAll, p_val_adj<0.05)
oligo_genes1$DEgenes <- rownames(oligo_genes1)

CvsE <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Early", logfc.threshold = 0.1)
head(CvsE)
oligo_genes2 <- subset(CvsE, p_val_adj<0.05)
oligo_genes2$DEgenes <- rownames(oligo_genes2)

EvsL <- FindMarkers(oligo, ident.1 = "Early", ident.2 = "Late", logfc.threshold = 0.1)
head(EvsL)
oligo_genes3 <- subset(EvsL, p_val_adj<0.05)
oligo_genes3$DEgenes <- rownames(oligo_genes3)

CvsL <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Late", logfc.threshold = 0.1)
head(CvsL)
oligo_genes4 <- subset(CvsL, p_val_adj<0.05)
oligo_genes4$DEgenes <- rownames(oligo_genes4)

ex_genes <- rbind(oligo_genes1, oligo_genes2, oligo_genes3, oligo_genes4)
ex_DE <- unique(ex_genes$DEgenes)



### Inhibitory neurons DE genes (control vs all pathology, control vs early, early vs late)
oligo <- subset(mathysF, reclassify=='In')
dim(oligo)
Idents(object = oligo) <- oligo@meta.data$Diagnosis
levels(oligo)
CvsAll <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = NULL, logfc.threshold = 0.1)
head(CvsAll)
oligo_genes1 <- subset(CvsAll, p_val_adj<0.05)
oligo_genes1$DEgenes <- rownames(oligo_genes1)

CvsE <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Early", logfc.threshold = 0.1)
head(CvsE)
oligo_genes2 <- subset(CvsE, p_val_adj<0.05)
oligo_genes2$DEgenes <- rownames(oligo_genes2)

EvsL <- FindMarkers(oligo, ident.1 = "Early", ident.2 = "Late", logfc.threshold = 0.1)
head(EvsL)
oligo_genes3 <- subset(EvsL, p_val_adj<0.05)
oligo_genes3$DEgenes <- rownames(oligo_genes3)

CvsL <- FindMarkers(oligo, ident.1 = "Cont", ident.2 = "Late", logfc.threshold = 0.1)
head(CvsL)
oligo_genes4 <- subset(CvsL, p_val_adj<0.05)
oligo_genes4$DEgenes <- rownames(oligo_genes4)

in_genes <- rbind(oligo_genes1, oligo_genes2, oligo_genes3, oligo_genes4)
in_DE <- unique(in_genes$DEgenes)


wilcox_DEgenes <- as.vector(c(in_DE, ex_DE, micro_DE, opc_DE, astro_DE, oligo_DE))
wilcox_DEgenes <- unique(wilcox_DEgenes)
wilcox_DEgenes <- as.data.frame(wilcox_DEgenes)





#subset the oligodendrocytes (#15)
DimPlot(glialcells, reduction="umap", label = TRUE, repel=TRUE)

oligo <- subset(glialcells, clusters==1|clusters==2|clusters==4)
oligo_meta <- oligo@meta.data
counts <- GetAssayData(object=oligo, slot="counts")


oligo <- CreateSeuratObject(counts = counts, project = "Clustering")

oligo <- AddMetaData(oligo, oligo_meta)
head(x=oligo[[]])


oligo <- NormalizeData(oligo)
oligo <- FindVariableFeatures(oligo, nfeatures = 3000)
oligo <- RunFastMNN(object.list = SplitObject(oligo, split.by = "batch"), features=3000)
oligo <- ScaleData(oligo, verbose = FALSE)
oligo <- RunUMAP(oligo, reduction = "mnn", dims = 1:30)
oligo <- FindNeighbors(oligo, reduction = "mnn", dims = 1:30)
oligo <- FindClusters(oligo, resolution=0.5)

DimPlot(oligo, reduction="umap", label = TRUE, repel=TRUE)
DimPlot(oligo, group.by = "batch", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(oligo, group.by = "ros_ids", reduction="umap", label = FALSE, repel=TRUE)
DimPlot(oligo, group.by = "Subcluster", reduction="umap", label = FALSE, repel=TRUE)

oligomarkers <- FindAllMarkers(oligo, only.pos=TRUE, min.pct=0.2, logfc.threshold=0.2)
oligomarkers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


oligo$oligo_clusters <- oligo$seurat_clusters
table(oligo$oligo_clusters)
#need to add these back to the metadata data frame
oligo$reclassify <- ifelse(oligo$oligo_clusters==0, 'Oligo0',
                           ifelse(oligo$oligo_clusters==1, 'Oligo1',
                                  ifelse(oligo$oligo_clusters==2, 'Oligo2', 'Oligo3')))


glial_meta <- glialcells@meta.data
glial_meta <- subset(glial_meta, reclassify!='Oligo')
oligo_meta <- oligo@meta.data
oligo_meta$RNA_snn_res.0.5<-NULL
oligo_meta$oligo_clusters<-NULL
glial_meta2 <- rbind(glial_meta, oligo_meta)
table(glial_meta2$reclassify)


glialcells <- GetAssayData(object=glialcells, slot="counts")
glialcells <- as.matrix(glialcells)
saveRDS(glialcells, file="~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/glialcells.RDS")
saveRDS(glial_meta2, file="~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/glialcells_metadata.RDS")







