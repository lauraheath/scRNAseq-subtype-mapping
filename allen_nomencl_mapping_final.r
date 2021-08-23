
#upload the gene expression matrix from allen institute's M1 reference data:

p1 <- synapser::synGet('syn24182839')
M1ref <- read.csv(p1$path)
#turn the csv file into a sparse matrix and transpose, so cells are columns and genes are rows
rownames(M1ref) <- M1ref$sample_name
M1ref$sample_name<-NULL
M1ref <- as.matrix(M1ref)
M1ref <- Matrix(M1ref, sparse=TRUE)
M1ref <- t(M1ref)
saveRDS(M1ref, file="~/scRNAseq-subtype-mapping/data/M1ref_counts.rds")
M1ref <- readRDS(file="~/scRNAseq-subtype-mapping/data/M1ref_counts.rds")

#read in metadata for cells:
p2 <- synapser::synGet('syn24182843')
metadata <- read.csv(p2$path)
rownames(metadata) <- colnames(M1ref)

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
#get all mathys metadata
#mathys_meta <- readRDS(file="~/scRNAseq-subtype-mapping/data/mathys_metadata.rds")
p2 <- synapser::synGet('syn24171852')
mathys_meta <- readRDS(p2$path)
#get short gene names list and make them into rownames on counts file: filtered_gene_row_names.txt
p3 <- synapser::synGet('syn18686382')
rownames(counts) <- readLines(p3$path)
colnames(counts) <- mathys_meta[,1]

#save the data structures:
saveRDS(counts, file="~/scRNAseq-subtype-mapping/data/mathys_scrannorm_counts.rds")
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
p1 + p2 + plot_layout(guides = "collect")
mathys3 <- merge(mathys2.batches[[1]], mathys2.batches[2:length(mathys2.batches)], merge.dr = "ref.umap")
DimPlot(mathys3, reduction = "ref.umap", group.by =  "predicted.subclass_label", label = FALSE, repel = TRUE, label.size = 1)
head(x=mathys3[[]])

hist(mathys3$predicted.subclass_label.score)
hist(mathys3$predicted.subclass_label.score)



#######################################
mathys4 <- merge(mathys2.batches[[1]], mathys2.batches[2:length(mathys2.batches)], merge.dr = "ref.pca")
dat$id <- 'reference'
mathys4$id <- 'query'
refquery <- merge(dat, mathys4)
refquery[["pca"]] <- merge(dat[["pca"]], mathys4[["ref.pca"]])
refquery <- RunUMAP(refquery, reduction = 'pca', dim = 1:20)
DimPlot(refquery, group.by = 'id')
DimPlot(refquery, group.by = 'predicted.subclass_label')
DimPlot(refquery, group.by = 'broad.cell.type')




#### repeat with control subjects ONLY ###

control <- subset(x = mathys2, subset = Diagnosis == "Cont")
head(x=control[[]])
control <- SCTransform(control, vars.to.regress = "ros_ids", do.scale = TRUE, verbose = TRUE)
head(x=control[[]])


anchors <- FindTransferAnchors(
  reference = dat,
  query = control,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:20
)

control <- MapQuery(
  anchorset = anchors,
  query =control,
  reference = dat,
  refdata = list(
    class_label = "class_label",
    subclass_label = "subclass_label",
    cluster_label = "cluster_label"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)
head(x=control[[]])

DimPlot(control, reduction = "ref.umap", group.by = "predicted.class_label", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(control, reduction = "ref.umap", group.by = "predicted.subclass_label", label = FALSE, label.size = 1, repel = TRUE)
hist(control$predicted.subclass_label.score)
summary(control$predicted.subclass_label.score)

#now with AD subjects only
ad <- subset(x = mathys2, subset = Diagnosis != "Cont")
head(x=ad[[]])
table(ad$Diagnosis)
ad <- SCTransform(ad, vars.to.regress = "ros_ids", do.scale = TRUE, verbose = TRUE)
head(x=ad[[]])


anchors <- FindTransferAnchors(
  reference = dat,
  query = ad,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:20
)

ad <- MapQuery(
  anchorset = anchors,
  query =ad,
  reference = dat,
  refdata = list(
    class_label = "class_label",
    subclass_label = "subclass_label",
    cluster_label = "cluster_label"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)
head(x=ad[[]])

DimPlot(ad, reduction = "ref.umap", group.by = "predicted.class_label", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(ad, reduction = "ref.umap", group.by = "predicted.subclass_label", label = FALSE, label.size = 1, repel = TRUE)
hist(ad$predicted.subclass_label.score)
summary(ad$predicted.subclass_label.score)
hist(control$predicted.subclass_label.score)



#### read in scran-normalized mathys data and attach predicted celltypes to metadata
cds <- readRDS(file="~/scAD_analysis/Mathys_scrannorm_cds.rds")
#get the metadata:
mathys_meta2 <- mathys2@meta.data

cds$pred.class_label <- mathys_meta2$predicted.class_label
cds$pred.class_label.score <- mathys_meta2$predicted.class_label.score
cds$pred.subclass_label <- mathys_meta2$predicted.subclass_label
cds$pred.subclass_label.score <- mathys_meta2$predicted.subclass_label.score
cds$pred.cluster_label <- mathys_meta2$predicted.cluster_label
cds$pred.cluster_label.score <- mathys_meta2$predicted.cluster_label.score
head(colData(cds))
#want to calculate subclass percentages in control and AD
df1 <- mathys_meta2 %>% group_by(Diagnosis, predicted.subclass_label) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

ggplot(df1, aes(fill = Diagnosis, y = percent, x = predicted.subclass_label))+
  geom_bar(position = "dodge", stat = "identity")+ theme_classic()

df2 <- mathys_meta2 %>% group_by(simpleDiagnosis, predicted.subclass_label) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

ggplot(df2, aes(fill = simpleDiagnosis, y = percent, x = predicted.subclass_label))+
  geom_bar(position = "dodge", stat = "identity")+ theme_classic()



cont_meta <- control@meta.data
ad_meta <- ad@meta.data
df3 <- rbind(cont_meta, ad_meta)
df3 <- df3 %>% group_by(Diagnosis, predicted.subclass_label) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

ggplot(df3, aes(fill = Diagnosis, y = percent, x = predicted.subclass_label))+
  geom_bar(position = "dodge", stat = "identity")+ theme_classic()


ggplot(mathys_meta2, aes(y = predicted.subclass_label.score, x = predicted.subclass_label, fill = Diagnosis))+
  geom_boxplot()+ theme_classic()+geom_point(position=position_jitter())



co <- subset(mathys_meta2, mathys_meta2$Diagnosis=='Cont')


