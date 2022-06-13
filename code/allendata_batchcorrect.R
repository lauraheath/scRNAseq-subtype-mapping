library(Seurat)
library(SingleCellExperiment)
library(scran)
library(batchelor)
library(biomaRt)

setwd("~/scRNAseq-subtype-mapping/")

#NOT sample-corrected data for Differential expression between ADNC groups
p1 <- synapser::synGet('syn31618644')
AllenAD <- readRDS(p1$path)
head(x=AllenAD[[]])
table(AllenAD$sex)
length(unique(AllenAD$`Donor ID`))
AllenAD

counts <- GetAssayData(object = AllenAD, slot = "counts")
metadata <- AllenAD@meta.data
metadata$TAG <- rownames(metadata)




###########################################################################################
# 
#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
metadata <- metadata[order(match(metadata$TAG, colnames_counts$columnnames)),]

#replace ensembl identifiers with gene short names
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='https://www.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}

Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}

ensembl <- row.names(counts)
ensembl <- as.data.frame(ensembl)
ensembl$gene_short_name <- Make.Gene.Symb(ensembl$ensembl)
gene_short_name <- ensembl$gene_short_name
ensembl$gene_short_name2 <- ifelse(ensembl$gene_short_name=="", ensembl$ensembl, ensembl$gene_short_name)

gene_short_name2 <- ensembl$gene_short_name2
rownames(counts) <- ensembl$gene_short_name2

#delete genes that did not get translated by biomaRt (virtually all are RNA genes, pseudogenes, and have low expression)
counts2 <- counts[-which(grepl("^ENSG",counts@Dimnames[[1]])),]
counts3 <- counts2[-which(grepl("^[0-9]",counts2@Dimnames[[1]])),]

dim(counts3)

#remake seurat object, which will fix column headers with spaces in the metadata and make rownames all unique
allen <- CreateSeuratObject(counts = counts3, meta.data = metadata, min.cells = 3)
dim(allen)
#calculate size factors for DE analysis later:
libsizes <- colSums(counts3)
allen <- AddMetaData(allen, metadata=libsizes, col.name='libsizes')
allen$size.factor <- allen$libsizes/mean(allen$libsizes)


#remove ADNC Reference group and non-AD individuals
allen <- subset(allen, ADNC!='Reference')
allen <- subset(allen, ADNC!='Not AD')
allen$ADNC<-as.character(allen$ADNC)
allen$ADNC<-as.factor(allen$ADNC)
length(unique(allen$Donor.ID))
table(allen$Donor.ID)
dim(allen)

saveRDS(allen, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroMF_Seurat.RDS") 
#save to synapse:
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroMF_Seurat.RDS', parentId='syn31924906')
file <- synapser::synStore(file)


#batch correct by sex to shrink sample size:
AstroFemale <- subset(allen, sex=='female')
counts <- GetAssayData(object = AstroFemale, slot = "counts")
metadata <- AstroFemale@meta.data

AstroF <- CreateSeuratObject(counts = counts, project = "AllenSCRNAseq", min.cells=3)
AstroF <- AddMetaData(AstroF, metadata)
head(x=AstroF[[]])
dim(AstroF)

#AllenAD2[["percent.mt"]] <- PercentageFeatureSet(AllenAD2, pattern = "^MT-")
#VlnPlot(AllenAD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#mitochondrial genes already removed
FeatureScatter(AstroF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

length(unique(AstroF$Donor.ID))
#AstroM <- subset(AllenAD, sex=='male')


#To visualize the cells, run Seurat pipeline to normalize data and perform dimension reduction
AstroF <- SCTransform(AstroF, do.scale = TRUE, verbose = TRUE, conserve.memory=TRUE)

# AstroF <- NormalizeData(AstroF, normalization.method = "LogNormalize", scale.factor = 10000)
# AstroF <- FindVariableFeatures(AstroF, selection.method = "vst", nfeatures = 2000)
# 
# top10 <- head(VariableFeatures(AstroF), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(AstroF)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# all.genes <- rownames(AstroF)
# AstroF <- ScaleData(AstroF, features = all.genes)




AstroF <- RunPCA(AstroF, npcs = 30)
AstroF <- RunUMAP(AstroF, dims = 1:30)
#for cleaner clusters, decrease resolution in FindNeighbors
AstroF <- FindNeighbors(AstroF, reduction = "umap", dims = 1:2)
AstroF <- FindClusters(AstroF, resolution = 0.5, verbose = TRUE)
DimPlot(AstroF, group.by = "Supertype", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(AstroF, reduction="umap", label = TRUE, repel=TRUE)
DimPlot(AstroF, group.by = "ADNC", reduction="umap", label = TRUE, repel=TRUE)




#For batch correction: extract counts and metadata from seurat object, and make a singlecellexperiment object
counts <- GetAssayData(object = AstroF, slot = "counts")
#counts <- GetAssayData(object = glialcellsM, slot = "counts")
counts <- as.matrix(counts)
metadata <- AstroF@meta.data

#create single single experiment object
sce <- SingleCellExperiment(list(counts=counts), 
                            colData=metadata)
sce
head(assay(sce[0:20,0:20]))
head(colData(sce))

clusters <- quickCluster(sce, min.size=100)
sce <- computeSumFactors(sce, cluster=clusters)
#sce <- computeSumFactors(sce, cluster=clusters, BPPARAM=MulticoreParam(8))
summary(sizeFactors(sce))

sce <- logNormCounts(sce)
head(logcounts(sce[0:20,0:20]))
head(rownames(logcounts(sce)))

dim(logcounts(sce))




sce$batch <- as.factor(sce$Donor.ID)

#sce4 <- fastMNN(sce, batch = sce$ros_ids, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
sce4 <- fastMNN(sce, batch = sce$Donor.ID, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)

#sce4 <- fastMNN(sce, batch = sce$ros_ids, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce5 <- fastMNN(sce, batch = sce$batch, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce4 <- fastMNN(sce, batch = sce$batch, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)


corrected_sample <- assay(sce4)
#corrected_batch <- assay(sce5)


saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected.RDS")
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

#save female metadata
metadata <- AstroF@meta.data
metadata$TAG <- rownames(metadata)
write.csv(metadata, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata.csv", row.names = FALSE)
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata.csv', parentId='syn31924906')
file <- synapser::synStore(file)


