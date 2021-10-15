BiocManager::install(c("scran", "scater"))
library(scran)
library(SingleCellExperiment)
library(batchelor)


#upload glial cells counts and metadata from seurat clustering:
#mathys <- readRDS(file="~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/mathysF_Seurat.RDS")
mathys <- readRDS(file="~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/mathysM_Seurat.RDS")

#metadata <- readRDS(file="~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/glialcellsF_metadata.RDS")

glialcells <- subset(mathys, reclassify!='other')
#pull metadata from seurat object
metadata <- glialcells@meta.data

counts <- GetAssayData(object = glialcells, slot = "counts")
counts <- as.matrix(counts)


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

chosen.hvgs <- getTopHVGs(sce, n=3000)

#sce$batch <- as.factor(sce$batch)
#sce3 <- mnnCorrect(sce, batch = sce$batch, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm.in=TRUE, cos.norm.out=FALSE)

sce4 <- fastMNN(sce, batch = sce$batch, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce4 <- fastMNN(sce, batch = sce$batch, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)


corrected <- assay(sce4)

table(metadata$batch, metadata$ros_ids)

#saveRDS(corrected, file="~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/glialcellsF_batchcorrected.RDS")
saveRDS(corrected, file="~/scRNAseq-subtype-mapping/scRNAseq-subtype-mapping/glialcellsM_batchcorrected.RDS")


