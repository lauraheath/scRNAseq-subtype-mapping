
library(SingleCellExperiment)
library(scran)
library(batchelor)

#upload relabeled seurat objects by sex from seurat clustering:
mathys <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/allcellsF_Seurat.RDS")
mathys <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/allcellsM_Seurat.RDS")
dim(mathys)

glialcells <- subset(mathys, reclassify=='Astro'|reclassify=='Micro-PVM'|reclassify=='OPC'|reclassify=='Oligo')
table(glialcells$reclassify)
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

#chosen.hvgs <- getTopHVGs(sce, n=3000)

sce$batch <- as.factor(sce$batch)
#sce3 <- mnnCorrect(sce, batch = sce$batch, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm.in=TRUE, cos.norm.out=FALSE)

sce4 <- fastMNN(sce, batch = sce$ros_ids, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
sce5 <- fastMNN(sce, batch = sce$batch, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce4 <- fastMNN(sce, batch = sce$batch, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)


corrected_sample <- assay(sce4)
corrected_batch <- assay(sce5)


table(metadata$batch, metadata$ros_ids)

saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_samplecorrected.RDS")
saveRDS(corrected_batch, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_batchcorrected.RDS")
saveRDS(metadata, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_metadata.RDS")

saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_samplecorrected.RDS")
saveRDS(corrected_batch, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_batchcorrected.RDS")
saveRDS(metadata, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_metadata.RDS")





