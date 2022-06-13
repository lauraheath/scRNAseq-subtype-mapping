BiocManager::install("scran")

library(SingleCellExperiment)
library(scran)
library(batchelor)


#raw counts matrix, not normalized, with labeled rownames and colnames:
p1 <- synapser::synGet('syn25871862')
counts <- readRDS(p1$path)
counts <- as(counts, "dgCMatrix")

#upload metadata with allen nomenclature mapping attached (mathys_meta_reclassifiedCellTypes.csv):
p <- synapser::synGet('syn29262930')
mathys_metadata <- read.csv(p$path)
rownames(mathys_metadata) <- mathys_metadata$TAG


mathys <- CreateSeuratObject(counts = counts)

mathys2 <- AddMetaData(mathys, mathys_metadata)
head(x=mathys2[[]])

#save this seurat object with attached metadata and reclassified cell types into the synapse space
saveRDS(mathys2, file="~/scRNAseq-subtype-mapping/data_objects/SeuratObj_mathys_reclassified.RDS")
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/SeuratObj_mathys_reclassified.RDS', parentId='syn25871777')
file <- synapser::synStore(file)




#upload relabeled seurat objects by sex from seurat clustering:
#mathys <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/allcellsF_Seurat.RDS")
# mathys <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/allcellsM_Seurat.RDS")
# dim(mathys)

glialcells <- subset(mathys2, reclassify=='Astro'|reclassify=='Micro-PVM'|reclassify=='OPC'|reclassify=='Oligo')
table(glialcells$reclassify)

glialcellsF <- subset(glialcells, sex=='female')
glialcellsM <- subset(glialcells, sex=='male')

#pull metadata from seurat object
metadata <- glialcellsF@meta.data
#metadata <- glialcellsM@meta.data

counts <- GetAssayData(object = glialcellsF, slot = "counts")
#counts <- GetAssayData(object = glialcellsM, slot = "counts")
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

mergeorder <- list(list('ROS1','ROS3','ROS5','ROS2','ROS10','ROS11','ROS12','ROS4','ROS8','ROS6',
                        'ROS7','ROS9'), list('ROS35','ROS28','ROS36','ROS26','ROS32','ROS25',
                                             'ROS27','ROS34','ROS33','ROS30','ROS29','ROS31'))
#chosen.hvgs <- getTopHVGs(sce, n=3000)

sce$batch <- as.factor(sce$batch)
sce$batch <- as.factor(sce$ros_ids)

#sce4 <- fastMNN(sce, batch = sce$ros_ids, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
sce4 <- fastMNN(sce, batch = sce$ros_ids, correct.all=TRUE, cos.norm=TRUE, merge.order = mergeorder, prop.k=0.05)

#sce4 <- fastMNN(sce, batch = sce$ros_ids, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce5 <- fastMNN(sce, batch = sce$batch, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce4 <- fastMNN(sce, batch = sce$batch, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)


corrected_sample <- assay(sce4)
#corrected_batch <- assay(sce5)


table(metadata$batch, metadata$ros_ids)

saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_samplecorrected.RDS")
saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_samplecorrected_orderedmerge.RDS")

#saveRDS(corrected_batch, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_batchcorrected.RDS")
saveRDS(metadata, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_metadata.RDS")

saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_samplecorrected.RDS")
#saveRDS(corrected_batch, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_batchcorrected.RDS")
saveRDS(metadata, file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_metadata.RDS")





