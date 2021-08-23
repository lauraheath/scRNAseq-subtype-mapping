library(scran)
library(scater)
library(Matrix)
library(monocle3)



#get allen celltype predictions
p <- synapser::synGet('syn25871851')
mathys_metadata <- read.csv(p$path)
rownames(mathys_metadata) <- mathys_metadata$TAG
#get mathys counts matrix:
p1 <- synapser::synGet('syn18686381')
counts <- readMM(p1$path)
#get original mathys metadata (with properly ordered rows to match counts matrix)
#meta_temp <- readRDS(file="~/celltype-mapping/data/mathys_metadata.rds")
p2 <- synapser::synGet('syn24171852')
meta_temp <- readRDS(p2$path)
#get short gene names list and make them into rownames on counts file: filtered_gene_row_names.txt
p3 <- synapser::synGet('syn18686382')

#assign short gene names to rownames of counts matrix
rownames(counts) <- readLines(p3$path)
#assign rownames of original metadata file to column names of counts matrix (cell names)
colnames(counts) <- rownames(meta_temp)

#create dataframe of gene list for monocle3:
gene_short_name <- readLines(p3$path)
gene_short_name <- (data.frame(gene_short_name))
rownames(gene_short_name) <- gene_short_name$gene_short_name

#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
mathys_metadata <- mathys_metadata[order(match(mathys_metadata$TAG, colnames_counts$columnnames)),]
rownames(mathys_metadat)

meta_temp<-0

#normalize counts matrix in scran:
sce <- SingleCellExperiment(list(counts=counts))
dim(sce)
clusters <- quickCluster(sce, min.size=100)
sce <- computeSumFactors(sce, cluster=clusters)
#sce <- computeSumFactors(sce, cluster=clusters, BPPARAM=MulticoreParam(8))
summary(sizeFactors(sce))

sce <- logNormCounts(sce)
head(logcounts(sce[0:20,0:20]))
head(rownames(logcounts(sce)))

dim(logcounts(sce))
counts2 = logcounts(sce)
dim(counts2)

#EXPORT SCRAN NORMALIZED COUNTS MATRIX TO SYNAPSE
saveRDS(counts2, file="~/celltype_mapping/data/mathys_scrannorm_counts.rds")
file <- synapser::File(path='~/celltype_mapping/data/mathys_scrannorm_counts.rds', parentId='syn25871777')
file <- synapser::synStore(file)

#scran normalized count matrix:
#counts <- readRDS(file="mathys_scran_normalized_counts2.rds")
p <- synapser::synGet('syn22363128')
counts_saved <- readRDS(p$path)
dim(counts_saved)
head(counts_saved)
#extract count matrix from the cds
counts <- exprs(cds)


cds <- new_cell_data_set(counts2,
                         cell_metadata = mathys_metadata,
                         gene_metadata = gene_short_name)


#cds <- readRDS(file="~/scAD_analysis/Mathys_scrannorm_cds.rds")


cds$educ = as.numeric(cds$educ)
cds$Education = cds$educ
cds$Sex = cds$sex
cds$CERAD = cds$ceradsc
cds$Braak = cds$braaksc
cds$APOE_genotype = cds$apoe_genotype
cds$batch = as.factor(cds$batch)
cds$Diagnosis[cds$Diagnosis=='Control']='Cont'
cds$Diagnosis[cds$Diagnosis=='Late-pathology AD']='Late'
cds$Diagnosis[cds$Diagnosis=='Early-pathology AD']='Early'
cds$simpleDiagnosis = cds$Diagnosis
cds$simpleDiagnosis[cds$simpleDiagnosis!='Cont'] <- "AD"


##Separate cds by sex for sex-specific analyses from here on:
#make sex-specific cds and clear the cds slots and rescale data
cdsF = cds[,cds$sex=='female']
dim(cdsF)
cdsF <- clear_cds_slots(cdsF)
cdsF$batch <- as.character(cdsF$batch)
cdsF$batch <- as.factor(cdsF$batch)
#Preprocessing the data using PCA. Data is already size-normalized and log transformed.
cdsF = preprocess_cds(cdsF, num_dim = 30,method="PCA", norm_method = "none")
plot_pc_variance_explained(cdsF)
cdsF <- align_cds(cdsF, 
                  preprocess_method="PCA",
                  #alignment_group="batch",
                  alignment_group="ros_ids",
                  residual_model_formula_str="~educ+pmi")
cdsF = reduce_dimension(cdsF)
cdsF = cluster_cells(cdsF, cluster_method="louvain")
plot_cells(cdsF, color_cells_by="predicted.subclass_label")
plot_cells(cdsF, color_cells_by="partition")
plot_cells(cdsF, color_cells_by="broad.cell.type")


#Male cds 
cdsM = cds[,cds$sex=='male']
dim(cdsM)
cdsM <- clear_cds_slots(cdsM)
cdsM$batch <- as.character(cdsM$batch)
cdsM$batch <- as.factor(cdsM$batch)
cdsM = preprocess_cds(cdsM, num_dim = 30,method="PCA", norm_method = "none")
plot_pc_variance_explained(cdsM)
cdsM <- align_cds(cdsM,
                  preprocess_method="PCA",
                  alignment_group="ros_ids",
                  residual_model_formula_str="~educ+pmi")
cdsM = reduce_dimension(cdsM)
cdsM = cluster_cells(cdsM, cluster_method="louvain")
plot_cells(cdsM, color_cells_by="predicted.subclass_label")

### start with female samples only, glial cells

cds_mic = cdsF[,cdsF$predicted.subclass_label=='Micro-PVM']
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

cds_mic$partition <- partitions(cds_mic)

cds_badmatch <- cds_mic[,cds_mic$partition!=10]
cds_badmatch <- cds_badmatch[,cds_badmatch$partition!=9]
dim(cds_badmatch)
head(colData(cds_badmatch))
hist(cds_badmatch$predicted.subclass_label.score )


cds_mic <- clear_cds_slots(cds_mic)
cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method="none")
#adjust for sample effects
cds_mic$ros_ids <- as.character(cds_mic$ros_ids)
cds_mic$ros_ids <- as.factor(cds_mic$ros_ids)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="ros_ids",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain", k=500)
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="broad.",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
