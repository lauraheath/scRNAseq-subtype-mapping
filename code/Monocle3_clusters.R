#DO NOT RUN IN SAME CONTAINER USING MONOCLE2 ANALYSIS. Monocle2 and monocle3 should not be installed in same container!

#use monocle3 to break the counts matrix into subclass labels, as determined by the allen mapping.
#Stratify by sex, recluster, and reassign labels into individual clusters for trajectory analysis 
#(reclassifying necessary for sublcass types that encompass more than one discrete cluster or 
#clusters that contain multiple subclass types according to the mapping)


#get scran-normalized counts matrix:
#counts2 <- readRDS(file="~/scRNAseq-subtype-mapping/data/mathys_scrannorm_counts.rds")



#raw counts matrix, not normalized, with labeled rownames and colnames:
p1 <- synapser::synGet('syn25871862')
counts <- readRDS(p1$path)


#get allen celltype predictions
mathys_metadata <- read.csv(file="~/scRNAseq-subtype-mapping/data/mathys_new_celltypes_seurat4.csv")
#p <- synapser::synGet('syn25871851')
#mathys_metadata <- read.csv(p$path)
rownames(mathys_metadata) <- mathys_metadata$TAG
mathys_metadata$X<-NULL
#get original mathys metadata (with properly ordered rows to match counts matrix)
#meta_temp <- readRDS(file="~/scRNAseq-subtype-mapping/data/mathys_metadata.rds")
#p2 <- synapser::synGet('syn24171852')
#meta_temp <- readRDS(p2$path)

#get short gene names list and make them into rownames on counts file: filtered_gene_row_names.txt
p3 <- synapser::synGet('syn18686382')
gene_short_name <- readLines(p3$path)
gene_short_name <- (data.frame(gene_short_name))
rownames(gene_short_name) <- gene_short_name$gene_short_name


#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
mathys_metadata <- mathys_metadata[order(match(mathys_metadata$TAG, colnames_counts$columnnames)),]



cds <- new_cell_data_set(counts,
                         cell_metadata = mathys_metadata,
                         gene_metadata = gene_short_name)


cds$educ = as.numeric(cds$educ)
cds$batch = as.factor(cds$batch)
cds$Diagnosis[cds$Diagnosis=='Control']='Cont'
cds$Diagnosis[cds$Diagnosis=='Late-pathology AD']='Late'
cds$Diagnosis[cds$Diagnosis=='Early-pathology AD']='Early'
cds$simpleDiagnosis = cds$Diagnosis
cds$simpleDiagnosis[cds$simpleDiagnosis!='Cont'] <- "AD"
cds$msex <- ifelse(cds$sex=='female',0,1)

cds = preprocess_cds(cds, num_dim = 30,method="PCA", norm_method = "log")
plot_pc_variance_explained(cds)
cds <- align_cds(cds, 
                 preprocess_method="PCA",
                 #alignment_group="batch",
                 alignment_group="ros_ids",
                 residual_model_formula_str="~educ+pmi+sex")
cds = reduce_dimension(cds)
cds = cluster_cells(cds, cluster_method="louvain")
plot_cells(cds, color_cells_by="predicted.subclass_label")
plot_cells(cds, color_cells_by="partition")




#subset cds by celltype, parse out matrix/metadata/gene short name
dim(cds)





############### refine microglia cluster for analysis ########
#isolate microglia: get all cells in partitions 9 and 10, along with all predicted microglia
cds$partition <- partitions(cds)
table(cds$partition)
table(cds$predicted.subclass_label)
mic_cells <- row.names(subset(pData(cds),
                              #partition==10 |
                              #partition==9 |
                              predicted.subclass_label=='Micro-PVM'))

cds_mic <- cds[,mic_cells]
dim(cds_mic)
plot_cells(cds_mic, color_cells_by="predicted.subclass_label")
plot_cells(cds_mic, color_cells_by="predicted.subclass_label.score", cell_size=2)
#recluster
cds_mic <- clear_cds_slots(cds_mic)
cds_mic$ros_ids <- as.character(cds_mic$ros_ids)
cds_mic$ros_ids <- as.factor(cds_mic$ros_ids)
table(cds_mic$ros_ids)
cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method = "log")
plot_pc_variance_explained(cds_mic)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="ros_ids",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain")
plot_cells(cds_mic, color_cells_by="predicted.subclass_label",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


plot_cells(cds_mic, color_cells_by="predicted.subclass_label.score",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


#keep partition 1
cds_mic$partition <- partitions(cds_mic)
table(cds_mic$partition)
cds_mic <- cds_mic[,cds_mic$partition==1]
plot_cells(cds_mic, color_cells_by="predicted.subclass_label",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="predicted.subclass_label.score",cell_size=2,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="partition",cell_size=2,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


#relabel all the cells in the cluster as microglial cells:
cds_mic$refined_subclass <- 'Micro-PVM'
#save the cds for monocle2 analysis:
saveRDS(cds_mic, file="~/scRNAseq-subtype-mapping/Glial_clusters/cds_mic.RDS")
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/cds_mic.RDS', parentId='syn25871777')
file <- synapser::synStore(file)


#save elements of the cds into synapse to upload into separate notebook for monocle2 analysis
microglia_counts <- exprs(cds_mic)
microglia_genes <- data.frame(fData(cds_mic))
microglia_meta <- data.frame(pData(cds_mic))

saveRDS(microglia_counts, file='~/scRNAseq-subtype-mapping/Glial_clusters/microglia_counts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/microglia_counts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

saveRDS(microglia_genes, file='~/scRNAseq-subtype-mapping/Glial_clusters/microglia_genes.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/microglia_genes.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

saveRDS(microglia_meta, file='~/scRNAseq-subtype-mapping/Glial_clusters/microglia_meta.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/microglia_meta.RDS', parentId='syn25871777')
file <- synapser::synStore(file)











############### refine astrocyte cluster for analysis ########
plot_cells(cds, color_cells_by="predicted.subclass_label")
plot_cells(cds, color_cells_by="partition")
#isolate astrocytes: get all cells in partition 12, along with all predicted astrocytes
table(cds$predicted.subclass_label)
ast_cells <- row.names(subset(pData(cds),
                              # partition==12 |
                              predicted.subclass_label=='Astro'))


cds_ast <- cds[,ast_cells]
dim(cds_ast)
plot_cells(cds_ast, color_cells_by="predicted.subclass_label")
#recluster
cds_ast <- clear_cds_slots(cds_ast)
cds_ast$ros_ids <- as.character(cds_ast$ros_ids)
cds_ast$ros_ids <- as.factor(cds_ast$ros_ids)
table(cds_ast$ros_ids)
cds_ast = preprocess_cds(cds_ast, num_dim = 30,method="PCA", norm_method = "log")
plot_pc_variance_explained(cds_ast)
cds_ast <- align_cds(cds_ast, 
                     preprocess_method="PCA",
                     alignment_group="ros_ids",
                     residual_model_formula_str="~educ+pmi")
cds_ast = reduce_dimension(cds_ast)
cds_ast = cluster_cells(cds_ast, cluster_method="louvain")
plot_cells(cds_ast, color_cells_by="predicted.subclass_label",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_ast, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


plot_cells(cds_ast, color_cells_by="predicted.subclass_label.score",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_ast, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#all in one partition
#relabel all the cells in the cluster as Astro cells:
cds_ast$refined_subclass <- 'Astro'
#save the cds for monocle2 analysis:
saveRDS(cds_ast, file="~/scRNAseq-subtype-mapping/Glial_clusters/cds_ast.RDS")
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/cds_ast.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

#save elements of the cds into synapse to upload into separate notebook for monocle2 analysis
astrocyte_counts <- exprs(cds_ast)
astrocyte_genes <- data.frame(fData(cds_ast))
astrocyte_meta <- data.frame(pData(cds_ast))

saveRDS(astrocyte_counts, file='~/scRNAseq-subtype-mapping/Glial_clusters/astrocyte_counts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/astrocyte_counts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

saveRDS(astrocyte_genes, file='~/scRNAseq-subtype-mapping/Glial_clusters/astrocyte_genes.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/astrocyte_genes.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

saveRDS(astrocyte_meta, file='~/scRNAseq-subtype-mapping/Glial_clusters/astrocyte_meta.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/astrocyte_meta.RDS', parentId='syn25871777')
file <- synapser::synStore(file)



############### refine oligodendrocyte cluster for analysis ########
plot_cells(cds, color_cells_by="predicted.subclass_label")
plot_cells(cds, color_cells_by="partition")
#isolate oligos: all predicted Oligos
table(cds$partition)
table(cds$predicted.subclass_label)
oli_cells <- row.names(subset(pData(cds),
                              #partition==6 |
                              predicted.subclass_label=='Oligo'))


cds_oli <- cds[,oli_cells]
dim(cds_oli)
plot_cells(cds_oli, color_cells_by="predicted.subclass_label")
#recluster
cds_oli <- clear_cds_slots(cds_oli)
cds_oli$ros_ids <- as.character(cds_oli$ros_ids)
cds_oli$ros_ids <- as.factor(cds_oli$ros_ids)
table(cds_oli$ros_ids)
cds_oli = preprocess_cds(cds_oli, num_dim = 30,method="PCA", norm_method = "log")
plot_pc_variance_explained(cds_oli)
cds_oli <- align_cds(cds_oli, 
                     preprocess_method="PCA",
                     alignment_group="ros_ids",
                     residual_model_formula_str="~educ+pmi")
cds_oli = reduce_dimension(cds_oli)
cds_oli = cluster_cells(cds_oli, cluster_method="louvain")
plot_cells(cds_oli, color_cells_by="predicted.subclass_label",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_oli, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


plot_cells(cds_oli, color_cells_by="predicted.subclass_label.score",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_oli, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#all in one partition
#relabel all the cells in the cluster as oligodendrocyte cells:
cds_oli$refined_subclass <- 'Oligo'
#save the cds for monocle2 analysis:
saveRDS(cds_oli, file="~/scRNAseq-subtype-mapping/Glial_clusters/cds_oli.RDS")
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/cds_oli.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

#save elements of the cds into synapse to upload into separate notebook for monocle2 analysis
oligo_counts <- exprs(cds_oli)
oligo_genes <- data.frame(fData(cds_oli))
oligo_meta <- data.frame(pData(cds_oli))

saveRDS(astrocyte_counts, file='~/scRNAseq-subtype-mapping/Glial_clusters/oligo_counts.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/oligo_counts.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

saveRDS(astrocyte_genes, file='~/scRNAseq-subtype-mapping/Glial_clusters/oligo_genes.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/oligo_genes.RDS', parentId='syn25871777')
file <- synapser::synStore(file)

saveRDS(astrocyte_meta, file='~/scRNAseq-subtype-mapping/Glial_clusters/oligo_meta.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/Glial_clusters/oligo_meta.RDS', parentId='syn25871777')
file <- synapser::synStore(file)
