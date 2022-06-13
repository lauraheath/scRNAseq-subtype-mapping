#Rerun differential expression on each broad cell type, using nebula
#library(Seurat)
devtools::install_github("lhe17/nebula")
library(nebula)
install.packages("VennDiagram")
library(VennDiagram)
remotes::install_github("satijalab/sctransform", ref="develop")
install.packages('Seurat')
#remotes::install_github("satijalab/seurat-wrappers")
library(Seurat)


#raw counts matrix, not normalized, with labeled rownames and colnames:
p1 <- synapser::synGet('syn25871862')
counts <- readRDS(p1$path)
counts <- as(counts, "dgCMatrix")

p <- synapser::synGet('syn29262930')
mathys_metadata <- read.csv(p$path)
rownames(mathys_metadata) <- mathys_metadata$TAG




#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
mathys_metadata <- mathys_metadata[order(match(mathys_metadata$TAG, colnames_counts$columnnames)),]


#create datasets by broad celltype (microglia, astrocyte, oligodendrocyte, OPC, excitatory neurons, inhibitory neurons)
#create a seurat object because it is easier to work with
mathys <- CreateSeuratObject(counts = counts, meta.data = mathys_metadata)

micro <- subset(mathys, reclassify=='Micro-PVM')
astro <- subset(mathys, reclassify=='Astro')
oligo <- subset(mathys, reclassify=='Oligo')
opc <- subset(mathys, reclassify=='OPC')
ex <- subset(mathys, reclassify=='Ex')
inter <- subset(mathys, reclassify=='In')



# counts2 <- GetAssayData(object=micro, slot="counts")
# counts2 <- as.matrix(counts2)
# dim(counts2)
# 
# counts2 <- GetAssayData(object=astro, slot="counts")
# counts2 <- as.matrix(counts2)
# dim(counts2)

#counts2 <- GetAssayData(object=ex2, slot="counts")
# counts2 <- GetAssayData(object=ex, slot="counts")
# counts2 <- as.matrix(counts2)
# dim(counts2)


#counts2 <- GetAssayData(object=micro, slot="counts")
counts2 <- GetAssayData(object=astro, slot="counts")
#counts2 <- GetAssayData(object=oligo, slot="counts")
#counts2 <- GetAssayData(object=opc, slot="counts")
counts2 <- as.matrix(counts2)
dim(counts2)
libsizes <- colSums(counts2)
astro <- AddMetaData(astro, metadata=libsizes, col.name='libsizes')
astro$size.factor <- astro$libsizes/mean(astro$libsizes)
#micro <- AddMetaData(micro, metadata=libsizes, col.name='libsizes')
#micro$size.factor <- micro$libsizes/mean(astro$libsizes)



#create categories for the nebula comparisons:
#test1: control vs. any pathology
#test2: control vs. early-pathology
#test3: early-pathology vs. late-pathology
test1 <- astro
test1$diag <- ifelse(test1$Diagnosis=='Control', 0, 1)
test2 <- subset(astro, Diagnosis!='Late-pathology AD')
test2$diag <- ifelse(test2$Diagnosis=='Control', 0, 1)
test3 <- subset(astro, Diagnosis!='Control')
test3$diag <- ifelse(test3$Diagnosis=='Early-pathology AD', 0, 1)


# test1 <- micro
# test1$diag <- ifelse(test1$Diagnosis=='Control', 0, 1)
# test2 <- subset(micro, Diagnosis!='Late-pathology AD')
# test2$diag <- ifelse(test2$Diagnosis=='Control', 0, 1)
# test3 <- subset(micro, Diagnosis!='Control')
# test3$diag <- ifelse(test3$Diagnosis=='Early-pathology AD', 0, 1)


###### test1 ##########
#model test1:

#order the metadata by subject ID, then reorder the counts matrix to match:
meta <- test1@meta.data
mat <- GetAssayData(object=test1, slot="counts")
mat <- as.matrix(mat)
meta <- meta[order(meta$ros_ids),]
col.order <- meta$TAG
mat <- mat[,col.order]
#set design matrix and run nebula
modmat = model.matrix(~diag+sex+pmi+batch, data=meta)
model1 <- nebula(mat, meta$ros_ids, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_diag, "fdr")
model1sum <- cbind(adjustp, model1sum)
model1sum$model <- 'nopath_vs_path'
head(model1sum)


#model test2:
#order the metadata by subject ID, then reorder the counts matrix to match:
meta <- test2@meta.data
mat <- GetAssayData(object=test2, slot="counts")
mat <- as.matrix(mat)
meta <- meta[order(meta$ros_ids),]
col.order <- meta$TAG
mat <- mat[,col.order]
#set design matrix and run nebulaa
modmat = model.matrix(~diag+sex+pmi+batch, data=meta)
model2 <- nebula(mat, meta$ros_ids, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model2sum <- model2$summary
adjustp <- p.adjust(model2sum$p_diag, "fdr")
model2sum <- cbind(adjustp, model2sum)
model2sum$model <- 'nopath_vs_earlypath'
head(model2sum)

#model test3
meta <- test3@meta.data
mat <- GetAssayData(object=test3, slot="counts")
mat <- as.matrix(mat)
meta <- meta[order(meta$ros_ids),]
col.order <- meta$TAG
mat <- mat[,col.order]
modmat = model.matrix(~diag+sex+pmi+batch, data=meta)
model3 <- nebula(mat, meta$ros_ids, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model3sum <- model3$summary
adjustp <- p.adjust(model3sum$p_diag, "fdr")
model3sum <- cbind(adjustp, model3sum)
model3sum$model <- 'earlypath_vs_latepath'
head(model3sum)



DEgenes <- rbind(model1sum, model2sum)
DEgenes <- rbind(DEgenes, model3sum)


write.csv(DEgenes, file='~/scRNAseq-subtype-mapping/dataframes/astro_DE_genesALL.csv', row.names=FALSE)
#write.csv(DEgenes, file='~/scRNAseq-subtype-mapping/dataframes/micro_DE_genesALL.csv', row.names=FALSE)


pvalue05 <- DEgenes[(DEgenes$p_diag<0.05),]
length(unique(pvalue05$gene))





#compare DEgenes from mathys to DE genes from nebula in astrocytes:
#upload the list of DE genes from my synapse stash:
DE1obj <- synapser::synGet('syn29791645')
mathysDEastro <- read.csv(DE1obj$path)

#what is the overlap between significant genes in mathys and significant genes from nebula?
#compare genes significant under unadjusted MixedModel.p to nebula p
mathysDEold <- mathysDEastro[(mathysDEastro$MixedModel.p<0.05),]
mathysDEold <- mathysDEold[!is.na(mathysDEold$MixedModel.p),]
table(mathysDEold$model)

oldgenes <- mathysDEold$gene
oldgenes <- unique(oldgenes)
newgenes <- pvalue05$gene
newgenes <- unique(newgenes)
#library(VennDiagram)
venn.diagram(
  x=list(oldgenes,newgenes),
  category.names = c("Original", "Update"),
  filename = 'DEgene_overlap.png',
  output=TRUE
)
















#write.csv(model1sum, file='~/scRNAseq-subtype-mapping/micro1.csv', row.names=FALSE)
#write.csv(model1sum, file='~/scRNAseq-subtype-mapping/astro1.csv', row.names=FALSE)
#write.csv(model1sum, file='~/scRNAseq-subtype-mapping/oligo1.csv', row.names=FALSE)
#write.csv(model1sum, file='~/scRNAseq-subtype-mapping/opc1.csv', row.names=FALSE)




###### test1 ##########
#model exadj1:
subjectID <- meta$ros_ids
#modmat = model.matrix(~test1, data=meta)
modmat = model.matrix(~test1+sex+pmi+batch, data=meta)
model1 <- nebula(counts2, meta$ros_ids, pred=modmat, offset=meta$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)

write.csv(model1sum, file='~/scRNAseq-subtype-mapping/micro1.csv', row.names=FALSE)

#model exadj2:
modmat = model.matrix(~test1, data=meta)
model1 <- nebula(counts2, meta$ros_ids, pred=modmat, offset=meta$libsizes)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)


write.csv(model1sum, file='~/scRNAseq-subtype-mapping/micro2.csv', row.names=FALSE)


#model exadj3:
modmat = model.matrix(~test1 + seurat_clusters, data=meta)
model1 <- nebula(counts2, meta$ros_ids, pred=modmat, offset=meta$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)



#model exadj4: test the zero seurat_cluster only
ex0 <- subset(ex, seurat_clusters==0)
counts0 <- GetAssayData(object=ex0, slot="counts")
counts0 <- as.matrix(counts0)
dim(counts0)

libsizes <- colSums(counts0)
size.factors <- as.data.frame(libsizes/mean(libsizes))
meta0 <- ex0@meta.data
meta0 <- cbind(meta0, size.factors)
names(meta0)[names(meta0) == "libsizes/mean(libsizes)"] <- "size.factor"

#order the metadata by subject ID, then reorder the counts matrix to match:
meta0 <- meta0[order(meta0$ros_ids),]
col.order <- meta0$TAG
counts0 <- counts0[,col.order]
modmat = model.matrix(~test1, data=meta0)
model1 <- nebula(counts0, meta0$ros_ids, pred=modmat, offset=meta0$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)


#model exadj5: test the mathys-assigned Ex0 cluster (the biggest cluster)
ex0 <- subset(mathys, Subcluster=='Ex0')
counts0 <- GetAssayData(object=ex0, slot="counts")
counts0 <- as.matrix(counts0)
dim(counts0)

libsizes <- colSums(counts0)
size.factors <- as.data.frame(libsizes/mean(libsizes))
meta0 <- ex0@meta.data
meta0 <- cbind(meta0, size.factors)
names(meta0)[names(meta0) == "libsizes/mean(libsizes)"] <- "size.factor"

#order the metadata by subject ID, then reorder the counts matrix to match:
meta0 <- meta0[order(meta0$ros_ids),]
col.order <- meta0$TAG
counts0 <- counts0[,col.order]
modmat = model.matrix(~test1, data=meta0)
model1 <- nebula(counts0, meta0$ros_ids, pred=modmat, offset=meta0$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)




#model exadj6: test the mathys-assigned Ex broad celltype cluster with size.factor offset
ex0 <- subset(mathys, broad.cell.type=='Ex')
counts0 <- GetAssayData(object=ex0, slot="counts")
counts0 <- as.matrix(counts0)
dim(counts0)
libsizes <- colSums(counts0)
size.factors <- as.data.frame(libsizes/mean(libsizes))
meta0 <- ex0@meta.data
meta0 <- cbind(meta0, size.factors)
names(meta0)[names(meta0) == "libsizes/mean(libsizes)"] <- "size.factor"
#order the metadata by subject ID, then reorder the counts matrix to match:
meta0 <- meta0[order(meta0$ros_ids),]
col.order <- meta0$TAG
counts0 <- counts0[,col.order]
modmat = model.matrix(~test1, data=meta0)
model1 <- nebula(counts0, meta0$ros_ids, pred=modmat, offset=meta0$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)

#model exadj7: test the mathys-assigned Ex broad celltype cluster withOUT size.factor offset
modmat = model.matrix(~test1, data=meta0)
model1 <- nebula(counts0, meta0$ros_ids, pred=modmat)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)



#try no path vs. early path ex (exadj8):
ex <- subset(mathys, broad.cell.type=='Ex')
ex <- subset(ex, test2==0|test2==1)
counts0 <- GetAssayData(object=ex, slot="counts")
counts0 <- as.matrix(counts0)
dim(counts0)
libsizes <- colSums(counts0)
size.factors <- as.data.frame(libsizes/mean(libsizes))
meta0 <- ex@meta.data
meta0 <- cbind(meta0, size.factors)
names(meta0)[names(meta0) == "libsizes/mean(libsizes)"] <- "size.factor"
#order the metadata by subject ID, then reorder the counts matrix to match:
meta0 <- meta0[order(meta0$ros_ids),]
col.order <- meta0$TAG
counts0 <- counts0[,col.order]
modmat = model.matrix(~test2, data=meta0)
model1 <- nebula(counts0, meta0$ros_ids, pred=modmat, offset=meta0$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)




#model in1:
modmat = model.matrix(~test1, data=meta)
model1 <- nebula(counts2, meta$ros_ids, pred=modmat, offset=meta$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)

#model in2:
modmat = model.matrix(~test1, data=meta)
model1 <- nebula(counts2, meta$ros_ids, pred=modmat, offset=meta$size.factor)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_test1, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)

write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj1.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj2.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj3.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj4.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj5.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj6.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj7.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exadj8.csv', row.names=FALSE)

write.csv(model1sum, file='~/scRNAseq-subtype-mapping/in1.csv', row.names=FALSE)

write.csv(model1sum, file='~/scRNAseq-subtype-mapping/microglia.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/microglia2PMM.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/microglia_adj3.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/microglia_adj4.csv', row.names=FALSE)

write.csv(model1sum, file='~/scRNAseq-subtype-mapping/astro_adj4.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/astro_adj5.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exNeuro_adj4.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exNeuro_adj5.csv', row.names=FALSE)
write.csv(model1sum, file='~/scRNAseq-subtype-mapping/exmathys_adj5.csv', row.names=FALSE)
