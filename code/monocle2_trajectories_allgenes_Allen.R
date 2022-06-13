install.packages('lme4')

library(monocle)
library(dplyr)
library(ggplot2)
library(lme4)
library(nebula)
library(biomaRt)
synapser::synLogin()


##############################################################################
#Differential expression by ADNC status

allen <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroMF_Seurat.RDS") 
#create categories for the nebula comparisons:
#test1: ADNC Low vs Intermediate/High combined
#test2: ADNC Low vs Intermediate
#test3: ADNC Intermediate vs High

#remove donors who do not have AD (some other pathology)
test1 <- allen
test1$diag <- ifelse(test1$ADNC=='Low', 0, 1)
test2 <- subset(allen2, ADNC!='High')
test2$diag <- ifelse(test2$ADNC=='Low', 0, 1)
test3 <- subset(allen2, ADNC!='Low')
test3$diag <- ifelse(test3$ADNC=='Intermediate', 0, 1)

table(test1$diag, test1$ADNC)
table(test2$diag, test2$ADNC)
table(test3$diag, test3$ADNC)

#order the metadata by subject ID, then reorder the counts matrix to match:
meta <- test1@meta.data
mat <- GetAssayData(object=test1, slot="counts")
#mat <- as.matrix(mat)
meta <- meta[order(meta$Donor.ID),]
col.order <- meta$TAG
mat <- mat[,col.order]
#set design matrix and run nebula
modmat = model.matrix(~diag+sex, data=meta)
model1 <- nebula(mat, meta$Donor.ID, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the comparison, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_diag, "fdr")
model1sum <- cbind(adjustp, model1sum)
model1sum$model <- 'Low_vs_All'
head(model1sum)


#model test2:
#order the metadata by subject ID, then reorder the counts matrix to match:
meta <- test2@meta.data
mat <- GetAssayData(object=test2, slot="counts")
meta <- meta[order(meta$Donor.ID),]
col.order <- meta$TAG
mat <- mat[,col.order]
#set design matrix and run nebulaa
modmat = model.matrix(~diag+sex, data=meta)
model2 <- nebula(mat, meta$Donor.ID, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model2sum <- model2$summary
adjustp <- p.adjust(model2sum$p_diag, "fdr")
model2sum <- cbind(adjustp, model2sum)
model2sum$model <- 'Low_vs_Intermediate'
head(model2sum)

#model test3
meta <- test3@meta.data
mat <- GetAssayData(object=test3, slot="counts")
meta <- meta[order(meta$Donor.ID),]
col.order <- meta$TAG
mat <- mat[,col.order]
modmat = model.matrix(~diag+sex, data=meta)
model3 <- nebula(mat, meta$Donor.ID, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model3sum <- model3$summary
adjustp <- p.adjust(model3sum$p_diag, "fdr")
model3sum <- cbind(adjustp, model3sum)
model3sum$model <- 'Int_vs_High'
head(model3sum)



DEgenes <- rbind(model1sum, model2sum)
DEgenes <- rbind(DEgenes, model3sum)


write.csv(DEgenes, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_DE_genesALL.csv', row.names=FALSE)
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_DE_genesALL.csv', parentId='syn31924906')
file <- synapser::synStore(file)

pvalue05 <- DEgenes[(DEgenes$p_diag<0.05),]
length(unique(pvalue05$gene))



#upload sample-corrected counts matrix
counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected.RDS")
counts <- as.matrix(counts)
metadata <- read.csv(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata.csv")
dim(counts)
dim(metadata)


DEgenes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_DE_genesALL.csv')
DEgenes05 <- subset(DEgenes, DEgenes$p_diag<0.05)
length(unique(DEgenes05$gene))


genes2 <- c()
for (gene in unique(c(as.vector(DEgenes05$gene)))){
  #for (gene in unique(c(as.vector(wilcox_DEgenes$genes)))){
  
  if (gene %in% rownames(counts)){
    genes2 <- c(genes2,which(rownames(counts)==gene))
  }
}

counts2 <- counts[genes2,]
dim(counts2)




gene_short_name <- rownames(counts)
dim(counts)
gene_short_name <- as.data.frame(gene_short_name)
rownames(gene_short_name)<-gene_short_name$gene_short_name

temp <- counts
temp2 <- metadata

rownames(temp) <- NULL
colnames(temp) <- NULL



pd <- new("AnnotatedDataFrame", data = temp2)
fd <- new("AnnotatedDataFrame", data = gene_short_name)



MonRun <- newCellDataSet(as.matrix(temp),
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = tobit())

MonRun <- reduceDimension(MonRun, max_components=2, reduction_method = 'DDRTree', norm_method = 'none', pseudo_expr = 0, verbose = TRUE)


#save the Monocle object with reduced dimensions:
saveRDS(MonRun, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_reduceddim.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_reduceddim.RDS', parentId='syn31924906')
file <- synapser::synStore(file)


MonRun <- orderCells(MonRun)





RunMonocleTobit2 <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                             gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily = gaussianff())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, norm_method = 'none', pseudo_expr = 0)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}



####### For non-batch-corrected matrix, run RunMonocleTobit. For batch-corrected matrix, run RunMonocleTobit2 ###

MonRun <- RunMonocleTobit2(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroF_samplecorrected_tree.tiff',height=125,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AstroM_samplecorrected_tree.tiff',height=125,width=200,units='mm',res=300)

#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroF_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroM_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis") + guides(colour=guide_legend(override.aes=list(size=3)))
g
dev.off()


g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State") + guides(colour=guide_legend(override.aes=list(size=3)))
g

g<- plot_cell_trajectory(MonRun,color_by = "Subcluster",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="subcluster")
g
plot_cell_trajectory(MonRun,color_by = "Subcluster",show_branch_points=F,use_color_gradient = F,cell_size = 1)


plot_cell_trajectory(MonRun,color_by = "ros_ids",show_branch_points=F,use_color_gradient = F,cell_size = 1)
plot_cell_trajectory(MonRun,color_by = "ros_ids",show_branch_points=F,use_color_gradient = F,cell_size = 1)
MonRun$batch <- as.factor(MonRun$batch)
plot_cell_trajectory(MonRun,color_by = "batch",show_branch_points=F,use_color_gradient = F,cell_size = 1)

#determine which state has highest proportion of control cells compared to AD cells
df <- list()
df$State <- MonRun$State
df$Diagnosis <- MonRun$Diagnosis
df <- as.data.frame(df)
table(df$State)
proptable <- with(df, table(Diagnosis, State)) %>% prop.table(margin = 2)
proptable

table(MonRun$State)

#If necessary, reset root state to the state with the highest proportion of control cells

#Astro F: no correction to root 
#microglia females sample-corrected: root at state 6
#MonRun <- orderCells(MonRun, root_state = 6)
#microglia males sample-corrected: root at state 3
#MonRun <- orderCells(MonRun, root_state = 3)


plot_cell_trajectory(MonRun,color_by = "Pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1)

plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)


#relabel states:
#AstroF 
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 9] <- 2
MonRun$State2[MonRun$State == 3] <- 3
MonRun$State2[MonRun$State == 4] <- 4
MonRun$State2[MonRun$State == 5] <- 5
MonRun$State2[MonRun$State == 6] <- 6
MonRun$State2[MonRun$State ==8] <- 7
MonRun$State2[MonRun$State == 7] <- 8


#AstroM
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 1
MonRun$State2[MonRun$State == 1] <- 2
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 9] <- 3
MonRun$State2[MonRun$State == 8] <- 4
MonRun$State2[MonRun$State == 3] <- 5
MonRun$State2[MonRun$State == 4] <- 5
MonRun$State2[MonRun$State == 5] <- 6

 

MonRun$State2 <- as.character(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)
plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1)
df <- list()
df$State2 <- MonRun$State2
df$Diagnosis <- MonRun$Diagnosis
df <- as.data.frame(df)
table(df$State2)
proptable <- with(df, table(Diagnosis, State2)) %>% prop.table(margin = 2)
proptable
plot_cell_trajectory(MonRun,color_by = "Subcluster",show_branch_points=F,use_color_gradient = F,cell_size = 1)

#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroF_state_tree.tiff',height=125,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AstroM_state_tree.tiff',height=125,width=200,units='mm',res=300)

#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroF_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroM_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_state_tree.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State") + guides(colour=guide_legend(override.aes=list(size=3)))
g
dev.off()


#tiff(file='~/scRNAseq-subtype-mapping/figures/astroF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/astroM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)

#tiff(file='~/scRNAseq-subtype-mapping/figures/microF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/microM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
MonRun$braaksc <- as.factor(MonRun$braaksc)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()

#tiff(file='~/scRNAseq-subtype-mapping/figures/astroF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/astroM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)

#tiff(file='~/scRNAseq-subtype-mapping/figures/microF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/microM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
MonRun$cogdx <- as.factor(MonRun$cogdx)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="COGDX",y="Pseudotime",x="cogdx")
g
dev.off()


#tiff(file='~/scRNAseq-subtype-mapping/figures/astroF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/astroM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)

#tiff(file='~/scRNAseq-subtype-mapping/figures/microF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/microM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
MonRun$ceradsc <- as.factor(MonRun$ceradsc)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="ceradsc")
g
dev.off()

head(pData(MonRun))
x <- list()
x$cellID <- MonRun$TAG
x$projid <- MonRun$projid
x$ros_ids <- MonRun$ros_ids
x$State <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$Diagnosis <- MonRun$Diagnosis
x$diag2 <- MonRun$simpleDiagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$amyloid <- MonRun$amyloid
x$plaq_n <- MonRun$plaq_n
x$nft <- MonRun$nft
x$tangles <- MonRun$tangles
x$gpath <- MonRun$gpath
x$amyloid.group <- MonRun$amyloid.group
x$apoe <- MonRun$apoe_genotype
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi
x$reclassify <- MonRun$reclassify

#females: rename and create a scaled pseudotime variable
Fvariables <- as.data.frame(x)
Fvariables$pseudotime_sc <- scale(Fvariables$Pseudotime, center=F)

#save variables file for later
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/MicroF_samplecorrected_pstimeStates.csv", row.names=FALSE)
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/MicroM_samplecorrected_pstimeStates.csv", row.names=FALSE)

write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AstroF_allgenes_pstimeStates.csv", row.names=FALSE)
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AstroM_samplecorrected_pstimeStates.csv", row.names=FALSE)

#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OligoF_samplecorrected_pstimeStates.csv", row.names=FALSE)
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OligoM_samplecorrected_pstimeStates.csv", row.names=FALSE)

#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OpcF_samplecorrected_pstimeStates.csv", row.names=FALSE)
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OpcM_samplecorrected_pstimeStates.csv", row.names=FALSE)


#run logistic regression comparing pseudotiem between cases and controls only
#casecontrolF <- subset(Fvariables, Fvariables$simpleDiagnosis=='AD'|Fvariables$simpleDiagnosis=='Cont')

Fvariables$diag3 <- ifelse(Fvariables$diag2=='AD', 1, 0)
Fvariables$ADvsNoPath <- ifelse(Fvariables$diag3==0, 'No pathology', 'AD pathology')
Fvariables$ADvsNoPath <- factor(Fvariables$ADvsNoPath, levels=c('No pathology', 'AD pathology'))

summary(glm(diag3 ~ pseudotime_sc,Fvariables,family='binomial'))

#control vs early
CvsE <- subset(Fvariables, Fvariables$Diagnosis=='Control' | Fvariables$Diagnosis=='Early-pathology AD')
CvsE$diag <- ifelse(CvsE$Diagnosis=='Control', 0, 1)
summary(glm(diag ~ pseudotime_sc,CvsE,family='binomial'))

#control vs late
CvsL <- subset(Fvariables, Fvariables$Diagnosis=='Control' | Fvariables$Diagnosis=='Late-pathology AD')
CvsL$diag <- ifelse(CvsL$Diagnosis=='Control', 0, 1)
summary(glm(diag ~ pseudotime_sc,CvsL,family='binomial'))

#early vs late
EvsL <- subset(Fvariables, Fvariables$Diagnosis=='Early-pathology AD' | Fvariables$Diagnosis=='Late-pathology AD')
EvsL$diag <- ifelse(EvsL$Diagnosis=='Early-pathology AD', 0, 1)
summary(glm(diag ~ pseudotime_sc,EvsL,family='binomial'))

#summary(lmer(diag3 ~ pseudotime_sc + pmi + (pseudotime_sc|ros_ids), Fvariables))
#summary(lmer(diag3 ~ pseudotime_sc + (1|ros_ids), Fvariables))


#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroF_bargraph_3cat_diagnosis.tiff',height=170,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AstroM_bargraph_3cat_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(Fvariables,aes(x=Diagnosis,
                             y=pseudotime_sc,
                             color=Diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:3]) 
g
dev.off()

#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroF_bargraph_2cat_diagnosis.tiff',height=170,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AstroM_bargraph_2cat_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(Fvariables,aes(x=ADvsNoPath,
                           y=pseudotime_sc,
                           color=ADvsNoPath)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:3])
g
dev.off()


Fvariables$Diagnosis2 <- ifelse(Fvariables$Diagnosis=='Cont', 0,
                                ifelse(Fvariables$Diagnosis=='Early', 1, 2))
Fvariables$Diagnosis2 <- factor(Fvariables$Diagnosis2,levels = c(0:2))
Fvariables$braaksc <- factor(Fvariables$braaksc,levels = c(1:6))
Fvariables$ceradsc <- factor(Fvariables$ceradsc,levels = c(1:4))
Fvariables$cogdx <- factor(Fvariables$cogdx,levels = c(1:6))


#run proportional odds logistic regression for neuropath/cognitive endpoints:
diagfit <- MASS::polr(Diagnosis2 ~ pseudotime_sc,Fvariables)
braakfit <- MASS::polr(braaksc ~ pseudotime_sc,Fvariables)
ceradfit <- MASS::polr(ceradsc ~ pseudotime_sc,Fvariables)
cogdxfit <- MASS::polr(cogdx ~ pseudotime_sc,Fvariables)

cat('diagnosis p-value: ',pt(abs(summary(diagfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')



#BEAM (branch point genes)
plot_cell_trajectory(MonRun, color_by = "State")
BEAM_res <- BEAM(MonRun, branch_point = , cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

#by patient
table(MonRun$ros_ids)
MonRun2 <- MonRun
MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS13', 'ROS13', 'OTHER')
plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)

MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS14', 'ROS14', 'OTHER')
plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)

MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS15', 'ROS15', 'OTHER')
plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)

MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS16', 'ROS16', 'OTHER')
plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)

MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS17', 'ROS17', 'OTHER')
plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
