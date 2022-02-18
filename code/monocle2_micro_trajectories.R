

library(monocle)
library(dplyr)
library(ggplot2)
synapser::synLogin()

##############################################################################
#use sample-corrected glial cell matrix from fastMNN as matrix:
#counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_batchcorrected.RDS")
counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_samplecorrected.RDS")
counts <- as.matrix(counts)
metadata <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/glialcellsF_metadata.RDS")
dim(counts)
dim(metadata)

#counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_batchcorrected.RDS")
counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_samplecorrected.RDS")
counts <- as.matrix(counts)
metadata <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/glialcellsM_metadata.RDS")
dim(counts)
dim(metadata)


#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(counts))
names(colnames_counts)[names(colnames_counts) == "colnames(counts)"] <- "columnnames"
metadata2 <- metadata[order(match(metadata$TAG, colnames_counts$columnnames)),]


###########################################################################################

#extract celltype of interest
#cells <- which(metadata$reclassify == 'Micro-PVM' & metadata$ros_ids != 'ROS45')
#cells <- which(metadata2$reclassify == 'Micro-PVM')
#cells <- which(metadata$reclassify == 'Astro')
#cells <- which(metadata$reclassify == 'Oligo')
cells <- which(metadata$reclassify == 'OPC')
counts2 <- counts[,cells]
micrometa2 <- metadata[cells,]
dim(counts2)
dim(micrometa2)


#upload differentially expressed genes curated from mathys supplementary table 2
p4 <- synapser::synGet('syn26136476')
mathys_DEG <- read.csv(p4$path, header=TRUE)
mathys_DEG2 <- unique(mathys_DEG$DEGgenes)
names(mathys_DEG2)[names(mathys_DEG2) == "mathys_DEG2"] <- "genes"


genes2 <- c()
for (gene in unique(c(as.vector(mathys_DEG$DEGgenes)))){
#for (gene in unique(c(as.vector(wilcox_DEgenes$genes)))){
  
  if (gene %in% rownames(counts2)){
    genes2 <- c(genes2,which(rownames(counts2)==gene))
  }
}

counts2 <- counts2[genes2,]
dim(counts2)

gene_short_name <- rownames(counts2)
dim(counts2)
gene_short_name <- as.data.frame(gene_short_name)
rownames(gene_short_name)<-gene_short_name$gene_short_name

temp <- counts2
temp2 <- micrometa2

rownames(temp) <- NULL
colnames(temp) <- NULL

####### For non-batch-corrected matrix, run RunMonocleTobit. For batch-corrected matrix, run RunMonocleTobit2 ###

#MonRun <- RunMonocleTobit2(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

MonRun <- RunMonocleTobit2(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroF_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroM_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroF_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroM_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_samplecorrected_tree.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
dev.off()


g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g

g<- plot_cell_trajectory(MonRun,color_by = "Subcluster",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="subcluster")
g
plot_cell_trajectory(MonRun,color_by = "Subcluster",show_branch_points=F,use_color_gradient = F,cell_size = 1)


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

#microglia females sample-corrected: root at state 6
MonRun <- orderCells(MonRun, root_state = 6)
#microglia males sample-corrected: root at state 3
MonRun <- orderCells(MonRun, root_state = 3)
#Astrocyes females: reverse 
MonRun <- orderCells(MonRun, reverse=TRUE)
#Astrocyes males: root at state 6
MonRun <- orderCells(MonRun, root_state = 6)
#Oligodendrocytes females: root at state 7
MonRun <- orderCells(MonRun, root_state = 7)
#Oligodendrocytes males: reverse
MonRun <- orderCells(MonRun, reverse=TRUE)
#OPCs females: root at state 3
MonRun <- orderCells(MonRun, root_state = 3)
#OPCs males: root at state 4
MonRun <- orderCells(MonRun, root_state = 4)

plot_cell_trajectory(MonRun,color_by = "Pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1)

plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)


#relabel states: mic female

#MicroF
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 1
MonRun$State2[MonRun$State == 4] <- 2
MonRun$State2[MonRun$State == 3] <- 2
MonRun$State2[MonRun$State == 2] <- 2
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 1] <- 4
MonRun$State2[MonRun$State == 5] <- 5


#MicroM
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 3] <- 1
MonRun$State2[MonRun$State == 4] <- 2
MonRun$State2[MonRun$State == 2] <- 3
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 1] <- 5


#AstroF
MonRun$State2 <- MonRun$State
#States are in order as is after reversal

#AstroM
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 1
MonRun$State2[MonRun$State == 7] <- 2
MonRun$State2[MonRun$State == 5] <- 3
MonRun$State2[MonRun$State == 4] <- 4
MonRun$State2[MonRun$State == 3] <- 3
MonRun$State2[MonRun$State == 2] <- 5
MonRun$State2[MonRun$State == 1] <- 5

#OligoF
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 7] <- 1
MonRun$State2[MonRun$State == 4] <- 2
MonRun$State2[MonRun$State == 6] <- 2
MonRun$State2[MonRun$State == 5] <- 3
MonRun$State2[MonRun$State == 3] <- 4
MonRun$State2[MonRun$State == 2] <- 4
MonRun$State2[MonRun$State == 1] <- 5

table(MonRun$State)
#OligoM
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 1] <- 1
MonRun$State2[MonRun$State == 2] <- 2
MonRun$State2[MonRun$State == 10] <- 2
MonRun$State2[MonRun$State == 11] <- 2
MonRun$State2[MonRun$State == 3] <- 2
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 9] <- 5
MonRun$State2[MonRun$State == 7] <- 6
MonRun$State2[MonRun$State == 8] <- 7

#OpcF
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 3] <- 1
MonRun$State2[MonRun$State == 1] <- 2
MonRun$State2[MonRun$State == 2] <- 3

#OpcM
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 4] <- 1
MonRun$State2[MonRun$State == 7] <- 2
MonRun$State2[MonRun$State == 6] <- 3
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 3] <- 4
MonRun$State2[MonRun$State == 2] <- 5
MonRun$State2[MonRun$State == 1] <- 6
 

MonRun$State2 <- as.character(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)

#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroF_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/MicroM_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroF_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/AstroM_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_state_tree.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_state_tree.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_state_tree.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()

#tiff(file='~/scRNAseq-subtype-mapping/figures/microF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/microM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/astroF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/astroM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
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

#tiff(file='~/scRNAseq-subtype-mapping/figures/microF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/microM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/astroF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/astroM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
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

#tiff(file='~/scRNAseq-subtype-mapping/figures/microF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/microM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/astroF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/astroM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OligoM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/OpcF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/OpcM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
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

#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AstroF_samplecorrected_pstimeStates.csv", row.names=FALSE)
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AstroM_samplecorrected_pstimeStates.csv", row.names=FALSE)

#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OligoF_samplecorrected_pstimeStates.csv", row.names=FALSE)
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OligoM_samplecorrected_pstimeStates.csv", row.names=FALSE)

#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OpcF_samplecorrected_pstimeStates.csv", row.names=FALSE)
write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/OpcM_samplecorrected_pstimeStates.csv", row.names=FALSE)


#run logistic regression comparing pseudotiem between cases and controls only
#casecontrolF <- subset(Fvariables, Fvariables$simpleDiagnosis=='AD'|Fvariables$simpleDiagnosis=='Cont')

Fvariables$diag3 <- ifelse(Fvariables$diag2=='AD', 1, 0)

summary(glm(diag3 ~ pseudotime_sc,Fvariables,family='binomial'))
summary(glm(Diagnosis ~ pseudotime_sc,Fvariables,family='binomial'))
#tiff(file='~/prot-lineage/figures/FEMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(Fvariables,aes(x=Diagnosis,
                             y=pseudotime_sc,
                             color=Diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:3])
g
#dev.off()

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




