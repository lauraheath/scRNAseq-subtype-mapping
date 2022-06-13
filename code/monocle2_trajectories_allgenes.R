install.packages('lme4')

library(monocle)
library(dplyr)
library(ggplot2)
library(lme4)
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
cells <- which(metadata$reclassify == 'Astro')
#cells <- which(metadata$reclassify == 'Oligo')
#cells <- which(metadata$reclassify == 'OPC')
counts2 <- counts[,cells]
micrometa2 <- metadata2[cells,]
dim(counts2)
dim(micrometa2)


gene_short_name <- rownames(counts2)
dim(counts2)
gene_short_name <- as.data.frame(gene_short_name)
rownames(gene_short_name)<-gene_short_name$gene_short_name

temp <- counts2
temp2 <- micrometa2

rownames(temp) <- NULL
colnames(temp) <- NULL

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
