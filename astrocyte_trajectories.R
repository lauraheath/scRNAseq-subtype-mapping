BiocManager::install("monocle")
library(monocle)
library(dplyr)
library(ggplot2)
synapser::synLogin()



#astrocyte data
#expression matrix:
p <- synapser::synGet('syn26148105')
counts <- readRDS(p$path)
counts <- as.matrix(counts)

#genes
p2 <- synapser::synGet('syn26148106')
gene_short_name <- readRDS(p2$path)

#metadata
p3 <- synapser::synGet('syn26148107')
metadata <- readRDS(p3$path)

#upload differentially expressed genes curated from mathys supplementary table 2
p4 <- synapser::synGet('syn26136476')
mathys_DEG <- read.csv(p4$path, header=TRUE)
mathys_DEG2 <- unique(mathys_DEG$DEGgenes)
#names(mathys_DEG2)[names(mathys_DEG2) == "DEGgenes"] <- "gene_short_name"

genes2 <- c()
for (gene in unique(c(as.vector(mathys_DEG$DEGgenes)))){
  
  
  if (gene %in% rownames(counts)){
    genes2 <- c(genes2,which(rownames(counts)==gene))
  }
}
counts2 <- counts[genes2,]

#delete genes with zero expression
counts2 <- as.data.frame(counts2)
counts2 <- counts2[rowSums(counts2[])>0,]
dim(counts2)
counts2 <- as.matrix(counts2)
gene_short_name <- rownames(counts2)
gene_short_name <- as.data.frame(gene_short_name)
rownames(gene_short_name)<-gene_short_name$gene_short_name



#female patients only
female <- which(metadata$sex == 'female')
counts3 <- counts2[,female]
meta2 <- metadata[female,]

#male patients only
# male <- which(metadata$sex == 'male')
# counts3 <- counts2[,male]
# meta2 <- metadata[male,]

#scale the counts in the matrix using ColNorm function
counts3 <- ColNorm(counts3)
temp <- counts3
temp2 <- meta2

rownames(temp) <- NULL
colnames(temp) <- NULL


MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g

g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g

MonRun$Pseudotime <- as.numeric(MonRun$Pseudotime)
plot_cell_trajectory(MonRun,color_by = "Pseudotime")


#root is in the wrong state, can't just reverse to correct. need to call the the correct state for the root
#determine state with highest proportion of control cells
table(MonRun$Diagnosis, MonRun$State)
df <- list()
df$State <- MonRun$State
df$Diagnosis <- MonRun$Diagnosis
df <- as.data.frame(df)
table(df$State)
proptable <- with(df, table(Diagnosis, State)) %>% prop.table(margin = 2)
proptable

#State 4 is more appropriate root (state 2 has higher proportion of control cells but only 7 total cells in state)
#Reorder cells and set root state = 4
MonRun <- orderCells(MonRun, root_state = 4)
plot_cell_trajectory(MonRun,color_by = "Pseudotime")
plot_cell_trajectory(MonRun,color_by = "State")

g<- plot_cell_trajectory(MonRun,color_by = "ros_ids",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

#relabel states so they are more intuitive (i.e. first state at the root, etc.), and combine small state with others due to small numbers
table(MonRun$State)
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 4] <- 1
MonRun$State2[MonRun$State == 5] <- 2
MonRun$State2[MonRun$State == 1] <- 3
MonRun$State2[MonRun$State == 3] <- 4

table(MonRun$State2)
MonRun$State2 <- as.character(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g

#tiff(file='~/prot-lineage/figures/MALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
MonRun$braaksc <- as.factor(MonRun$braaksc)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
#dev.off()

MonRun$cogdx <- as.factor(MonRun$cogdx)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="cogdx")
g


MonRun$ceradsc <- as.factor(MonRun$ceradsc)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="ceradsc")
g


#Save the monocle object:
saveRDS(MonRun, file="~/scRNAseq-subtype-mapping/female_astro_MonObject.RDS")
#save matrix for DE analysis:
saveRDS(counts3, file="~/scRNAseq-subtype-mapping/female_astro_matrix.RDS")


x <- list()
x$ros_ids <- MonRun$ros_ids
x$State <- MonRun$State
x$Pseudotime <- MonRun$Pseudotime
x$Diagnosis <- MonRun$Diagnosis
x$diag2 <- MonRun$simpleDiagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$apoe <- MonRun$apoe_genotype
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi

#females: rename and create a scaled pseudotime variable
Fvariables <- as.data.frame(x)
Fvariables$pseudotime_sc <- scale(Fvariables$Pseudotime, center=F)

#save variables file for later
#write.csv(prot_pstime_covars_F, file="~/prot-lineage/results/prot_pstime_covars_F.csv", row.names=FALSE)

#for males:
#prot_pstime_covars_M <- as.data.frame(x)
#prot_pstime_covars_M$pseudotime_sc <- scale(prot_pstime_covars_M$Pseudotime, center=F)
#write.csv(prot_pstime_covars_M, file="~/prot-lineage/results/prot_pstime_covars_M.csv", row.names=FALSE)



#run logistic regression comparing pseudotiem between cases and controls only
#casecontrolF <- subset(Fvariables, Fvariables$simpleDiagnosis=='AD'|Fvariables$simpleDiagnosis=='Cont')
Fvariables$diag3 <- ifelse(Fvariables$diag2=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,Fvariables,family='binomial'))
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




