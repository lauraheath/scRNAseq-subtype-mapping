library(dplyr)
library(ggplot2)
library(broom)
library(pROC)
library(caret)
install.packages("e1071")
library(tidyverse)
library(viridis)

#upload mathys data with reclassified celltypes & metadata as a seurat object
p <- synapser::synGet('syn30070853')
seurat <- readRDS(p$path)



#need original, non-normalized counts matrix with integers (can upload Seurat object and export counts)
#seuratF <- readRDS(file='~/scRNAseq-subtype-mapping/data_objects/allcellsF_Seurat.RDS')
#seuratM <- readRDS(file='~/scRNAseq-subtype-mapping/data_objects/allcellsM_Seurat.RDS')


#seurat2 <- subset(seurat, reclassify=='Micro-PVM')
#seurat2 <- subset(seurat, reclassify=='Astro' & sex=='female')
# seurat2 <- subset(seurat, reclassify=='Astro' & sex=='male')
# counts <- GetAssayData(object=seurat2, slot="counts")
# counts <- as.matrix(counts)
# dim(counts)

#upload dataframe with pseudotimes and states from monocle2_micro_trajectories.R
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/MicroF_samplecorrected_pstimeStates.csv', row.names=F)
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/MicroM_samplecorrected_pstimeStates.csv', row.names=F)
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AstroF_samplecorrected_pstimeStates.csv')
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AstroM_samplecorrected_pstimeStates.csv')

pstimes <- read.csv(file="~/scRNAseq-subtype-mapping/data_objects/AstroF_allgenes_pstimeStates.csv")

#order the metadata by subject ID, then reorder the counts matrix to match:
# rownames(pstimes) <- pstimes$cellID
# pstimes <- pstimes[order(pstimes$ros_ids),]
# col.order <- pstimes$cellID
# counts <- counts[,col.order]



#percentage of broad celltypes by disease state
# allmeta <- seurat@meta.data
# df1 <- allmeta %>% group_by(Diagnosis, reclassify) %>%
#   summarise(Nb = n()) %>%
#   mutate(C = sum(Nb)) %>%
#   mutate(percent = Nb/C*100)
# names(df1)[names(df1) == "reclassify"] <- "broad_cell_type"
# df1 <- as.data.frame(df1)
# 
# 
# ggplot(df1, aes(fill = Diagnosis, y = percent, x = broad_cell_type))+
#   geom_bar(position = "dodge", stat = "identity")+ theme_classic()
# 
# #are there sex differences?
# allmetaF <- subset(allmeta, allmeta$sex=='female')
# df1 <- allmetaF %>% group_by(Diagnosis, reclassify) %>%
#   summarise(Nb = n()) %>%
#   mutate(C = sum(Nb)) %>%
#   mutate(percent = Nb/C*100)
# names(df1)[names(df1) == "reclassify"] <- "broad_cell_type"
# df1 <- as.data.frame(df1)
# female_p <- ggplot(df1, aes(fill = Diagnosis, y = percent, x = broad_cell_type))+
#   geom_bar(position = "dodge", stat = "identity")+ theme_classic()
# female_p
# 
# allmetaM <- subset(allmeta, allmeta$sex=='male')
# df2 <- allmetaM %>% group_by(Diagnosis, reclassify) %>%
#   summarise(Nb = n()) %>%
#   mutate(C = sum(Nb)) %>%
#   mutate(percent = Nb/C*100)
# names(df2)[names(df2) == "reclassify"] <- "broad_cell_type"
# df2 <- as.data.frame(df2)
# male_p <- ggplot(df2, aes(fill = Diagnosis, y = percent, x = broad_cell_type))+
#   geom_bar(position = "dodge", stat = "identity")+ theme_classic()
# male_p



#distribution of pseudotimes by patient/disease state
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AstroF_samplecorrected_pstimeStates.csv')
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AstroM_samplecorrected_pstimeStates.csv')

ggplot(data=pstimes, aes(x=Pseudotime)) + geom_histogram(binwidth=1)


control <- subset(pstimes, pstimes$Diagnosis=='Control')
p1 <- ggplot(data=control, aes(x=Pseudotime)) + geom_histogram(binwidth=1)
p1
p1 + facet_wrap(~ros_ids)

early <- subset(pstimes, pstimes$Diagnosis=='Early-pathology AD')
p1 <- ggplot(data=early, aes(x=Pseudotime)) + geom_histogram(binwidth=1)
p1
p1 + facet_wrap(~ros_ids)

late <- subset(pstimes, pstimes$Diagnosis=='Late-pathology AD')
p1 <- ggplot(data=late, aes(x=Pseudotime)) + geom_histogram(binwidth=1)
p1
p1 + facet_wrap(~ros_ids)




#stacked barchart by patient:
#need a data frame with one row per patient, state, #of cells per state, and diagnosis
#start with the pstime dataframe

df <- subset(pstimes, select=c(cellID, ros_ids, State, Diagnosis))

statecounts <- as.data.frame(table(df$ros_ids, df$State))
#get diagnosis data for each ros_id
diagnos <- as.data.frame(table(pstimes$ros_ids, pstimes$Diagnosis))
diagnos <- diagnos[diagnos$Freq!=0, ]
diagnos$Freq<-NULL
names(diagnos)[names(diagnos) == "Var2"] <- "Diagnosis"
statecounts2 <- merge(statecounts, diagnos, by='Var1')
names(statecounts2)[names(statecounts2) == "Var1"] <- "ros_ids"
names(statecounts2)[names(statecounts2) == "Var2"] <- "State"
names(statecounts2)[names(statecounts2) == "Freq"] <- "Cells"


# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(statecounts2$State))
statecounts2$id <- rep( seq(1, nrow(statecounts2)/nObsType) , each=nObsType)
label_data <- statecounts2 %>% group_by(id, ros_ids) %>% summarize(tot=sum(Cells))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- statecounts2 %>% 
  group_by(Diagnosis) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

#draw a stacked barchart by diagnosis
p <- ggplot(statecounts2) +    
  geom_bar(aes(x=as.factor(ros_ids), y=Cells, fill=State), stat="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  facet_grid(~Diagnosis, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90))

p

#plot distributions of pseudotime in box plots by patient
p1 <- ggplot(pstimes, aes(x=ros_ids, y=Pseudotime, fill=Diagnosis)) +
  geom_boxplot(outlier.size=0.5) +
  facet_wrap(~Diagnosis, scale="free") +
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_point(size=0.5, position=ggplot2::position_jitterdodge(), color = "dimgray")
p1




