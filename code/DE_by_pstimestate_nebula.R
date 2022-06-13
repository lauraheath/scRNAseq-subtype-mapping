library(nebula)


#upload mathys data with reclassified celltypes & metadata as a seurat object
p <- synapser::synGet('syn30070853')
seurat <- readRDS(p$path)



#need original, non-normalized counts matrix with integers (can upload Seurat object and export counts)
#seuratF <- readRDS(file='~/scRNAseq-subtype-mapping/data_objects/allcellsF_Seurat.RDS')
#seuratM <- readRDS(file='~/scRNAseq-subtype-mapping/data_objects/allcellsM_Seurat.RDS')


#seurat2 <- subset(seurat, reclassify=='Micro-PVM')
#seurat2 <- subset(seurat, reclassify=='Astro' & sex=='female')
seurat2 <- subset(seurat, reclassify=='Astro' & sex=='male')
counts <- GetAssayData(object=seurat2, slot="counts")
counts <- as.matrix(counts)
dim(counts)

#upload dataframe with pseudotimes and states from monocle2_micro_trajectories.R
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/MicroF_samplecorrected_pstimeStates.csv', row.names=F)
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/MicroM_samplecorrected_pstimeStates.csv', row.names=F)
#pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AstroF_samplecorrected_pstimeStates.csv')
pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AstroM_samplecorrected_pstimeStates.csv')



#order the metadata by subject ID, then reorder the counts matrix to match:
rownames(pstimes) <- pstimes$cellID
pstimes <- pstimes[order(pstimes$ros_ids),]
col.order <- pstimes$cellID
counts <- counts[,col.order]

#test <- as.data.frame(colnames(counts))

#restrict matrix to differentially expressed by AD genes used to construct trajectory:
# p4 <- synapser::synGet('syn26136476')
# mathys_DEG <- read.csv(p4$path, header=TRUE)
# mathys_DEG2 <- unique(mathys_DEG$DEGgenes)
# names(mathys_DEG2)[names(mathys_DEG2) == "mathys_DEG2"] <- "genes"
# 
# 
# genes2 <- c()
# for (gene in unique(c(as.vector(mathys_DEG$DEGgenes)))){
#   #for (gene in unique(c(as.vector(wilcox_DEgenes$genes)))){
#   
#   if (gene %in% rownames(counts)){
#     genes2 <- c(genes2,which(rownames(counts)==gene))
#   }
# }
# 
# counts2 <- counts[genes2,]
# dim(counts2)



###### first pass: differentially expressed genes with State 1 as the reference state #####
################for Microglia F and M: five total states, run nebula 4 times #####################
################for Astrocytes F and M: six total states, run nebula 5 times #####################

#need to run nebula for each State vs State 1
#first label all state1 as '1', rest is 0
table(pstimes$State)
pstimes$state <- pstimes$State
pstimes["state"][pstimes["state"] == 1] <- 0
pstimes["state"][pstimes["state"] > 1] <- 1
table(pstimes$state)

counts2 <- counts

###### State 2 v State 1 ##########
state1v2 <- which(pstimes$State==1|pstimes$State==2)
counts1v2 <- counts2[,state1v2]
df_1v2 <- pstimes[state1v2,]
table(df_1v2$State)
table(df_1v2$state)
#perform nebula on state1 vs state2
subjectID <- df_1v2$ros_ids
modmat = model.matrix(~state, data=df_1v2)
model1 <- nebula(counts1v2, subjectID, pred=modmat)
#model1
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
model1sum$state <- '2'
adjustp <- p.adjust(model1sum$p_state, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)

###### State 3 v State 1 ##########
state1v2 <- which(pstimes$State==1|pstimes$State==3)
counts1v2 <- counts2[,state1v2]
df_1v2 <- pstimes[state1v2,]
table(df_1v2$State)
table(df_1v2$state)
#perform nebula on state1 vs state2
subjectID <- df_1v2$ros_ids
modmat = model.matrix(~state, data=df_1v2)
model1 <- nebula(counts1v2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model2sum <- model1$summary
model2sum$state <- '3'
adjustp <- p.adjust(model2sum$p_state, "fdr")
model2sum <- cbind(adjustp, model2sum)


###### State 4 v State 1 ##########
state1v2 <- which(pstimes$State==1|pstimes$State==4)
counts1v2 <- counts2[,state1v2]
df_1v2 <- pstimes[state1v2,]
table(df_1v2$State)
table(df_1v2$state)
#perform nebula on state1 vs state2
subjectID <- df_1v2$ros_ids
modmat = model.matrix(~state, data=df_1v2)
model1 <- nebula(counts1v2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model3sum <- model1$summary
model3sum$state <- '4'
adjustp <- p.adjust(model3sum$p_state, "fdr")
model3sum <- cbind(adjustp, model3sum)


###### State 5 v State 1 ##########
state1v2 <- which(pstimes$State==1|pstimes$State==5)
counts1v2 <- counts2[,state1v2]
df_1v2 <- pstimes[state1v2,]
table(df_1v2$State)
table(df_1v2$state)
#perform nebula on state1 vs state2
subjectID <- df_1v2$ros_ids
modmat = model.matrix(~state, data=df_1v2)
model1 <- nebula(counts1v2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model4sum <- model1$summary
model4sum$state <- '5'
adjustp <- p.adjust(model4sum$p_state, "fdr")
model4sum <- cbind(adjustp, model4sum)

###### State 6 v State 1 ##########
state1v2 <- which(pstimes$State==1|pstimes$State==6)
counts1v2 <- counts2[,state1v2]
df_1v2 <- pstimes[state1v2,]
table(df_1v2$State)
table(df_1v2$state)
#perform nebula on state1 vs state2
subjectID <- df_1v2$ros_ids
modmat = model.matrix(~state, data=df_1v2)
model1 <- nebula(counts1v2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model5sum <- model1$summary
model5sum$state <- '6'
adjustp <- p.adjust(model5sum$p_state, "fdr")
model5sum <- cbind(adjustp, model5sum)

#now join the results dataframe and save:
results <- rbind(model1sum, model2sum, model3sum, model4sum, model5sum)


###### To characterize State 1 (reference state), run nebula for state 1 vs all other states
#reverse the 0/1 coding (so State 1 is 1, all others are reference)
table(pstimes$State)
pstimes$state2 <- pstimes$State
pstimes["state2"][pstimes["state2"] > 1] <- 0
table(pstimes$state2)

#perform nebula on state1 vs allstates
subjectID <- pstimes$ros_ids
modmat = model.matrix(~state2, data=pstimes)
model1 <- nebula(counts2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
model1sum$state <- '1'
adjustp <- p.adjust(model1sum$p_state, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)

####state 2 vs all others
table(pstimes$State)
pstimes$state2 <- pstimes$State
pstimes["state2"][pstimes["state2"] != 2] <- 0
pstimes["state2"][pstimes["state2"] == 2] <- 1
table(pstimes$state2)
#perform nebula on state1 vs allstates
subjectID <- pstimes$ros_ids
modmat = model.matrix(~state2, data=pstimes)
model1 <- nebula(counts2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model2sum <- model1$summary
model2sum$state <- '2'
adjustp <- p.adjust(model2sum$p_state, "fdr")
model2sum <- cbind(adjustp, model2sum)
head(model2sum)

####state 3 vs all others
table(pstimes$State)
pstimes$state2 <- pstimes$State
pstimes["state2"][pstimes["state2"] != 3] <- 0
pstimes["state2"][pstimes["state2"] == 3] <- 1
table(pstimes$state2)
#perform nebula on state1 vs allstates
subjectID <- pstimes$ros_ids
modmat = model.matrix(~state2, data=pstimes)
model1 <- nebula(counts2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model3sum <- model1$summary
model3sum$state <- '3'
adjustp <- p.adjust(model3sum$p_state, "fdr")
model3sum <- cbind(adjustp, model3sum)
head(model3sum)


####state 4 vs all others
table(pstimes$State)
pstimes$state2 <- pstimes$State
pstimes["state2"][pstimes["state2"] != 4] <- 0
pstimes["state2"][pstimes["state2"] == 4] <- 1
table(pstimes$state2)
#perform nebula on state1 vs allstates
subjectID <- pstimes$ros_ids
modmat = model.matrix(~state2, data=pstimes)
model1 <- nebula(counts2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model4sum <- model1$summary
model4sum$state <- '4'
adjustp <- p.adjust(model4sum$p_state, "fdr")
model4sum <- cbind(adjustp, model4sum)
head(model4sum)


####state 5 vs all others
table(pstimes$State)
pstimes$state2 <- pstimes$State
pstimes["state2"][pstimes["state2"] != 5] <- 0
pstimes["state2"][pstimes["state2"] == 5] <- 1
table(pstimes$state2)
#perform nebula on state1 vs allstates
subjectID <- pstimes$ros_ids
modmat = model.matrix(~state2, data=pstimes)
model1 <- nebula(counts2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model5sum <- model1$summary
model5sum$state <- '5'
adjustp <- p.adjust(model5sum$p_state, "fdr")
model5sum <- cbind(adjustp, model5sum)
head(model5sum)

####state 6 vs all others
table(pstimes$State)
pstimes$state2 <- pstimes$State
pstimes["state2"][pstimes["state2"] != 6] <- 0
pstimes["state2"][pstimes["state2"] == 6] <- 1
table(pstimes$state2)
#perform nebula on state1 vs allstates
subjectID <- pstimes$ros_ids
modmat = model.matrix(~state2, data=pstimes)
model1 <- nebula(counts2, subjectID, pred=modmat)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model6sum <- model1$summary
model6sum$state <- '6'
adjustp <- p.adjust(model6sum$p_state, "fdr")
model6sum <- cbind(adjustp, model6sum)
head(model6sum)


results2 <- rbind(model1sum, model2sum, model3sum, model4sum, model5sum, model6sum)


#write.csv(results, file="~/scRNAseq-subtype-mapping/DEresults/microgliaF_samplecorrected_DEresults.csv")
#write.csv(model1sum, file="~/scRNAseq-subtype-mapping/DEresults/microgliaF_samplecorrected_ReferenceResults.csv")

#write.csv(results, file="~/scRNAseq-subtype-mapping/DEresults/microgliaM_samplecorrected_DEresults.csv")
#write.csv(model1sum, file="~/scRNAseq-subtype-mapping/DEresults/microgliaM_samplecorrected_ReferenceResults.csv")

#write.csv(results, file="~/scRNAseq-subtype-mapping/DEresults/AstroF_pstimeDE_State1Reference.csv", row.names=FALSE)
#write.csv(results2, file="~/scRNAseq-subtype-mapping/DEresults/AstroF_pstimeDE_eachState_vs_allstates.csv", row.names=FALSE)

write.csv(results, file="~/scRNAseq-subtype-mapping/DEresults/AstroM_pstimeDE_State1Reference.csv", row.names=FALSE)
write.csv(results2, file="~/scRNAseq-subtype-mapping/DEresults/AstroM_pstimeDE_eachState_vs_allstates.csv", row.names=FALSE)

