library(dplyr)
library(ggplot2)
library(broom)
library(pROC)
library(caret)
install.packages("e1071")

mathys_meta3 <- mathys3@meta.data

df1 <- mathys_meta3 %>% group_by(Diagnosis, predicted.subclass_label) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)
names(df1)[names(df1) == "predicted.subclass_label"] <- "subclass_label"
df1 <- as.data.frame(df1)

df2 <- metadata %>% group_by(subclass_label) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)
df2$Diagnosis = "M1reference"
col_order <- c("Diagnosis", "subclass_label", "Nb", "C", "percent")
df2 <- df2[,col_order]
df2 <- as.data.frame(df2)

df3 <- rbind(df1, df2)


ggplot(df3, aes(fill = Diagnosis, y = percent, x = subclass_label))+
  geom_bar(position = "dodge", stat = "identity")+ theme_classic()

#no maps to L5 ET and Sst Chodl subclasses in cont, early, or late. add those classes with zeroes
Diagnosis <- c("Cont", "Cont", "Early", "Early", "Late", "Late")
subclass_label <- c("L5 ET", "Sst Chodl", "L5 ET", "Sst Chodl", "L5 ET", "Sst Chodl")
Nb <- c(0, 0, 0, 0, 0, 0)
C <- c(35140, 35140, 23072, 23072, 12422, 12422)
percent <- c(0,0,0,0,0,0)
addons <- data.frame(Diagnosis, subclass_label, Nb, C, percent)
df3 <- rbind(df3, addons)

ggplot(df3, aes(fill = Diagnosis, y = percent, x = subclass_label))+
  geom_bar(position = "dodge", stat = "identity")+ theme_classic()



#create a confusion matrix by comparing the mathys broad cell type calls to the predictions. First recode the subclass labels as broad cell types:

mathys_meta3$predicted.broadcelltypes <- ifelse(mathys_meta3$predicted.subclass_label=='Astro', 'Ast',
                                                ifelse(mathys_meta3$predicted.subclass_label=='Endo', 'End',
                                                       ifelse(mathys_meta3$predicted.subclass_label=='L2/3 IT', 'Ex',
                                                              ifelse(mathys_meta3$predicted.subclass_label=='L5 ET', 'Ex',
                                                                     ifelse(mathys_meta3$predicted.subclass_label=='L5 IT', 'Ex',
                                                                            ifelse(mathys_meta3$predicted.subclass_label=='L5/6 NP', 'Ex',
                                                                                   ifelse(mathys_meta3$predicted.subclass_label=='L6 CT', 'Ex',
                                                                                          ifelse(mathys_meta3$predicted.subclass_label=='L6 IT', 'Ex',
                                                                                                 ifelse(mathys_meta3$predicted.subclass_label=='L6 IT Car3', 'Ex',
                                                                                                        ifelse(mathys_meta3$predicted.subclass_label=='L6b', 'Ex',
                                                                                                               ifelse(mathys_meta3$predicted.subclass_label=='Lamp5', 'In',
                                                                                                                      ifelse(mathys_meta3$predicted.subclass_label=='Micro-PVM', 'Mic',
                                                                                                                             ifelse(mathys_meta3$predicted.subclass_label=='Oligo', 'Oli',
                                                                                                                                    ifelse(mathys_meta3$predicted.subclass_label=='OPC', 'Opc',
                                                                                                                                           ifelse(mathys_meta3$predicted.subclass_label=='Pvalb','In',
                                                                                                                                                  ifelse(mathys_meta3$predicted.subclass_label=='Sncg', 'In',
                                                                                                                                                         ifelse(mathys_meta3$predicted.subclass_label=='Sst', 'In',
                                                                                                                                                                ifelse(mathys_meta3$predicted.subclass_label=='Sst Chodl', 'In',
                                                                                                                                                                       ifelse(mathys_meta3$predicted.subclass_label=='Vip', 'In', 'Per')))))))))))))))))))
table(mathys_meta3$predicted.broadcelltypes)
table(mathys_meta3$predicted.broadcelltypes, mathys_meta3$broad.cell.type)

mathys_meta3$broad.cell.type <- as.factor(mathys_meta3$broad.cell.type)
mathys_meta3$predicted.broadcelltypes <- as.factor(mathys_meta3$predicted.broadcelltypes)
conf <- confusionMatrix(mathys_meta3$predicted.broadcelltypes, mathys_meta3$broad.cell.type)
conf


#create heatmaps of the confusion matrices:
conf_mat <- as.matrix(conf)
conf_mat <- t(conf_mat)
heatmap(conf_mat, Rowv = NA, Colv = NA, scale="column")
library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(9, "Blues"))(100)
heatmap(conf_mat, Rowv = NA, Colv = NA, col = cols, scale="column")
heatmap(conf_mat, Rowv = NA, Colv = NA, col = cols, scale="row")
heatmap(conf_mat, Rowv = NA, Colv = NA, col = cols)





#create a misclassification variable to flag the cells that were misclassified:
mathys_meta3$misclass <- ifelse(mathys_meta3$predicted.broadcelltypes!=mathys_meta3$broad.cell.type, 1, 0)
table(mathys_meta3$misclass)

misclass <- subset(mathys_meta3, mathys_meta3$misclass==1)
table(mathys_meta3$Diagnosis, mathys_meta3$misclass)
table(misclass$predicted.subclass_label)
table(misclass$predicted.broadcelltypes)
table(misclass$broad.cell.type)
mean(misclass$predicted.subclass_label.score)
not_misclass <- subset(mathys_meta3, mathys_meta3$misclass==0)
mean(not_misclass$predicted.subclass_label.score)

#statistical testing to compare prediction scores
cont_early <- subset(mathys_meta3, mathys_meta3$Diagnosis!="Late")
cont_late <- subset(mathys_meta3, mathys_meta3$Diagnosis!="Early")
early_late <- subset(mathys_meta3, mathys_meta3$Diagnosis!="Cont")

Astro1 <- subset(cont_early, cont_early$predicted.subclass_label=='Astro')
Endo1 <- subset(cont_early, cont_early$predicted.subclass_label=='Endo')
L23_IT1 <- subset(cont_early, cont_early$predicted.subclass_label=='L2/3 IT')
L5_ET1 <- subset(cont_early, cont_early$predicted.subclass_label=='L5 ET')
L5_IT1 <- subset(cont_early, cont_early$predicted.subclass_label=='L5 IT')
L56_NP1 <- subset(cont_early, cont_early$predicted.subclass_label=='L5/6 NP')
L6_CT1 <- subset(cont_early, cont_early$predicted.subclass_label=='L6 CT')
L6_IT1 <- subset(cont_early, cont_early$predicted.subclass_label=='L6 IT')
L6_ITCar31 <- subset(cont_early, cont_early$predicted.subclass_label=='L6 IT Car3')
L6b1 <- subset(cont_early, cont_early$predicted.subclass_label=='L6b')
Lamp51 <- subset(cont_early, cont_early$predicted.subclass_label=='Lamp5')
Micro_PVM1 <- subset(cont_early, cont_early$predicted.subclass_label=='Micro-PVM')
Oligo1 <- subset(cont_early, cont_early$predicted.subclass_label=='Oligo')
OPC1 <- subset(cont_early, cont_early$predicted.subclass_label=='OPC')
Pvalb1 <- subset(cont_early, cont_early$predicted.subclass_label=='Pvalb')
Sncg1 <- subset(cont_early, cont_early$predicted.subclass_label=='Sncg')
Sst1 <- subset(cont_early, cont_early$predicted.subclass_label=='Sst')
Sst_Chodl1 <- subset(cont_early, cont_early$predicted.subclass_label=='Sst Chodl')
Vip1 <- subset(cont_early, cont_early$predicted.subclass_label=='Vip')
VLMC1 <- subset(cont_early, cont_early$predicted.subclass_label=='VLMC')

stats <- 0
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Astro1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Astro1'
stats <- df
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Endo1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Endo1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L23_IT1)
df <- as.data.frame(df$p.value)
df$celltype <- 'L23_IT1'
stats <- rbind(stats, df)
#df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L5_ET1)
#df <- as.data.frame(df$p.value)
#df$celltype <- 'L5_ET1'
#stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L5_IT1)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5_IT1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L56_NP1)
df <- as.data.frame(df$p.value)
df$celltype <- 'L56_NP1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_CT1)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_CT1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_IT1 )
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_IT1 '
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_ITCar31)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_ITCar31'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6b1)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6b1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Lamp51)
df <- as.data.frame(df$p.value)
df$celltype <- 'Lamp51'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Micro_PVM1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Micro_PVM1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Oligo1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Oligo1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=OPC1)
df <- as.data.frame(df$p.value)
df$celltype <- 'OPC1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Pvalb1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Pvalb1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sncg1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sncg1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sst1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sst1'
stats <- rbind(stats, df)
#df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sst_Chodl1)
#df <- as.data.frame(df$p.value)
#df$celltype <- 'Sst_Chodl1'
#stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Vip1)
df <- as.data.frame(df$p.value)
df$celltype <- 'Vip1'
stats <- rbind(stats, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=VLMC1)
df <- as.data.frame(df$p.value)
df$celltype <- 'VLMC1'
stats <- rbind(stats, df)

ggplot(L5_IT1, aes(y = predicted.subclass_label.score, x = predicted.subclass_label, fill = Diagnosis))+
  geom_violin(trim=TRUE)

#control vs late
Astro2 <- subset(cont_late, cont_late$predicted.subclass_label=='Astro')
Endo2 <- subset(cont_late, cont_late$predicted.subclass_label=='Endo')
L23_IT2 <- subset(cont_late, cont_late$predicted.subclass_label=='L2/3 IT')
L5_ET2 <- subset(cont_late, cont_late$predicted.subclass_label=='L5 ET')
L5_IT2 <- subset(cont_late, cont_late$predicted.subclass_label=='L5 IT')
L56_NP2 <- subset(cont_late, cont_late$predicted.subclass_label=='L5/6 NP')
L6_CT2 <- subset(cont_late, cont_late$predicted.subclass_label=='L6 CT')
L6_IT2 <- subset(cont_late, cont_late$predicted.subclass_label=='L6 IT')
L6_ITCar32 <- subset(cont_late, cont_late$predicted.subclass_label=='L6 IT Car3')
L6b2 <- subset(cont_late, cont_late$predicted.subclass_label=='L6b')
Lamp52 <- subset(cont_late, cont_late$predicted.subclass_label=='Lamp5')
Micro_PVM2 <- subset(cont_late, cont_late$predicted.subclass_label=='Micro-PVM')
Oligo2 <- subset(cont_late, cont_late$predicted.subclass_label=='Oligo')
OPC2 <- subset(cont_late, cont_late$predicted.subclass_label=='OPC')
Pvalb2 <- subset(cont_late, cont_late$predicted.subclass_label=='Pvalb')
Sncg2 <- subset(cont_late, cont_late$predicted.subclass_label=='Sncg')
Sst2 <- subset(cont_late, cont_late$predicted.subclass_label=='Sst')
Sst_Chodl2 <- subset(cont_late, cont_late$predicted.subclass_label=='Sst Chodl')
Vip2 <- subset(cont_late, cont_late$predicted.subclass_label=='Vip')
VLMC2 <- subset(cont_late, cont_late$predicted.subclass_label=='VLMC')

stats2 <- 0
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Astro2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Astro2'
stats2 <- df
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Endo2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Endo2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L23_IT2)
df <- as.data.frame(df$p.value)
df$celltype <- 'L23_IT2'
stats2 <- rbind(stats2, df)
# df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L5_ET2)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'L5_ET2'
# stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L5_IT2)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5_IT2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L56_NP2)
df <- as.data.frame(df$p.value)
df$celltype <- 'L56_NP2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_CT2)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_CT2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_IT2 )
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_IT2 '
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_ITCar32)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_ITCar32'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6b2)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6b2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Lamp52)
df <- as.data.frame(df$p.value)
df$celltype <- 'Lamp52'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Micro_PVM2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Micro_PVM2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Oligo2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Oligo2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=OPC2)
df <- as.data.frame(df$p.value)
df$celltype <- 'OPC2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Pvalb2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Pvalb2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sncg2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sncg2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sst2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sst2'
stats2 <- rbind(stats2, df)
# df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sst_Chodl2)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'Sst_Chodl2'
# stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Vip2)
df <- as.data.frame(df$p.value)
df$celltype <- 'Vip2'
stats2 <- rbind(stats2, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=VLMC2)
df <- as.data.frame(df$p.value)
df$celltype <- 'VLMC2'
stats2 <- rbind(stats2, df)


ggplot(Lamp52, aes(y = predicted.subclass_label.score, x = predicted.subclass_label, fill = Diagnosis))+
  geom_violin(trim=TRUE)


### early vs. late
Astro3 <- subset(early_late, early_late$predicted.subclass_label=='Astro')
Endo3 <- subset(early_late, early_late$predicted.subclass_label=='Endo')
L23_IT3 <- subset(early_late, early_late$predicted.subclass_label=='L2/3 IT')
L5_ET3 <- subset(early_late, early_late$predicted.subclass_label=='L5 ET')
L5_IT3 <- subset(early_late, early_late$predicted.subclass_label=='L5 IT')
L56_NP3 <- subset(early_late, early_late$predicted.subclass_label=='L5/6 NP')
L6_CT3 <- subset(early_late, early_late$predicted.subclass_label=='L6 CT')
L6_IT3 <- subset(early_late, early_late$predicted.subclass_label=='L6 IT')
L6_ITCar33 <- subset(early_late, early_late$predicted.subclass_label=='L6 IT Car3')
L6b3 <- subset(early_late, early_late$predicted.subclass_label=='L6b')
Lamp53 <- subset(early_late, early_late$predicted.subclass_label=='Lamp5')
Micro_PVM3 <- subset(early_late, early_late$predicted.subclass_label=='Micro-PVM')
Oligo3 <- subset(early_late, early_late$predicted.subclass_label=='Oligo')
OPC3 <- subset(early_late, early_late$predicted.subclass_label=='OPC')
Pvalb3 <- subset(early_late, early_late$predicted.subclass_label=='Pvalb')
Sncg3 <- subset(early_late, early_late$predicted.subclass_label=='Sncg')
Sst3 <- subset(early_late, early_late$predicted.subclass_label=='Sst')
Sst_Chodl3 <- subset(early_late, early_late$predicted.subclass_label=='Sst Chodl')
Vip3 <- subset(early_late, early_late$predicted.subclass_label=='Vip')
VLMC3 <- subset(early_late, early_late$predicted.subclass_label=='VLMC')



stats3 <- 0
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Astro3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Astro3'
stats3 <- df
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Endo3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Endo3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L23_IT3)
df <- as.data.frame(df$p.value)
df$celltype <- 'L23_IT3'
stats3 <- rbind(stats3, df)
# df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L5_ET3)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'L5_ET3'
# stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L5_IT3)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5_IT3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L56_NP3)
df <- as.data.frame(df$p.value)
df$celltype <- 'L56_NP3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_CT3)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_CT3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_IT2 )
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_IT2 '
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6_ITCar33)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6_ITCar33'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=L6b3)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6b3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Lamp53)
df <- as.data.frame(df$p.value)
df$celltype <- 'Lamp53'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Micro_PVM3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Micro_PVM3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Oligo3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Oligo3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=OPC3)
df <- as.data.frame(df$p.value)
df$celltype <- 'OPC3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Pvalb3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Pvalb3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sncg3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sncg3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sst3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sst3'
stats3 <- rbind(stats3, df)
# df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Sst_Chodl3)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'Sst_Chodl3'
# stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=Vip3)
df <- as.data.frame(df$p.value)
df$celltype <- 'Vip3'
stats3 <- rbind(stats3, df)
df <- kruskal.test(predicted.subclass_label.score ~ Diagnosis, data=VLMC3)
df <- as.data.frame(df$p.value)
df$celltype <- 'VLMC3'
stats3 <- rbind(stats3, df)

stats$comparison <- 'cont_early'
stats2$comparison <- 'cont_late'
stats3$comparison <- 'early_late'
KWstats <- rbind(stats, stats2)
KWstats <- rbind(KWstats, stats3)
names(KWstats)[names(KWstats) == "df$p.value"] <- "pvalue"
head(KWstats)
KWstats$fdr <- p.adjust(KWstats$pvalue, method="fdr")
write.csv(KWstats, file="~/celltype_mapping/KWpredscores_stats.csv")


#chisquare test for proportions
dat1 <- cont_early
dat1$diag01 <- ifelse(dat1$simpleDiagnosis=='Cont', 0, 1)

dat1$astro01 <- ifelse(dat1$predicted.subclass_label=='Astro', 0, 1)
dat1$endo01 <- ifelse(dat1$predicted.subclass_label=='Endo', 0, 1)
dat1$L23it01 <- ifelse(dat1$predicted.subclass_label=='L2/3 IT', 0, 1)
dat1$L5et01 <- ifelse(dat1$predicted.subclass_label=='L5 ET', 0, 1)
dat1$L5it01 <- ifelse(dat1$predicted.subclass_label=='L5 IT', 0, 1)
dat1$L56np01 <- ifelse(dat1$predicted.subclass_label=='L5/6 NP', 0, 1)
dat1$L6ct01 <- ifelse(dat1$predicted.subclass_label=='L6 CT', 0, 1)
dat1$L6it01 <- ifelse(dat1$predicted.subclass_label=='L6 IT', 0, 1)
dat1$L6itcar301 <- ifelse(dat1$predicted.subclass_label=='L6 IT Car3', 0, 1)
dat1$L6b01 <- ifelse(dat1$predicted.subclass_label=='L6b', 0, 1)
dat1$lamp501 <- ifelse(dat1$predicted.subclass_label=='Lamp5', 0, 1)
dat1$micro01 <- ifelse(dat1$predicted.subclass_label=='Micro-PVM', 0, 1)
dat1$oligo01 <- ifelse(dat1$predicted.subclass_label=='Oligo', 0, 1)
dat1$opc01 <- ifelse(dat1$predicted.subclass_label=='OPC', 0, 1)
dat1$pvalb01 <- ifelse(dat1$predicted.subclass_label=='Pvalb', 0, 1)
dat1$sncg01 <- ifelse(dat1$predicted.subclass_label=='Sncg', 0, 1)
dat1$sst01 <- ifelse(dat1$predicted.subclass_label=='Sst', 0, 1)
dat1$sstchodl01 <- ifelse(dat1$predicted.subclass_label=='Sst Chodl', 0, 1)
dat1$vip01 <- ifelse(dat1$predicted.subclass_label=='Vip', 0, 1)
dat1$vlmc01 <- ifelse(dat1$predicted.subclass_label=='VLMC', 0, 1)

stats<-0
df <- chisq.test(dat1$astro01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Astro'
stats <- df
df <- chisq.test(dat1$endo01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Endo'
stats <- rbind(stats, df)
df <- chisq.test(dat1$L23it01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L2/3 IT'
stats <- rbind(stats, df)
# df <- chisq.test(dat1$L5et01, dat1$diag01)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'L5 ET'
# stats <- rbind(stats, df)
df <- chisq.test(dat1$L5it01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5 IT'
stats <- rbind(stats, df)
df <- chisq.test(dat1$L56np01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5/6 NP'
stats <- rbind(stats, df)
df <- chisq.test(dat1$L6ct01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 CT'
stats <- rbind(stats, df)
df <- chisq.test(dat1$L6it01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 IT'
stats <- rbind(stats, df)
df <- chisq.test(dat1$L6itcar301, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 IT Car3'
stats <- rbind(stats, df)
df <- chisq.test(dat1$L6b01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6b'
stats <- rbind(stats, df)
df <- chisq.test(dat1$lamp501, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Lamp5'
stats <- rbind(stats, df)
df <- chisq.test(dat1$micro01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Micro PVM'
stats <- rbind(stats, df)
df <- chisq.test(dat1$oligo01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Oligo'
stats <- rbind(stats, df)
df <- chisq.test(dat1$opc01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'OPC'
stats <- rbind(stats, df)
df <- chisq.test(dat1$pvalb01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Pvalb'
stats <- rbind(stats, df)
df <- chisq.test(dat1$sncg01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sncg'
stats <- rbind(stats, df)
df <- chisq.test(dat1$sst01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sst'
stats <- rbind(stats, df)
# df <- chisq.test(dat1$sstchodl01, dat1$diag01)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'Sst Chodl'
# stats <- rbind(stats, df)
df <- chisq.test(dat1$vip01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Vip'
stats <- rbind(stats, df)
df <- chisq.test(dat1$vlmc01, dat1$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'VLMC'
stats <- rbind(stats, df)



#control vs late
dat2 <- cont_late
dat2$diag01 <- ifelse(dat2$simpleDiagnosis=='Cont', 0, 1)

dat2$astro01 <- ifelse(dat2$predicted.subclass_label=='Astro', 0, 1)
dat2$endo01 <- ifelse(dat2$predicted.subclass_label=='Endo', 0, 1)
dat2$L23it01 <- ifelse(dat2$predicted.subclass_label=='L2/3 IT', 0, 1)
dat2$L5et01 <- ifelse(dat2$predicted.subclass_label=='L5 ET', 0, 1)
dat2$L5it01 <- ifelse(dat2$predicted.subclass_label=='L5 IT', 0, 1)
dat2$L56np01 <- ifelse(dat2$predicted.subclass_label=='L5/6 NP', 0, 1)
dat2$L6ct01 <- ifelse(dat2$predicted.subclass_label=='L6 CT', 0, 1)
dat2$L6it01 <- ifelse(dat2$predicted.subclass_label=='L6 IT', 0, 1)
dat2$L6itcar301 <- ifelse(dat2$predicted.subclass_label=='L6 IT Car3', 0, 1)
dat2$L6b01 <- ifelse(dat2$predicted.subclass_label=='L6b', 0, 1)
dat2$lamp501 <- ifelse(dat2$predicted.subclass_label=='Lamp5', 0, 1)
dat2$micro01 <- ifelse(dat2$predicted.subclass_label=='Micro-PVM', 0, 1)
dat2$oligo01 <- ifelse(dat2$predicted.subclass_label=='Oligo', 0, 1)
dat2$opc01 <- ifelse(dat2$predicted.subclass_label=='OPC', 0, 1)
dat2$pvalb01 <- ifelse(dat2$predicted.subclass_label=='Pvalb', 0, 1)
dat2$sncg01 <- ifelse(dat2$predicted.subclass_label=='Sncg', 0, 1)
dat2$sst01 <- ifelse(dat2$predicted.subclass_label=='Sst', 0, 1)
dat2$sstchodl01 <- ifelse(dat2$predicted.subclass_label=='Sst Chodl', 0, 1)
dat2$vip01 <- ifelse(dat2$predicted.subclass_label=='Vip', 0, 1)
dat2$vlmc01 <- ifelse(dat2$predicted.subclass_label=='VLMC', 0, 1)

stats2<-0
df <- chisq.test(dat2$astro01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Astro'
stats2 <- df
df <- chisq.test(dat2$endo01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Endo'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$L23it01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L2/3 IT'
stats2 <- rbind(stats2, df)
# df <- chisq.test(dat2$L5et01, dat2$diag01)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'L5 ET'
# stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$L5it01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5 IT'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$L56np01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5/6 NP'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$L6ct01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 CT'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$L6it01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 IT'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$L6itcar301, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 IT Car3'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$L6b01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6b'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$lamp501, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Lamp5'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$micro01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Micro PVM'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$oligo01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Oligo'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$opc01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'OPC'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$pvalb01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Pvalb'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$sncg01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sncg'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$sst01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sst'
stats2 <- rbind(stats2, df)
# df <- chisq.test(dat2$sstchodl01, dat2$diag01)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'Sst Chodl'
# stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$vip01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Vip'
stats2 <- rbind(stats2, df)
df <- chisq.test(dat2$vlmc01, dat2$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'VLMC'
stats2 <- rbind(stats2, df)


#early vs. late
dat3 <- early_late
dat3$diag01 <- ifelse(dat3$Diagnosis=='Early', 0, 1)

dat3$astro01 <- ifelse(dat3$predicted.subclass_label=='Astro', 0, 1)
dat3$endo01 <- ifelse(dat3$predicted.subclass_label=='Endo', 0, 1)
dat3$L23it01 <- ifelse(dat3$predicted.subclass_label=='L2/3 IT', 0, 1)
dat3$L5et01 <- ifelse(dat3$predicted.subclass_label=='L5 ET', 0, 1)
dat3$L5it01 <- ifelse(dat3$predicted.subclass_label=='L5 IT', 0, 1)
dat3$L56np01 <- ifelse(dat3$predicted.subclass_label=='L5/6 NP', 0, 1)
dat3$L6ct01 <- ifelse(dat3$predicted.subclass_label=='L6 CT', 0, 1)
dat3$L6it01 <- ifelse(dat3$predicted.subclass_label=='L6 IT', 0, 1)
dat3$L6itcar301 <- ifelse(dat3$predicted.subclass_label=='L6 IT Car3', 0, 1)
dat3$L6b01 <- ifelse(dat3$predicted.subclass_label=='L6b', 0, 1)
dat3$lamp501 <- ifelse(dat3$predicted.subclass_label=='Lamp5', 0, 1)
dat3$micro01 <- ifelse(dat3$predicted.subclass_label=='Micro-PVM', 0, 1)
dat3$oligo01 <- ifelse(dat3$predicted.subclass_label=='Oligo', 0, 1)
dat3$opc01 <- ifelse(dat3$predicted.subclass_label=='OPC', 0, 1)
dat3$pvalb01 <- ifelse(dat3$predicted.subclass_label=='Pvalb', 0, 1)
dat3$sncg01 <- ifelse(dat3$predicted.subclass_label=='Sncg', 0, 1)
dat3$sst01 <- ifelse(dat3$predicted.subclass_label=='Sst', 0, 1)
dat3$sstchodl01 <- ifelse(dat3$predicted.subclass_label=='Sst Chodl', 0, 1)
dat3$vip01 <- ifelse(dat3$predicted.subclass_label=='Vip', 0, 1)
dat3$vlmc01 <- ifelse(dat3$predicted.subclass_label=='VLMC', 0, 1)

stats3<-0
df <- chisq.test(dat3$astro01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Astro'
stats3 <- df
df <- chisq.test(dat3$endo01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Endo'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$L23it01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L2/3 IT'
stats3 <- rbind(stats3, df)
# df <- chisq.test(dat3$L5et01, dat3$diag01)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'L5 ET'
# stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$L5it01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5 IT'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$L56np01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L5/6 NP'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$L6ct01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 CT'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$L6it01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 IT'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$L6itcar301, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6 IT Car3'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$L6b01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'L6b'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$lamp501, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Lamp5'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$micro01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Micro PVM'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$oligo01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Oligo'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$opc01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'OPC'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$pvalb01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Pvalb'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$sncg01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sncg'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$sst01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Sst'
stats3 <- rbind(stats3, df)
# df <- chisq.test(dat3$sstchodl01, dat3$diag01)
# df <- as.data.frame(df$p.value)
# df$celltype <- 'Sst Chodl'
# stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$vip01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'Vip'
stats3 <- rbind(stats3, df)
df <- chisq.test(dat3$vlmc01, dat3$diag01)
df <- as.data.frame(df$p.value)
df$celltype <- 'VLMC'
stats3 <- rbind(stats3, df)

stats$comparison <- 'cont_early'
stats2$comparison <- 'cont_late'
stats3$comparison <- 'early_late'
chi2stats <- rbind(stats, stats2)
chi2stats <- rbind(chi2stats, stats3)
names(chi2stats)[names(chi2stats) == "df$p.value"] <- "pvalue"
head(chi2stats)
chi2stats$fdr <- p.adjust(chi2stats$pvalue, method="fdr")
write.csv(chi2stats, file="~/celltype_mapping/chi2proportion_stats.csv")




#save the metadata file with predictions, scores, and misclass flag
write.csv(mathys_meta3, file="~/celltype_mapping/data/mathys_new_celltypes_seurat4.csv")

#EXPORT CDS SYNAPSE FOLDER TO IDENTIFY INDIVIDUAL BRANCHES FOR PSEUDOTIME ANALYSIS
#saveRDS(mathys_meta3, file="~/celltype_mapping/data/mathys_new_celltypes_seurat4.csv")
file <- synapser::File(path='~/celltype_mapping/data/mathys_new_celltypes_seurat4.csv', parentId='syn25871777')
file <- synapser::synStore(file)
