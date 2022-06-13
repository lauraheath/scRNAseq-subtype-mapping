# This script runs GSEA enrichment using branch-specific pseudotime effect sizes to rank genes and plots the results


# setup -------------------------------------------------------------------


# libraries
library(synapser)
library(fgsea)
library(tidyverse)


# plotting theme
theme_set(theme_bw())



#This function quantifies the number of significantly enriched GO terms by AD biological domain
bd.tally <- function( enrVct, biodomDefTbl){
  bdt <- bind_cols(
    domain = unique(biodomDefTbl$Biodomain),
    n_term = map_dbl( unique(biodomDefTbl$Biodomain),
                      ~ biodomDefTbl %>% filter(Biodomain == .x) %>% 
                        dplyr::select(GOterm_Name) %>% distinct() %>% nrow()),
    n_sig_term = map_dbl( unique(biodomDefTbl$Biodomain),
                          ~ enrVct[ enrVct %in% 
                                      biodomDefTbl$GOterm_Name[
                                        biodomDefTbl$Biodomain == .x]] %>% 
                            length()) ) %>% 
    mutate(domain = fct_reorder(domain, n_sig_term, .desc = F)) %>% 
    arrange(domain)
  bdt$prop <- bdt$n_sig_term / bdt$n_term
  return(bdt)
}


# data: biodomain definitions, biodomain plotting colors, and braak state-gene table
biodom <- readRDS( synapser::synGet('syn25428992')$path  )
biodom.annotated <- biodom %>% filter(!is.na(n_symbol)) %>% pull(symbol, name=GOterm_Name)
dom.cols <- read_csv( synGet('syn26856828')$path )
#braak.pt <- read_csv( synGet('syn27260460')$path )

#degs <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/AstroF_pstimeDE_eachState_vs_allstates.csv")
degs <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/AstroM_pstimeDE_eachState_vs_allstates.csv")


#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/microgliaF_DEGS_StatesVsAll.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/microgliaM_DEGS_StatesVsAll.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/astrocyteM_DEGS_StatesVsAll.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/oligoF_DEGS_StatesVsAll.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/oligoM_DEGS_StatesVsAll.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/opcF_DEGS_StatesVsAll.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/opcM_DEGS_StatesVsAll.csv')



# Run enrichments by pseudotime state -------------------------------------
# enrichment for each state; all significant genes (pval < 0.05), arranged by effect size
enr <- map_dfr(
  unique(degs$state),
  ~ degs %>% 
    filter(state == .x
           , p_state2 < 0.05
           # , !(gene_names %in% (dlpfc.pt %>% filter(pvalue <= 0.05) %>% pull(gene_names))) # use this to consider genes unique to the neuropath
    ) %>% #
    arrange(desc(logFC_state2)) %>% 
    pull(logFC_state2, name = gene) %>%
    fgseaMultilevel(biodom.annotated, ., eps = 0, scoreType = 'std') %>% 
    mutate(state = .x)
)


# count up the number of significant terms per biodomain
bdt <- bd.tally(enr$pathway[enr$pval <= 0.05], biodom) %>% 
  mutate(domain = fct_reorder(domain, n_sig_term, .desc = T)) %>% 
  arrange(domain)


# add biodomain annotations and plotting colors to enrichment results
enr <- enr %>% 
  left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain), by='pathway') %>%
  full_join(., dom.cols, by = c('Biodomain'='domain')) %>%
  mutate( Biodomain = fct_relevel(Biodomain, as.character(bdt$domain) )) %>% 
  mutate( state = fct_reorder( as.factor(state), as.numeric(state), .desc = T ))



# plot! -------------------------------------------------------------------

#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/MicroF_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/MicroM_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/astroF_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/astroM_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoF_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoM_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcF_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcM_stateVsAll_biodomains.tiff',height=150,width=200,units='mm',res=300)
enr %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enr$pval) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enr, pval > 0.01),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enr, pval < 0.01),# & !(pathway %in% idx)),
    aes(color = color, size = -log10(pval) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enr, pval < 0.01 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enr, pval < 0.01 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Biodomain)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Astro male: each state vs all; state pval < 0.01; term pval < 0.01')

#dev.off()

#save as csv file for state characterization later, if desired
enr$leadingEdge <- as.character(enr$leadingEdge)

#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/astroF_statevsall_enrichment.csv', row.names=FALSE)
write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/astroM_statevsall_enrichment.csv', row.names=FALSE)


#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/microF_statevsall_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/microM_statevsall_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/astroM_statevsall_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoF_statevsall_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoM_statevsall_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcF_statevsall_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcM_statevsall_enrichment.csv', row.names=FALSE)







##### Each state vs reference state 1 #####################################################

#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/microgliaF_DEGS_StatesVsRef.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/microgliaM_DEGS_StatesVsRef.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/astrocyteF_DEGS_StatesVsRef.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/astrocyteM_DEGS_StatesVsRef.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/oligoF_DEGS_StatesVsRef.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/oligoM_DEGS_StatesVsRef.csv')
#degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/opcF_DEGS_StatesVsRef.csv')
degs <- read.csv(file='~/scRNAseq-subtype-mapping/DEresults/opcM_DEGS_StatesVsRef.csv')

#run this for StatesVsRef
enr <- map_dfr(
  unique(degs$state),
  ~ degs %>% 
    filter(state == .x
           , p_state < 0.05
           # , !(gene_names %in% (dlpfc.pt %>% filter(pvalue <= 0.05) %>% pull(gene_names))) # use this to consider genes unique to the neuropath
    ) %>% #
    arrange(desc(logFC_state)) %>% 
    pull(logFC_state, name = gene) %>%
    fgseaMultilevel(biodom.annotated, ., eps = 0, scoreType = 'std') %>% 
    mutate(state = .x)
)

# count up the number of significant terms per biodomain
bdt <- bd.tally(enr$pathway[enr$pval <= 0.05], biodom) %>% 
  mutate(domain = fct_reorder(domain, n_sig_term, .desc = T)) %>% 
  arrange(domain)


# add biodomain annotations and plotting colors to enrichment results
enr <- enr %>% 
  left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain), by='pathway') %>%
  full_join(., dom.cols, by = c('Biodomain'='domain')) %>%
  mutate( Biodomain = fct_relevel(Biodomain, as.character(bdt$domain) )) %>% 
  mutate( state = fct_reorder( as.factor(state), as.numeric(state), .desc = T ))



# plot! -------------------------------------------------------------------
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/MicroF_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/MicroM_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/astroF_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/astroM_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoF_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoM_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcF_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcM_stateVsRef_biodomains.tiff',height=150,width=200,units='mm',res=300)

enr %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enr$pval) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enr, pval > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enr, pval < 0.05),# & !(pathway %in% idx)),
    aes(color = color, size = -log10(pval) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enr, pval < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enr, pval < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Biodomain)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('OPC male each state vs state 1; state pval < 0.05; term pval < 0.05')

dev.off()



#save as csv file for state characterization later, if desired
enr$leadingEdge <- as.character(enr$leadingEdge)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/microF_stateVsRef_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/microM_stateVsRef_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/astroF_stateVsRef_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/astroM_stateVsRef_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoF_stateVsRef_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/oligoM_stateVsRef_enrichment.csv', row.names=FALSE)
#write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcF_stateVsRef_enrichment.csv', row.names=FALSE)
write.csv(enr, file='~/scRNAseq-subtype-mapping/DEresults/To_export/opcM_stateVsRef_enrichment.csv', row.names=FALSE)
