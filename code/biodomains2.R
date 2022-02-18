#biodomains

BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)

#setup
BiocManager::install("fgsea")
BiocManager::install("org.Hs.eg.db")
devtools::install_github('https://github.com/YuLab-SMU/yulab.utils')
BiocManager::install( "clusterProfiler")

# setup: packages, data ---------------------------------------------------

require(fgsea)
require(synapser)
library(tidyverse)

synapser::synLogin()








# data -- biological domain definitions

biodom <- readRDS( synapser::synGet('syn25428992')$path )
dom.cols <- c(
  'APP Metabolism' = '#ff6600',
  'Autophagy' = '#9933ff',
  'Endolysosome' = '#3366cc',
  'Epigenetic' = '#cc3333',
  'Immune Response' = '#99cccc',
  'Lipid Metabolism' = '#999999', 
  'Mitochondrial Metabolism'= '#99cc99', 
  'Myelination' = '#996633',
  'Oxidative Stress' = '#ffcc66',
  'Proteostasis' = '#C7B169',
  'Apoptosis'= '#663399', 
  'RNA Spliceosome' = '#0099ff',
  'Structural Stabilization'= '#ff9999', 
  'Synapse' = '#339933',
  'Tau Homeostasis' = '#cc99cc',
  'Vasculature'= '#cccc00'
)


domains <- biodom %>% pull(Biodomain) %>% unique() %>% sort()
sub.dom <- list()
for(i in 1:length(domains)){ 
  sub.dom[[i]] <- biodom %>% dplyr::filter(Biodomain == domains[i], !is.na(n_ensGene))
}

#dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/microgliaF_samplecorrected_DEresults.csv")
dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/microgliaM_samplecorrected_DEresults.csv")
#dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/astrocyteF_samplecorrected_DEresults.csv")
#dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/astrocyteM_samplecorrected_DEresults.csv")


#dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/microgliaF_samplecorrected_ReferenceResults.csv")
#dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/microgliaM_samplecorrected_ReferenceResults.csv")
#dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/astrocyteF_samplecorrected_ReferenceResults.csv")
#dat <- read.csv(file="~/scRNAseq-subtype-mapping/DEresults/astrocyteM_samplecorrected_ReferenceResults.csv")

dat$X<-NULL
names(dat)[names(dat) == "gene"] <- "gene_short_name"

library(EnsDb.Hsapiens.v79)
ids2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= dat$gene_short_name, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
ids2 <- ids2[- grep("LRG", ids2$GENEID),]
dat <- ids2 %>%
  left_join(dat, ., by = c('gene_short_name'='SYMBOL')) %>%
  relocate('GENEID', .after='gene_short_name')
names(dat)[names(dat) == "GENEID"] <- "ENSEMBL"

# ids <- clusterProfiler::bitr(
#   dat$gene_short_name, fromType = 'SYMBOL', toType = c('ENTREZID','ENSEMBL'), OrgDb='org.Hs.eg.db' )
# dat <- ids %>%
#   left_join(dat, ., by = c('gene_short_name'='SYMBOL')) %>%
#   relocate('ENTREZID', 'ENSEMBL', .after='gene_short_name')


detach("package:EnsDb.Hsapiens.v79", unload = TRUE)
detach("package:ensembldb", unload = TRUE)


#effect distributions by state
dat %>% 
  ggplot(., aes(logFC_state))+
  geom_histogram(binwidth =.1)+
  facet_wrap(~state)+
  theme_minimal()

#pval distributions by state
dat %>% 
  ggplot(., aes(-log10(adjustp)))+
  geom_histogram(binwidth = .5)+
  facet_wrap(~state)+
  theme_minimal()



#run gsea and generate plots

enr_plots <- list()
lol_plots <- list()
biodom.annotated <- biodom %>% dplyr::filter(!is.na(n_ensGene)) %>% pull(ensembl_id, name=GOterm_Name)
for(s in unique(dat$state) ){
  # set up genelist
  gl <- dat %>% 
    filter( state==s, adjustp < 0.05 ) %>% 
    filter( !is.na(ENSEMBL), !duplicated(ENSEMBL), !is.na(logFC_state) ) %>% 
    arrange( desc(logFC_state) ) %>% pull( logFC_state, name=ENSEMBL )
  
  # run GSEA analysis
  bd.fgsea <-  fgseaMultilevel( pathways   = biodom.annotated,
                                stats      = gl,
                                minSize    = 1,
                                maxSize    = Inf,
                                #scoreType  = 'std',
                                eps        = 0,
                                nproc      = 4 )
  
  gsea.pathways <- bd.fgsea %>% filter(padj <= 0.05) %>%
    arrange(desc(NES)) %>% distinct() %>% head(n = 10) %>% pull(pathway)
  gsea.pathways <- bd.fgsea %>% filter(padj <= 0.05) %>%
    arrange(desc(NES)) %>% distinct() %>% tail(n = 10) %>% pull(pathway) %>% c(gsea.pathways, .) %>% unique()
  
  enr_plots[[s]] <- plotGseaTable(biodom.annotated[gsea.pathways], gl, bd.fgsea, render = FALSE) 
  
  sig <- bd.fgsea %>% filter(padj <= 0.05)
  bd.tally <- data.frame(domain=rep(domains,2), 
                         direction = c(rep('up', length(domains)),rep('dn', length(domains))) , 
                         n_sig_term=NaN, n_term = NaN)
  for(i in 1:nrow(bd.tally)){
    sd <- biodom %>% filter(Biodomain == bd.tally$domain[i]) 
    if(bd.tally$direction[i]=='up'){
      bd.tally$n_sig_term[i] <- sig %>% filter(pathway %in% sd$GOterm_Name, NES > 0) %>% 
        dplyr::select(pathway) %>% distinct() %>% nrow()
    }
    if(bd.tally$direction[i]=='dn'){
      bd.tally$n_sig_term[i] <- sig %>% filter(pathway %in% sd$GOterm_Name, NES < 0) %>% 
        dplyr::select(pathway) %>% distinct() %>% nrow()
    }
    bd.tally$n_term[i] <- sd %>% dplyr::select(GOterm_Name) %>% distinct() %>% nrow()
  }
  bd.tally$prop <- bd.tally$n_sig_term / bd.tally$n_term
  
  x <- bd.tally %>% filter(direction == 'up')
  x$domain <- with(x, reorder(domain, dplyr::desc(n_sig_term))) 
  x <- x %>% arrange(domain)
  p1 <- ggplot(x, aes(x = domain, y=n_sig_term)) +
    geom_segment( aes(xend=domain, yend=0) ) +
    geom_point( size = 4, color = dom.cols[as.character(sort(x$domain, decreasing =F))] ) + 
    theme_minimal() + labs(x='Biological Domain',y='# of significantly enriched GO terms') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle('pos NES')
  
  x <- bd.tally %>% filter(direction == 'dn')
  x$domain <- with(x, reorder(domain, dplyr::desc(n_sig_term))) 
  x <- x %>% arrange(domain)
  p2 <- ggplot(x, aes(x = domain, y=n_sig_term)) +
    geom_segment( aes(xend=domain, yend=0) ) +
    geom_point( size = 4, color = dom.cols[as.character(sort(x$domain, decreasing =F))] ) + 
    theme_minimal() + labs(x='Biological Domain',y='# of significantly enriched GO terms') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle('neg NES')
  
  lol_plots[[s]] <- cowplot::plot_grid(p1,p2,nrow=1)
}



#how many states? doublecheck:
table(dat$state)

#output individual ranked GO term plot for state 2
try(gridExtra::grid.arrange(enr_plots[[2]]), silent=T)
try(gridExtra::grid.arrange(enr_plots[[3]]), silent=T)
try(gridExtra::grid.arrange(enr_plots[[4]]), silent=T)
try(gridExtra::grid.arrange(enr_plots[[5]]), silent=T)

#plot biodomains dotplot for state 2 only
plot(lol_plots[[2]])
plot(lol_plots[[3]])
plot(lol_plots[[4]])
plot(lol_plots[[5]])



