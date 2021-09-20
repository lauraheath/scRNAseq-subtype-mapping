#biodomains

#install old version of RcppAnnoy in order to successfully install later packages from bioconductor:
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/RcppAnnoy/RcppAnnoy_0.0.14.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
#setup
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("fgsea")
devtools::install_github('https://github.com/YuLab-SMU/yulab.utils')
BiocManager::install( "clusterProfiler")

library(org.Hs.eg.db)
library(fgsea)
library(tidyverse)
install.packages("yulab.utils")
BiocManager::install( "yulab.utils")
# data -- biological domain definitions
biodomObj <- synapser::synGet('syn25764758')
biodom <- readRDS(biodomObj$path)

domains <- biodom %>% pull(Biodomain) %>% unique() %>% sort()
sub.dom <- list()
for(i in 1:length(domains)){ 
  sub.dom[[i]] <- biodom %>% filter(Biodomain == domains[i], !is.na(n_ensGene))
}
# Colors for plotting domains
dom.cols <- c(
  'APP Metabolism' = '#ff6600',
  'Autophagy' = '#9933ff',
  'Endolysosomal' = '#3366cc',
  'Epigenetic Regulation' = '#cc3333',
  'Epigenetic' = '#cc3333',
  'Immune Response' = '#99cccc',
  'Lipid Metabolism' = '#999999', 
  'Mitochondria Metabolism'= '#99cc99', 
  'Mitochondrial Metabolism'= '#99cc99', 
  'Myelination' = '#996633',
  'Oxidative Stress' = '#ffcc66',
  'Apoptosis Regulation'= '#663399', 
  'Regulation of Apoptosis'= '#663399', 
  'RNA Spliceosome' = '#0099ff',
  'Structural Stabilization'= '#ff9999', 
  'Synaptic Dysfunction' = '#339933',
  'Tau Homeostasis' = '#cc99cc',
  'Vascular Function'= '#cccc00'
)



# data -- pseudotime ANOVA results tables
#dlpfc.pt <- read_csv(file="~/prot-lineage/results/female_DEanova_stats.csv") 
#DEobj <- synapser::synGet('syn25764747')
dlpfc.pt <- read.csv(file="~/scRNAseq-subtype-mapping/female_astro_DEanova.csv")
#run over again from here for males
#dlpfc.pt <- read_csv(file="~/prot-lineage/results/male_DEanova_stats.csv") 
#DEobj <- synapser::synGet('syn25764751')
#dlpfc.pt <- read.csv(DEobj$path) 
#RNASeq results uploaded into data folder from Mukherjee et al supplementary material
# dlpfc.pt <- read.csv(file="~/prot-lineage/data/DLPFC_DEgenes_rnaseq.csv")
# names(dlpfc.pt)[names(dlpfc.pt) == "gene_names"] <- "gene_short_name"
# 
# dlpfc.pt <- read.csv(file="~/prot-lineage/data/TCX_DEgenes_rnaseq.csv")
names(dlpfc.pt)[names(dlpfc.pt) == "gene_names"] <- "gene_short_name"
# ID re-mapping -- biodomains are defined by ENSEMBL and/or ENTREZ IDs
# (a small percentage will fail to map)
ids <- clusterProfiler::bitr( 
  dlpfc.pt$gene_short_name, fromType = 'SYMBOL', toType = c('ENTREZID','ENSEMBL'), OrgDb='org.Hs.eg.db' )
dlpfc.pt <- ids %>% 
  left_join(dlpfc.pt, ., by = c('gene_short_name'='SYMBOL')) %>% 
  relocate('ENTREZID', 'ENSEMBL', .after='gene_short_name')

#effect distributions by state
dlpfc.pt %>% 
  ggplot(., aes(effect))+
  geom_histogram(binwidth =.1)+
  facet_wrap(~state)+
  theme_minimal()

#pval distributions by state
dlpfc.pt %>% 
  ggplot(., aes(-log10(pvalue)))+
  geom_histogram(binwidth = .5)+
  facet_wrap(~state)+
  theme_minimal()



#run gsea and generate plots

enr_plots <- list()
lol_plots <- list()
biodom.annotated <- biodom %>% filter(!is.na(n_ensGene)) %>% pull(ensembl_id, name=GOterm_Name)
for(s in unique(dlpfc.pt$state) ){
  # set up genelist
  gl <- dlpfc.pt %>% 
    filter( state==s, pvalue < 0.05 ) %>% 
    filter( !is.na(ENSEMBL), !duplicated(ENSEMBL), !is.na(effect) ) %>% 
    arrange( desc(effect) ) %>% pull( effect, name=ENSEMBL )
  
  # run GSEA analysis
  bd.fgsea <-  fgseaMultilevel( pathways   = biodom.annotated,
                                stats      = gl,
                                minSize    = 1,
                                maxSize    = Inf,
                                #scoreType  = 'std',
                                #eps        = 0,
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





for(s in unique(dlpfc.pt$state)) {
  cat("  \n### state:",  s, "  \n")
  try(gridExtra::grid.arrange(enr_plots[[s]]), silent=T)
  plot(lol_plots[[s]])
  plot.new()
  dev.off()
  #knit_expand(text = paste0("\n## ", {{names(combined.plots[i])}},  "\n  ")
  cat("  \n  \n")  
}



#output individual ranked GO term plot for state 2
try(gridExtra::grid.arrange(enr_plots[[2]]), silent=T)
try(gridExtra::grid.arrange(enr_plots[[3]]), silent=T)
try(gridExtra::grid.arrange(enr_plots[[4]]), silent=T)

#plot biodomains dotplot for state 2 only
plot(lol_plots[[2]])
plot(lol_plots[[3]])
plot(lol_plots[[4]])


