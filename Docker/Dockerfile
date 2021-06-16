FROM rocker/rstudio:4.1.0

RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev zlib1g-dev libbz2-dev liblzma-dev build-essential libglpk40\
&& apt-get install -y curl libcurl4-openssl-dev libxml2-dev\
&& apt-get install -y git\
&& R -e "install.packages('BiocManager')"\
&& R -e "BiocManager::install('biomaRt')"\
&& R -e "install.packages('caret')"\
&& R -e "install.packages('data.table')"\
&& R -e "install.packages('ggplot2')"\
&& R -e "install.packages('synapser', repos = c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"\
&& R -e "install.packages('devtools')"\
&& R -e "devtools::install_github('Sage-Bionetworks/knit2synapse')"\
&& R -e "install.packages('knitr')"\
&& R -e "BiocManager::install('limma')"\
&& R -e "install.packages('Matrix')"\
&& R -e "install.packages('patchwork')"\
&& R -e "install.packages('RColorBrewer')"\
&& R -e "install.packages('Seurat')"\
&& R -e "install.packages('tidyverse')"


run mkdir -p /home/lheath/celltype_mapping \ 
&& git clone https://github.com/lauraheath/scRNAseq-subtype-mapping.git /home/lheath/celltype_mapping
