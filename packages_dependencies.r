
#UPDATED:
#dependencies to run in terminal window before installing packages:
#sudo apt-get update -y
#sudo apt-get install -y libudunits2-dev libgdal-dev

## also set up synapse login configuration in terminal:
#nano .synapseConfig
##in nano window type the following:
#[authentication]
#username = <username>
#password = <password>

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

#install old version of RcppAnnoy in order to successfully install later packages from bioconductor:
packageurl <- "https://cran.r-project.org/src/contrib/Archive/RcppAnnoy/RcppAnnoy_0.0.14.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

BiocManager::install(c("Biobase", "BiocGenerics","DelayedArray","DelayedMatrixStats","limma","S4Vectors","SingleCellExperiment","SummarizedExperiment", "batchelor", "Matrix.utils"))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
remotes::install_github("jlmelville/uwot")
remotes::install_github("mojaveazure/seurat-disk")
devtools::install_github('satijalab/seurat-data')
install.packages("caret")

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(monocle3)
library(limma)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(caret)
