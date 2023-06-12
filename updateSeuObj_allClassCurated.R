setwd("/data/HBCC_analysis/ajeet_HBCC_analysis/snRNASeq_dlpfc-sgacc_HBCCproject/")
setwd("/Volumes/HBCC_analysis/ajeet_HBCC_analysis/snRNASeq_project")

library(scran)

library(Seurat)
clustObject <- readRDS(file = "/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/source_seur_obj.RDS")

clustObject <- readRDS(file = "/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.RDS")

clustObject <- readRDS(file = "/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.scale.clust.tSNE.RDS")

clustObject <- readRDS(file = "/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.scale.clust.RDS")

clustObject <- readRDS(file = "/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_neurons_curated.RDS")

clustObject <- readRDS(file = "/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_classes_curated.RDS")



## check that Seurat package is at least v3.0
utils::packageVersion('Seurat') < 3 
## check version of Seurat object 
clustObject@version < 3

## UpdateSeuratObject: Update old Seurat object to accommodate new features
clustObject_updated <- UpdateSeuratObject(clustObject)
clustObject_updated

saveRDS(clustObject_updated, "updatedSeuObj_all_classes_curated.RDS")
###########
class(clustObject)
sce <- as.SingleCellExperiment(clustObject_updated)
sce@colData


