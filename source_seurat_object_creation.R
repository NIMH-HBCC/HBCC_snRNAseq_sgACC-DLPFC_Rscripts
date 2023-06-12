# 11/02/2018
# creates source seurat object for processing and analysis
# 



# read in sample IDs
sample_IDs<-read.table("/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/sample_IDs", sep="\t")
# create Seurat objects for each sample and label cells with sample ID
# keep all genes and all cells
for (sample in seq((length(sample_IDs$V1)))){
  sample<-as.character(sample_IDs$V1[sample])
  print(paste("creating Seurat object for sample:",sample))
  # read in the dataset from 10X Genomics output from path
  temp_path<-"/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/snRNA_09132018/copy_to_local/filtered_matrices/GRCh38-1.2.0_premrna/"
  path2dir<-gsub("filtered_matrices", paste(sample,"_filtered_matrices",sep=""),temp_path)
  data<-Read10X(data.dir=path2dir)
  # Add sample_ID to each cell name (column name) because
  # sometimes cell barcodes are not unique across many library prepartaions
  # Also, this way, you can always know which sample the cell comes from.
  rep=substr(sample, nchar(sample), nchar(sample))
  colnames(data)<-paste(sample,"_",(colnames(data)),sep="")
  # Create Seurat object for each sample
  seur_obj<-CreateSeuratObject(raw.data=data, min.cells=0, min.genes=0)
  assign(paste("s",sample,sep="_"), seur_obj)
  rm(seur_obj,data)
}

## Merging Seurat objects
# for each subject, merge replicate objects,
# add sample_ID, replicate number, and batch number to metadata table
# 

# BR1443
s_1443_sgACC<-MergeSeurat(s_1443_sgACC_1,s_1443_sgACC_2)
s_1443_sgACC@meta.data$sample_ID<-(substr(rownames(s_1443_sgACC@meta.data),1,12))
s_1443_sgACC@meta.data$rep<-(substr(rownames(s_1443_sgACC@meta.data),12,12))
rm(s_1443_sgACC_1,s_1443_sgACC_2)
s_1443_DLPFC<-MergeSeurat(s_1443_DLPFC_1,s_1443_DLPFC_2)
s_1443_DLPFC@meta.data$sample_ID<-(substr(rownames(s_1443_DLPFC@meta.data),1,12))
s_1443_DLPFC@meta.data$rep<-(substr(rownames(s_1443_DLPFC@meta.data),12,12))
rm(s_1443_DLPFC_1, s_1443_DLPFC_2)
s_1443_sgACC@meta.data$batch<-"1"
s_1443_DLPFC@meta.data$batch<-"1"

# BR2543
s_2543_sgACC<-MergeSeurat(s_2543_sgACC_1,s_2543_sgACC_2)
s_2543_sgACC@meta.data$sample_ID<-(substr(rownames(s_2543_sgACC@meta.data),1,12))
s_2543_sgACC@meta.data$rep<-(substr(rownames(s_2543_sgACC@meta.data),12,12))
rm(s_2543_sgACC_1,s_2543_sgACC_2)
s_2543_DLPFC<-MergeSeurat(s_2543_DLPFC_1,s_2543_DLPFC_2)
s_2543_DLPFC@meta.data$sample_ID<-(substr(rownames(s_2543_DLPFC@meta.data),1,12))
s_2543_DLPFC@meta.data$rep<-(substr(rownames(s_2543_DLPFC@meta.data),12,12))
rm(s_2543_DLPFC_1, s_2543_DLPFC_2)
s_2543_sgACC@meta.data$batch<-"1"
s_2543_DLPFC@meta.data$batch<-"1"

# BR1281
s_1281_sgACC<-MergeSeurat(s_1281_sgACC_1,s_1281_sgACC_2)
s_1281_sgACC@meta.data$sample_ID<-(substr(rownames(s_1281_sgACC@meta.data),1,12))
s_1281_sgACC@meta.data$rep<-(substr(rownames(s_1281_sgACC@meta.data),12,12))
rm(s_1281_sgACC_1,s_1281_sgACC_2)
s_1281_DLPFC<-MergeSeurat(s_1281_DLPFC_1,s_1281_DLPFC_2)
s_1281_DLPFC@meta.data$sample_ID<-(substr(rownames(s_1281_DLPFC@meta.data),1,12))
s_1281_DLPFC@meta.data$rep<-(substr(rownames(s_1281_DLPFC@meta.data),12,12))
rm(s_1281_DLPFC_1, s_1281_DLPFC_2)
s_1281_sgACC@meta.data$batch<-"2"
s_1281_DLPFC@meta.data$batch<-"2"

# BR2521
s_2521_sgACC<-MergeSeurat(s_2521_sgACC_1,s_2521_sgACC_2)
s_2521_sgACC@meta.data$sample_ID<-(substr(rownames(s_2521_sgACC@meta.data),1,12))
s_2521_sgACC@meta.data$rep<-(substr(rownames(s_2521_sgACC@meta.data),12,12))
rm(s_2521_sgACC_1,s_2521_sgACC_2)
s_2521_DLPFC<-MergeSeurat(s_2521_DLPFC_1,s_2521_DLPFC_2)
s_2521_DLPFC@meta.data$sample_ID<-(substr(rownames(s_2521_DLPFC@meta.data),1,12))
s_2521_DLPFC@meta.data$rep<-(substr(rownames(s_2521_DLPFC@meta.data),12,12))
rm(s_2521_DLPFC_1, s_2521_DLPFC_2)
s_2521_sgACC@meta.data$batch<-"2"
s_2521_DLPFC@meta.data$batch<-"2"

# BR2587
s_2587_sgACC<-MergeSeurat(s_2587_sgACC_1,s_2587_sgACC_2)
s_2587_sgACC@meta.data$sample_ID<-(substr(rownames(s_2587_sgACC@meta.data),1,12))
s_2587_sgACC@meta.data$rep<-(substr(rownames(s_2587_sgACC@meta.data),12,12))
rm(s_2587_sgACC_1,s_2587_sgACC_2)
s_2587_DLPFC<-MergeSeurat(s_2587_DLPFC_1,s_2587_DLPFC_2)
s_2587_DLPFC@meta.data$sample_ID<-(substr(rownames(s_2587_DLPFC@meta.data),1,12))
s_2587_DLPFC@meta.data$rep<-(substr(rownames(s_2587_DLPFC@meta.data),12,12))
rm(s_2587_DLPFC_1, s_2587_DLPFC_2)
s_2587_sgACC@meta.data$batch<-"3"
s_2587_DLPFC@meta.data$batch<-"3"

# BR1735
s_1735_sgACC<-MergeSeurat(s_1735_sgACC_1,s_1735_sgACC_2)
s_1735_sgACC@meta.data$sample_ID<-(substr(rownames(s_1735_sgACC@meta.data),1,12))
s_1735_sgACC@meta.data$rep<-(substr(rownames(s_1735_sgACC@meta.data),12,12))
rm(s_1735_sgACC_1,s_1735_sgACC_2)
s_1735_DLPFC<-MergeSeurat(s_1735_DLPFC_1,s_1735_DLPFC_2)
s_1735_DLPFC@meta.data$sample_ID<-(substr(rownames(s_1735_DLPFC@meta.data),1,12))
s_1735_DLPFC@meta.data$rep<-(substr(rownames(s_1735_DLPFC@meta.data),12,12))
rm(s_1735_DLPFC_1, s_1735_DLPFC_2)
s_1735_sgACC@meta.data$batch<-"3"
s_1735_DLPFC@meta.data$batch<-"3"

# BR2285
s_2285_sgACC<-MergeSeurat(s_2285_sgACC_1,s_2285_sgACC_2)
s_2285_sgACC@meta.data$sample_ID<-(substr(rownames(s_2285_sgACC@meta.data),1,12))
s_2285_sgACC@meta.data$rep<-(substr(rownames(s_2285_sgACC@meta.data),12,12))
rm(s_2285_sgACC_1,s_2285_sgACC_2)
s_2285_DLPFC<-MergeSeurat(s_2285_DLPFC_1,s_2285_DLPFC_2)
s_2285_DLPFC@meta.data$sample_ID<-(substr(rownames(s_2285_DLPFC@meta.data),1,12))
s_2285_DLPFC@meta.data$rep<-(substr(rownames(s_2285_DLPFC@meta.data),12,12))
rm(s_2285_DLPFC_1, s_2285_DLPFC_2)
s_2285_sgACC@meta.data$batch<-"4"
s_2285_DLPFC@meta.data$batch<-"4"

# BR2582
s_2582_sgACC<-MergeSeurat(s_2582_sgACC_1,s_2582_sgACC_2)
s_2582_sgACC@meta.data$sample_ID<-(substr(rownames(s_2582_sgACC@meta.data),1,12))
s_2582_sgACC@meta.data$rep<-(substr(rownames(s_2582_sgACC@meta.data),12,12))
rm(s_2582_sgACC_1,s_2582_sgACC_2)
s_2582_DLPFC<-MergeSeurat(s_2582_DLPFC_1,s_2582_DLPFC_2)
s_2582_DLPFC@meta.data$sample_ID<-(substr(rownames(s_2582_DLPFC@meta.data),1,12))
s_2582_DLPFC@meta.data$rep<-(substr(rownames(s_2582_DLPFC@meta.data),12,12))
rm(s_2582_DLPFC_1, s_2582_DLPFC_2)
s_2582_sgACC@meta.data$batch<-"4"
s_2582_DLPFC@meta.data$batch<-"4"


## Merge Seurat objects by brain region
# add brain region to meta data

# sgACC
sgACC1<-MergeSeurat(s_1443_sgACC,s_2543_sgACC)
sgACC2<-MergeSeurat(s_1281_sgACC,s_2521_sgACC)
sgACC3<-MergeSeurat(s_2587_sgACC,s_1735_sgACC)
sgACC4<-MergeSeurat(s_2285_sgACC,s_2582_sgACC)
sgACC5<-MergeSeurat(sgACC1,sgACC2)
rm(sgACC1,sgACC2)
sgACC6<-MergeSeurat(sgACC3,sgACC4)
rm(sgACC3,sgACC4)
sgACC<-MergeSeurat(sgACC5,sgACC6)
rm(sgACC5,sgACC6)
sgACC@meta.data$br<-"sgACC"

# DLPFC
DLPFC1<-MergeSeurat(s_1443_DLPFC,s_2543_DLPFC)
DLPFC2<-MergeSeurat(s_1281_DLPFC,s_2521_DLPFC)
DLPFC3<-MergeSeurat(s_2587_DLPFC,s_1735_DLPFC)
DLPFC4<-MergeSeurat(s_2285_DLPFC,s_2582_DLPFC)
DLPFC5<-MergeSeurat(DLPFC1,DLPFC2)
rm(DLPFC1,DLPFC2)
DLPFC6<-MergeSeurat(DLPFC3,DLPFC4)
rm(DLPFC3,DLPFC4)
DLPFC<-MergeSeurat(DLPFC5,DLPFC6)
rm(DLPFC5,DLPFC6)
DLPFC@meta.data$br<-"DLPFC"


## Merge both brain regions together
# all
all<-MergeSeurat(sgACC,DLPFC)

################################## Calculate Meta Data #################################
# frac mito
exp.mat<-all@raw.data
mgenes<-grep(pattern="^MT-", x=rownames(exp.mat), value=T)
all@meta.data$frac.mito<-colSums(all@raw.data[mgenes,])/colSums(all@raw.data)
# frac ribo
rpgenes<-grep(pattern="^MRPL|^MRPS|^RPS|^RPL",x=rownames(exp.mat),value=T)
all@meta.data$frac.ribo<-colSums(all@raw.data[rpgenes,])/colSums(all@raw.data)

################################## save file ##################################
saveRDS(all, file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/source_seur_obj.RDS")
