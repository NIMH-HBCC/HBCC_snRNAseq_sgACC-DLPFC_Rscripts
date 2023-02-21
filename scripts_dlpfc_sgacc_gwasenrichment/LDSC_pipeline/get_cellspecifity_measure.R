# script to generate cell specificity 
# adapted from https://github.com/jbryois/scRNA_disease 

library(Seurat)
library(tidyverse)
dat <- readRDS('seuobj.RDS')

# filter for acc 
table(dat$br)
dat <- subset(dat, subset = br == 'DLPFC')
table(dat$br)

dat
#dat@assays$RNA@counts
dim(dat@assays$RNA@counts)

exp <- as.data.frame((dat@assays$RNA@counts[,]))
exp$Gene <- rownames(exp)

# Add Cell type names
cell_types <- as.data.frame(dat@meta.data)
cell_types$column <- rownames(cell_types)

cell_types$Lvl5 <- Idents(dat)
cell_types$Lvl4 <- cell_types$broad.class
cell_types$Lvl5 <- paste0(cell_types$br, cell_types$Lvl5)


toclust <- unique(cell_types$Lvl5); toclust
sum.exp <- exp[,1:(length(toclust)+1)]
for (i in 1:length(toclust)) { 
  cellidforcelltype <- cell_types$column[which(cell_types$Lvl5 == toclust[[i]])]
  cat(length(cellidforcelltype), "\n")
  sum.exp[,i+1] <- rowSums(exp[,which(names(exp) %in% cellidforcelltype)])
  names(sum.exp)[i+1] <- paste("V",i,sep = "")
}
sum.exp[,1] <- exp$Gene; names(sum.exp)[1] <- 'Gene'
exp <- sum.exp

# filtering duplicate gene names
exp <- exp %>% add_count(Gene) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as.tibble()

gene_coordinates <- 
  read_tsv("../NCBI37.3.gene.loc.extendedMHCexcluded.txt",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-100000<0,0,X3-100000),end=X4+100000) %>%
  select(X2,start,end,1) %>% 
  rename(chr="X2", ENTREZ="X1") %>% 
  mutate(chr=paste0("chr",chr))


cell_types <- cbind(column=as.character(paste0("V",1:length(toclust))),
                    Lvl4=as.character(toclust),
                    Description=as.character(toclust)) %>%  
  as.tibble()
cell_types

exp_lvl5 <- inner_join(exp,cell_types,by="column") %>% ungroup() %>% rename(Expr_sum_mean=Expr)
head(exp_lvl5)

# Remove not expressed genes
not_expressed <- exp_lvl5 %>% 
  group_by(Gene) %>% 
  summarise(total_sum=sum(Expr_sum_mean)) %>% 
  filter(total_sum==0) %>% 
  select(Gene) %>% unique() 

exp_lvl5 <- filter(exp_lvl5,!Gene%in%not_expressed$Gene)

# Scale data 
exp_lvl5 <- exp_lvl5 %>% 
  group_by(Lvl4) %>% 
  mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))

# specificity cal
# lvl5 
exp_lvl5 <- exp_lvl5 %>% 
  group_by(Gene) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>% 
  ungroup()

# magma
# Keep only ENTREZ genes protein coding and MAGMA tested genes
entrez2symbol <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL2EG) %>% rename(Gene="symbol",ENTREZ="gene_id")

exp_lvl5 <- inner_join(exp_lvl5,entrez2symbol,by="Gene") 
exp_lvl5 <- inner_join(exp_lvl5,gene_coordinates,by="ENTREZ")

# no of genes 
n_genes <- length(unique(exp_lvl5$ENTREZ))
n_genes_to_keep <- (n_genes * 0.1) %>% round()

# save expression profiling for other processing 
save(exp_lvl5,file = "expression.ready.Rdata")

# functions to save 
magma_top10 <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group_magma(.,Cell_type))
}

write_group_magma  = function(df,Cell_type) {
  df <- select(df,Lvl4,ENTREZ)
  df_name <- make.names(unique(df[1]))
  colnames(df)[2] <- df_name  
  dir.create(paste0("MAGMA/"), showWarnings = FALSE)
  select(df,2) %>% t() %>% as.data.frame() %>% rownames_to_column("Cat") %>%
    write_tsv("MAGMA/top10.txt",append=T)
  return(df)
}

# get LDSC input 10% 
write_group  = function(df,Cell_type) {
  df <- select(df,Lvl4,chr,start,end,ENTREZ)
  dir.create(paste0("LDSC/Bed"), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0("LDSC/Bed/",make.names(unique(df[1])),".bed"),col_names = F)
  return(df)
}

ldsc_bedfile <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group(.,Cell_type))
}

# write MAGMA/LDSC input files 
exp_lvl5 %>% filter(Expr_sum_mean>1) %>% magma_top10("Lvl4")

## Warning: group_by_() is deprecated. 
## Please use group_by() instead
## 
## The 'programming' vignette or the tidyeval book can help you
## to program with group_by() : https://tidyeval.tidyverse.org
## This warning is displayed once per session.

exp_lvl5 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("Lvl4")





