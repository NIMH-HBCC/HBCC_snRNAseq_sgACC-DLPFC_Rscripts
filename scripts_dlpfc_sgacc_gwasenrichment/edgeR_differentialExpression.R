# edgeR pseudobulk DEX between brain regions 
# implemented using code from 'Libra' R package 
# corrected for donor status 


library(Libra) # tool kit ... 
library(limma)

dat <- readRDS('seuobj.RDS')

pseudobulks = to_pseudobulk( # from Libra R package 
  input = dat@assays$RNA@counts[,],
  meta = dat@meta.data, 
  min_features = 0
)


# run DEX between regions using edgeR LRT 
# libra code modified to include donor as a covariate 
results = map(pseudobulks, function(x) {
  # create targets matrix
  targets = data.frame(group_sample = colnames(x)) %>%
    mutate(group = gsub(".*\\:", "", group_sample),
           donor = gsub("\\:.*", "", group_sample)
    )
  ## optionally, carry over factor levels from entire dataset
  if (is.factor(meta$label)) {
    targets$group %<>% factor(levels = levels(meta$label))
  }
  if (n_distinct(targets$group) > 2)
    return(NULL)
  # create design
  design = model.matrix(~ donor + group, data = targets)
  
  DE = switch(de_method,
              edgeR = {
                tryCatch({
                  y = DGEList(counts = x, group = targets$group) %>%
                    calcNormFactors(method = 'TMM') %>%
                    estimateDisp(design)
                  test = switch(de_type,
                                QLF = {
                                  fit = glmQLFit(y, design)
                                  test = glmQLFTest(fit, coef = 'groupsgACC')
                                },
                                LRT = {
                                  fit = glmFit(y, design = design)
                                  test = glmLRT(fit)
                                })
                  res = topTags(test, n = Inf) %>%
                    as.data.frame() %>%
                    rownames_to_column('gene') %>%
                    # flag metrics in results
                    mutate(de_family = 'pseudobulk',
                           de_method = de_method,
                           de_type = de_type)
                }, error = function(e) {
                  message(e)
                  data.frame()
                })
              }
  )
}) 
results %<>% bind_rows(.id = 'cell_type')
DE = results 

# clean up the output
suppressWarnings(
  colnames(DE) %<>%
    fct_recode('p_val' = 'p.value',  ## DESeq2
               'p_val' = 'pvalue',  ## DESeq2
               'p_val' = 'p.value',  ## t/wilcox
               'p_val' = 'P.Value',  ## limma
               'p_val' = 'PValue'  , ## edgeR
               'p_val_adj' = 'padj', ## DESeq2/t/wilcox
               'p_val_adj' = 'adj.P.Val',      ## limma
               'p_val_adj' = 'FDR',            ## edgeER
               'avg_logFC' = 'log2FoldChange', ## DESEeq2
               'avg_logFC' = 'logFC', ## limma/edgeR
               'avg_logFC' = 'avg_log2FC' # Seurat V4
    )
) %>%
  as.character()

DE %<>%
  # calculate adjusted p values
  group_by(cell_type) %>%
  mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
  # make sure gene is a character not a factor
  mutate(gene = as.character(gene)) %>%
  # invert logFC to match Seurat level coding
  mutate(avg_logFC = avg_logFC * -1) %>%
  dplyr::select(cell_type,
                gene,
                avg_logFC,
                p_val,
                p_val_adj,
                de_family,
                de_method,
                de_type
  ) %>%
  ungroup() %>%
  arrange(cell_type, gene)
head(DE)

DE_edger = DE 

# END 

