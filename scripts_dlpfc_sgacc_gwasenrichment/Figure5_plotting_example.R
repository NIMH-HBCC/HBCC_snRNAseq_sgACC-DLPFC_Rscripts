# Figure 5 - plotting 
# requires LDSC GWAS enrichment, and Hypergeometric DEX gene overlap enrichment 
# a) LDSC enrichment calculated previously as Skene et al 2018. 
# b) Hypergeometric enrichment between top 10% specific genes and Bulk RNA DEX genes 
# function used to calculate presented below 
gethyperpval <- function(list1, list2 ) {
  # http://nemates.org/MA/progs/representation.stats.html 
  both <- length(unique(intersect(list1, list2))) # overlap between two lists
  g1 = length(unique(list1)) # size of g1
  g2 = length(unique(list2)) # size of g2
  
  pval <- phyper(both-1, g1, totalgenome-g1, g2, lower.tail = FALSE, log.p = FALSE)
  rf <- (both-1)/((g1*g2)/totalgenome) # representation factor (analgous to effect size)
  return(pval)
}

# results plotted via following function 
library(ggplot2)
library(ggnewscale)
library(ggh4x)

ggplot(dat, mapping = aes(Dx, Cell)) +
  geom_text( aes(label = mylab), size = 0, color = 'black', nudge_y = -0) + 
  geom_tile(color = "black", size=0.25, data = dat[which(dat$type == 'GWAS'),], aes(fill = -log10(p_adj2))) +
  scale_fill_distiller(palette = 'Reds', direction = 1,name = "FDR_P \n (-log10)", values = c(0.001,1), na.value = 'white') +
  new_scale_fill() + # used to have two different color scales 
  geom_tile(color = "black", size=0.25, data = dat[which( dat$type %like% 'DEG'),], aes(fill = -log10(p_adj2))) +
  geom_text(data = dat[which( dat$type %like% 'DEG'),], aes(label = mylab), size = 1.7, color = 'black', nudge_y = -0) + 
  geom_text(data = dat[which(dat$type == 'GWAS'),], aes(label = mylab), size = 1.7, color = 'black', nudge_y = -0) + 
  # Color scale applied to geoms added after new_scale_color()
  scale_fill_distiller(palette = 'Blues', direction = 1,name = "FDR_P \n (-log10)", values = c(0.001,1), na.value = 'white') +
  facet_nested(Subcat ~ newDX + type, scales='free', space='free') +
  theme(axis.text.y=element_text(angle = 0, size = 8, hjust = 1),
        axis.text.x=element_text(size = 4, angle = 90,hjust=1,vjust=0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 8),
        panel.spacing.x = unit(c(0), "lines"),
        aspect.ratio = 0.75) 

ggsave(last_plot(), file = 'Figure5_mixedscales_09072022.pdf', height = 10, width = 7)

#