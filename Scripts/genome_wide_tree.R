library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(ape)
library(tibble)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))
load("Processed_Data/admix_K16.RData")

####### Fig. 1b - Genome-wide tree colored by unadmixed strains (grouped by admixture ancestry)

## tree generation

# Generate phylogeny files

#system(glue::glue("python Scripts/vcf2phylip-master/vcf2phylip.py -i Processed_Data/Ce328_complete_sites.vcf.gz -m 312 --fasta --nexus --nexus-binary"))
#WI_phy <- read.phyDat("Processed_Data/Ce328_complete_sites.min312.phy", format = "interleaved")
#save(WI_phy, file = "Processed_Data/WI_phy.RData")

load("Processed_Data/WI_phy.RData")
load("Processed_Data/strain_geo.RData")

dm <- dist.ml(WI_phy)
treeNJ <- NJ(dm)
tree_pt<-ggtree(treeNJ,
                branch.length="rate")
plot_tree <- ggtree(treeNJ, layout="fan", open.angle=60, branch.length="rate")

## genome-wide tree with labels for geographic origins

isotype_info_geo <- indep_strain_info_geo %>%
  dplyr::select(strain, geo) %>%
  dplyr::select(label=strain, geo)

df_node <- plot_tree[[1]] %>%
  dplyr::select(node, label)

df_node_geo <- df_node %>% 
  dplyr::left_join(., isotype_info_geo, by = 'label') %>%
  dplyr::distinct(node, geo) 

df_node_geo$node <- as.numeric(df_node_geo$node)

tbl_tree_geo <- as.tibble(treeNJ) %>%
  dplyr::full_join(., df_node_geo, by = 'node') %>%
  dplyr::rename(label = geo, strain = label) %>%
  dplyr::mutate(label = ifelse(is.na(label), "Z_link", as.character(label))) %>%
  dplyr::mutate(branch.length = ifelse(branch.length < 0, 0, branch.length)) ## to control negative branch length

tree_ph_geo <- as.phylo(tbl_tree_geo, use.labels = T)

#save(tree_ph_geo, file = "Processed_Data/tree_genome_wide.RData")
load("Processed_Data/tree_genome_wide.RData")

plot_tree_NJ_geo <- ggtree(tree_ph_geo, layout="fan", open.angle=60, branch.length="rate") + 
  aes(color = label, size = label) + 
  theme(legend.position = c(0.68,0.4), 
        legend.title = element_blank(), 
        legend.text = element_text(size=9, color = "black"),
        legend.direction = "horizontal", 
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.1,"in"),
        plot.margin = margin(b=-7, l=-8, r=-1.6, t=-8, unit = "in")) +
  #theme(legend.position = 'none', legend.text = element_text(size = 10), legend.title = element_blank(), plot.margin = margin(b=-8, l=-10, r=-6, t=-10, unit = "in")) +
  scale_color_manual(values = c(geo.colours, "Z_link"="gray")) + 
  guides(col= guide_legend(ncol=9))  +
  scale_size_manual(values=c(.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.2))

plot_tree_NJ_geo_rotate <- rotate_tree(plot_tree_NJ_geo,90) %>% flip(639, 642) %>% flip(579,585)

plot_tree_NJ_geo_rotate

ggsave("tree_geo.pdf", plot_tree_NJ_geo_rotate, width = 7.5, height = 5)

#ggsave("Plots/plot_tree_NJ_admix_rotate.pdf", plot_tree_NJ_admix_rotate, width = 30, height = 10)