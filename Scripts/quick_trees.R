library(tidyverse)
library(SNPRelate)
library(ggtree)
options(scipen = 999)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

# --------------------------------------------------------------------------------------------------------- # GLOBAL
gds_location <- "~/Dropbox/andersen_lab/divergent_regions/divergent_haplotype_analysis/Ce328_clean.gds"

# load up variant data in gds format
(genofile <- snpgdsOpen(gds_location))

# all snp information
snp_info <- data.frame(chrom = read.gdsn(index.gdsn(genofile, "snp.chromosome")),
                       pos = read.gdsn(index.gdsn(genofile, "snp.position")),
                       snpid = read.gdsn(index.gdsn(genofile, "snp.id")),
                       snp_allele = read.gdsn(index.gdsn(genofile, "snp.allele")))


assembled_strains <- c("N2","CB4856","DL238","ECA36","ECA396","EG4725","JU1400","JU2526","JU2600","JU310","MY2147","MY2693","NIC2","NIC526","QX1794", "XZ1516")

# --------------------------------------------------------------------------------------------------------- # REGION SPECIIFIC 

test_subset <- c("V:20193463-20267244")

quick_tree <- function(reg_oi = test_subset, 
                       strain_set = assembled_strains,
                       save_dir = "tree_outputs",
                       maf_cut = 0.006){
  
  # region info
  rchrom <- strsplit(reg_oi, split = ":")[[1]][1]
  swind <- as.numeric(strsplit(strsplit(reg_oi, split = ":")[[1]][2], split = "-")[[1]][1])
  ewind <- as.numeric(strsplit(strsplit(reg_oi, split = ":")[[1]][2], split = "-")[[1]][2])
  
  # get snp ids
  region_snppos <- snp_info %>%
    dplyr::filter(chrom == rchrom, pos > swind, pos < ewind)
  
  region_snpid <- region_snppos %>%
    dplyr::pull(snpid)
  
  number_of_snps <- length(region_snpid)
  
  # identity by state clustering
  hc_f <- snpgdsHCluster(snpgdsIBS(genofile, 
                                   num.thread=2, 
                                   autosome.only= F, 
                                   maf = maf_cut,
                                   sample.id = strain_set, 
                                   snp.id = region_snpid))
  
  rv_f <- snpgdsCutTree(hc_f, z.threshold = 30)
  
  # dendro to phylo object
  region_phylo_f <- ape::as.phylo(as.hclust(rv_f$dendrogram))
  
  # plot and save
  reg_tree_f <- ggtree::ggtree(region_phylo_f) +
    geom_tiplab(size=2) + 
    theme_tree2()
  
  ggsave(reg_tree_f, 
         filename = glue::glue("{rchrom}_{swind}-{ewind}_{number_of_snps}.pdf"), 
         height = 4, 
         width = 8)
  
  ape::write.tree(region_phylo_f, 
                  file = glue::glue("{rchrom}_{swind}-{ewind}_{number_of_snps}.tree"), 
                  append = FALSE,
                  digits = 10, 
                  tree.names = FALSE)
  
  return(list(region_phylo_f, reg_tree_f))
}

quick_tree()

# close file
closefn.gds(genofile)
