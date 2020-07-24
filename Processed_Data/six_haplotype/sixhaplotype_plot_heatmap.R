#!/usr/bin/env Rscript

# set up
library(scales)
library(ggplot2)
library(reshape2)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# define heatmap plotting order based on newick strain tree
strain_list <- c("N2", "NIC2", "JU2600", "JU1400", "MY2147", "JU2526", "JU310", "EG4725", "MY2693", "ECA36", "DL238", "ECA396", "CB4856", "QX1794", "NIC526", "XZ1516")

# read in file
aa_dist <- read.table(args[1], row.names = 1)
aa_dist_m <- as.matrix(aa_dist)
colnames(aa_dist_m) <- rownames(aa_dist_m)

# reformat to long 
aa_dist_m_long <- melt(aa_dist_m)

# get a list of sequences
seq_list <- row.names(aa_dist_m)

# sort sequence list based on strain list
strain_list %>% 
  map(~.x %>% str_detect(unlist(seq_list), .) %>% which() %>% seq_list[[.]]) -> sorted_seq_list

# reverse it
sorted_seq_list <- rev(sorted_seq_list)

# order input file by sorted list
aa_dist_m_long$Var1 <- factor(aa_dist_m_long$Var1, levels = sorted_seq_list)
aa_dist_m_long$Var2 <- factor(aa_dist_m_long$Var2, levels = sorted_seq_list)

# create plot file
pdf(file=paste(args[1], ".pdf", sep = ""), width = 10, height = 10)

# plot the heatmap
ggplot(data = aa_dist_m_long, aes(x = Var1, y = Var2)) + 
  geom_tile(aes(fill=value)) + scale_fill_gradientn(colours=c("white", "black"), na.value = "grey", limits=c(94.9,100), oob=squish) + 
  coord_fixed() +
  theme_void() + 
  scale_y_discrete(expand = c(.02, 0)) + scale_x_discrete(expand = c(.02, 0)) + 
  theme(panel.border=element_rect(fill = NA, colour="black",size=2))

dev.off()
