#!/usr/bin/env Rscript

# set up
library(ggplot2)
library(ggtree)
args = commandArgs(trailingOnly=TRUE)

# read in tree
tree <- read.tree(args[1])

# create plot file
pdf(file=paste(args[1], ".pdf", sep = ""), width = 10, height = 10)

# plot tree 
ggtree(tree, layout="equal_angle") + 
  coord_cartesian(xlim = c(-0.15, 0.15), ylim=c(-0.15,0.15)) + 
  geom_tiplab(aes(angle=angle)) + theme_tree2()

dev.off()
