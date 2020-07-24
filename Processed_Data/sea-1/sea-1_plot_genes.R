library(tidyverse)
library(scales)

N2_records <- read.table("N2_haplotype.tsv", header=TRUE)
CB4856_records <- read.table("CB4856_haplotype.tsv", header=TRUE)
DL238_records <- read.table("DL238_haplotype.tsv", header=TRUE)

DL238_transform = -2.415
CB4856_transform = 0.13

pdf(file="II_3667179-3701405_gene_plot.pdf", width=13.03, height=8.53)

plot <- 
  ggplot() +
  geom_hline(yintercept=-0.5) + 
  geom_hline(yintercept=-2) +
  geom_hline(yintercept=-3.5) +
  geom_rect(data=N2_records, aes(xmin=start/1e6, xmax=stop/1e6, ymax=-1, ymin=0, fill=sc_ortho_status), color="black") + 
  geom_rect(data=CB4856_records, aes(xmin=start/1e6-CB4856_transform, xmax=stop/1e6-CB4856_transform, ymax=-2.5, ymin=-1.5, fill=sc_ortho_status), color="black") + 
  geom_rect(data=DL238_records, aes(xmin=start/1e6-DL238_transform, xmax=stop/1e6-DL238_transform, ymax=-4, ymin=-3, fill=sc_ortho_status), color="black") +
  theme_void() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("grey", "#303030"))

plot

dev.off()
