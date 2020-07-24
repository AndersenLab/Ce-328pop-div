library(tidyverse)
library(scales)

N2_haplotype <- read.table("N2_haplotype.tsv", header=TRUE)
JU2526_haplotype <- read.table("JU2526_haplotype.tsv", header=TRUE)
EG4725_haplotype <- read.table("EG4725_haplotype.tsv", header=TRUE)
QX1794_haplotype <- read.table("QX1794_haplotype.tsv", header=TRUE)
DL238_haplotype <- read.table("DL238_haplotype.tsv", header=TRUE) 
ECA36_haplotype <- read.table("ECA36_haplotype.tsv", header=TRUE)
NIC526_haplotype <- read.table("NIC526_haplotype.tsv", header=TRUE)

JU2526_transform = -20.1905
EG4725_transform = -20.186
ECA36_transform = -20.189
DL238_transform = -20.188
QX1794_transform = -20.190
NIC526_transform = -20.1885

#pdf(file="I_2318291-2381851_gene_plot.pdf", width=13.03, height=4.53)

plot <- 
  ggplot() +
  geom_hline(yintercept=-0.5) + 
  geom_hline(yintercept=-2) +
  geom_hline(yintercept=-3.5) +
  geom_hline(yintercept=-5) +
  geom_hline(yintercept=-6.5) + 
  geom_hline(yintercept=-8) + 
  geom_hline(yintercept=-9.5) + 
  geom_rect(data=N2_haplotype, aes(xmin=start/1e6, xmax=stop/1e6, ymax=-1, ymin=0, fill=sc_ortho_status), color="black") + 
  # divergent regions
  geom_rect(aes(xmin=20187000/1e6, xmax=20234000/1e6, ymax=-1, ymin=0), alpha=0.5) + 
  geom_rect(data=JU2526_haplotype, aes(xmin=start/1e6-JU2526_transform, xmax=stop/1e6-JU2526_transform, ymax=-2.5, ymin=-1.5, fill=sc_ortho_status), color="black") +
  geom_rect(data=EG4725_haplotype, aes(xmin=start/1e6-EG4725_transform, xmax=stop/1e6-EG4725_transform, ymax=-4, ymin=-3, fill=sc_ortho_status), color="black") +
  geom_rect(data=ECA36_haplotype, aes(xmin=start/1e6-ECA36_transform, xmax=stop/1e6-ECA36_transform, ymax=-5.5, ymin=-4.5, fill=sc_ortho_status), color="black") +
  geom_rect(data=DL238_haplotype, aes(xmin=start/1e6-DL238_transform, xmax=stop/1e6-DL238_transform, ymax=-7, ymin=-6, fill=sc_ortho_status), color="black") +
  geom_rect(data=QX1794_haplotype, aes(xmin=start/1e6-QX1794_transform, xmax=stop/1e6-QX1794_transform, ymax=-8.5, ymin=-7.5, fill=sc_ortho_status), color="black") +
  geom_rect(data=NIC526_haplotype, aes(xmin=start/1e6-NIC526_transform, xmax=stop/1e6-NIC526_transform, ymax=-10, ymin=-9, fill=sc_ortho_status), color="black") +
  theme_void() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("grey", "ORANGE"))
plot

#dev.off()
