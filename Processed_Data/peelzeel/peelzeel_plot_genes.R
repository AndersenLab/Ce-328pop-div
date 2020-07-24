library(tidyverse)
library(scales)

N2_records <- read.table("N2_haplotype.tsv", header=TRUE)
CB4856_records <- read.table("CB4856_haplotype.tsv", header=TRUE)

CB4856_transform = -2.32

#pdf(file="I_2318291-2381851_gene_plot.pdf", width=13.03, height=4.53)

plot <- 
  ggplot() +
  geom_hline(yintercept=-0.5) + 
  geom_hline(yintercept=-2) +
  geom_rect(data=N2_records, aes(xmin=start/1e6, xmax=stop/1e6, ymax=-1, ymin=0, fill=sc_ortho_status), color="black") + 
  geom_rect(data=CB4856_records, aes(xmin=start/1e6-CB4856_transform, xmax=stop/1e6-CB4856_transform, ymax=-2.5, ymin=-1.5, fill=sc_ortho_status), color="black") +
  geom_rect(aes(xmin=2329000/1e6, xmax=2361000/1e6, ymax=-0.9, ymin=-0.1), fill="black") +
  theme_void() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("grey", "ORANGE"))

plot

#dev.off()
