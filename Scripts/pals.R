library(dplyr)
library(tidyr)
library(ggplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))
load("Processed_Data/divergent_classification.RData")

df_pals <- read.table(file="Processed_Data/pals.bed") 

system("bedtools intersect -a Processed_Data/pals.bed -b Processed_Data/divergent_regions/All_divergent_regions_clustered.bed -wao > Processed_Data/pals_div.bed")

df_pals_div <- read.table(file="Processed_Data/pals_div.bed") %>%
  dplyr::mutate(size=V3-V2) %>%
  dplyr::select(CHROM=V1, start=V2, stop=V3, size, div_overlap=V14, div_start=V12, div_stop=V13) %>%
  dplyr::mutate(frac_overlap=div_overlap/size)

save(df_pals_div, file="Processed_Data/pals_div.RData")

### 33/38(86.8% of pals genes are in the divergent region)

df_all_cluster_freq <- df_divergent_final %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, Frequency=freq) %>%
  dplyr::arrange(CHROM, START_BIN)

plot_div_pop_freq_pals <- df_all_cluster_freq %>%
  na.omit() %>%
  ggplot(.) +
  geom_segment(data=dplyr::distinct(df_pals_div,div_start, .keep_all=T), aes(x=(start+stop)/2e6, xend=(start+stop)/2e6, y=-0.01, yend=0.85), color='red', alpha=0.7, size=0.1) +
  geom_segment(data=data.frame(CHROM="I", start=2342215, stop=2356238), aes(x=start/1e6, xend=stop/1e6, y=-0.01, yend=0.85), color='blue', alpha=0.7, size=0.1) + ## peel-1, zeel-1
  geom_segment(data=data.frame(CHROM="III", start=11116993, stop=11136380), aes(x=start/1e6, xend=stop/1e6, y=-0.01, yend=0.85), color='blue', alpha=0.7, size=0.1) + ## sup-35, pha-1
  geom_rect(aes(xmin =  START_BIN/1e6, xmax = END_BIN/1e6, ymin = Frequency , ymax = Frequency+0.02), color='black', fill='black', size = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, color = "black"), axis.text.y = element_text(size=10, color = "black"), strip.text = element_text(size=10, color = "black"), legend.position = 'none', axis.title.x = element_text(size=11, color = "black"), axis.title.y = element_text(size=11, color = "black"), panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in"), panel.grid = element_blank()) +
  labs(x="Genomic position (Mb)", y="Frequency") +
  scale_y_continuous(expand = c(0, 0), limits=c(-0.01,0.85)) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  facet_grid(~CHROM, scales="free", space = 'free')

plot_div_pop_freq_pals
