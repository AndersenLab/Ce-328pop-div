library(dplyr)
library(tidyr)
library(ggplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/divergent_classification.RData")

div_CB4856_all <- df_div_cluster_gapjoined_size %>%
  dplyr::filter(STRAIN == "CB4856")

df_pacbio <- read_tsv(file="Processed_Data/longread_align/CB4856_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

df_CB4856_freq <- df_divergent_final %>%
  dplyr::filter(STRAIN == "CB4856") %>%
  dplyr::select(CHROM, START_BIN, END_BIN, freq, mwf, freq_bin)

div_CB4856_all_select <- div_CB4856_all %>%
  dplyr::select(CHROM, START_BIN, END_BIN, div_id) %>%
  dplyr::mutate(div_id=ifelse(is.na(div_id), "Pass", as.character(div_id))) %>%
  dplyr::left_join(., df_CB4856_freq, by=c('CHROM', 'START_BIN', 'END_BIN'))

div_CB4856_all_select$div_id <- factor(div_CB4856_all_select$div_id, levels = c("Count", "LowCoverage","CountLowCoverage","TwoFlank", "gap","Pass"), labels = c("Count", "Low coverage","Count & Low coverage","Two flank", "Gap","Pass"))

class.colors <- c("Pass"="#D9FFFC", "Count"="#FF0009", "Low coverage"="#0047ab", "Count & Low coverage"="#C24CF6","Two flank"="maroon", "Gap"="yellow3")

div_longread_freq <- function(df1=df_pacbio, df2_1=div_CB4856_all_select, df2_2=df_CB4856_freq, chr="V",left=3331000, right=3446000, ymin=55, guide_rownum=1) {
  
  df2_1_filter <- df2_1 %>%
    dplyr::filter(CHROM == chr, START_BIN >=left-3e4, END_BIN <= right+3e4)
  
  plot <- df1 %>%
    dplyr::filter(CHROM == chr) %>%
    ggplot(.) +
    geom_rect(data=dplyr::filter(df1, N2_length > 1e3 & CHROM == chr), aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity, ymin = Identity-1)) +
    geom_rect(data=df2_1_filter, aes(xmin=START_BIN/1e6, xmax=END_BIN/1e6, ymax=102, ymin = 101, fill = div_id), color = 'black', size = 0.02) +
    geom_rect(data=dplyr::filter(df2_2, CHROM == chr), aes(xmin=START_BIN/1e6, xmax=END_BIN/1e6, ymax=103.3, ymin = 103), size=0.3, color = "darkorange1") +
    scale_fill_manual(values = class.colors) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
          axis.title=element_text(size=11, color='black'), 
          legend.text = element_text(size=9, color='black'), legend.margin = margin(0,0,0,0), 
          legend.direction = "horizontal", legend.box = "horizontal", legend.position = c(0.5, 0.15), 
          legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.05, "in")) +
    guides(fill= guide_legend(nrow=guide_rownum, order = 1, title.position = "left", direction = "horizontal", byrow=T)) +
    labs(x="Genomic position (Mb)", y="Identity (%)", fill = "") +
    scale_x_continuous(expand = c(0,0)) +
    ylim(ymin,103.5) +
    facet_grid(~CHROM) +
    coord_cartesian(xlim=c((left-3e4)/1e6, (right+3e4)/1e6))
  
  return(plot)
}

div_longread_freq()

CB4856_divergent <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  dplyr::filter(STRAIN == "CB4856")

save(df_pacbio, div_CB4856_all_select,df_CB4856_freq,class.colors,div_longread_freq, file = "Processed_Data/CB4856_div_longread.RData")

## plot fraction of common divergent regions in CB4856

plot_window_freq_bar_CB4856 <- df_divergent_final %>%
  dplyr::filter(STRAIN == "CB4856") %>%  
  dplyr::distinct(., window_ID, .keep_all = T) %>%
  dplyr::group_by(window_ID) %>%
  dplyr::ungroup() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=freq_bin) +
  theme_bw() +
  theme(axis.text=element_text(size=11, color='black'), axis.title=element_text(size=12, color='black')) +
  labs(x="Frequency", y="Size in CB4856 (kb)")

plot_window_freq_bar_CB4856

ggsave("Plots/plot_div_freq_CB4856.pdf", plot_window_freq_bar_CB4856, width =  4, height = 3)


for (i in 1:nrow(CB4856_divergent)) {
  chr <- as.character(CB4856_divergent$CHROM[i])
  left <- as.numeric(CB4856_divergent$cluster_start[i])
  right <- as.numeric(CB4856_divergent$cluster_end[i])
  div_longread_freq(chr=chr, left=left, right=right)
  ggsave(glue::glue("Plots/CB4856_divergent_region/CB4856_freq/{chr}_{left}-{right}.pdf"), width = 9, height = 5.5, units = "in")
}


div_longread_freq2 <- function(df1=df_pacbio, df2=div_CB4856_all_select, chr = "V", left=3331000, right=3446000) {
  
  plot1 <- df1 %>%
    dplyr::filter(CHROM == chr) %>%
    ggplot(.) +
    geom_rect(data=dplyr::filter(df1, N2_length > 1e3 & CHROM == chr), aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity, ymin = Identity-1)) +
    geom_rect(data=dplyr::filter(df2, CHROM == chr), aes(xmin=START_BIN/1e6, xmax=END_BIN/1e6, ymax=102, ymin = 101, fill = div_id)) +
    scale_fill_manual(values = class.colors) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size=11, color='black'), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y=element_text(size=12, color='black'), legend.text = element_text(size=9, color='black'), legend.direction = "horizontal", legend.box = "horizontal", legend.position = c(0.5,0.08), legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.2, "in")) +
    guides(fill= guide_legend(nrow=1, order = 1, title.position = "left")) +
    labs(x="Genomic position (Mb)", y="Identity (%)", fill = "") +
    ylim(60,102.5) +
    facet_grid(~CHROM) +
    coord_cartesian(xlim=c((left-5e4)/1e6, (right+5e4)/1e6))
  
  plot2 <- df2 %>%
    dplyr::filter(CHROM == chr) %>%
    dplyr::filter(!is.na(freq_bin)) %>%
    ggplot(.) +
    geom_rect(aes(xmin=START_BIN/1e6, xmax=END_BIN/1e6, ymax=mwf+0.01, ymin = mwf, fill = freq_bin), size=0.5) +
    scale_fill_manual(values = freq.colors) +
    geom_vline(xintercept = left/1e6, color = 'red') +
    geom_vline(xintercept = right/1e6, color = 'red') +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(size=11, color='black'), axis.title=element_text(size=12, color='black'), strip.text = element_blank(), strip.background = element_blank(), legend.title=element_blank(), legend.text = element_text(size=9, color='black'), legend.direction = "vertical", legend.box = "vertical", legend.position = c(0.9,0.87), legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.2, "in"), legend.background = element_rect(fill = "white", colour = "white")) +
    labs(x="Genomic position (Mb)", y="Fraction of strains with \ndivergent genotype") +
    scale_y_continuous(limits = c(0,0.53), expand = c(0,0.01), breaks = c(0,0.1,0.2,0.3,0.4,0.5)) +
    facet_grid(~CHROM) +
    coord_cartesian(xlim=c((left-5e4)/1e6, (right+5e4)/1e6))
  
  plot <- cowplot::plot_grid(plot1, plot2, ncol=1, rel_heights = c(3,2), labels = c("a", "b"), align = "v", label_size = 12)  
  
  return(plot)
}

div_longread_freq2()

for (i in 1:nrow(CB4856_divergent)) {
  chr <- as.character(CB4856_divergent$CHROM[i])
  left <- as.numeric(CB4856_divergent$cluster_start[i])
  right <- as.numeric(CB4856_divergent$cluster_end[i])
  div_longread_freq2(chr=chr, left=left, right=right)
  ggsave(glue::glue("Plots/CB4856_divergent_region/CB4856_freq2/{chr}_{left}-{right}.pdf"), width = 9, height = 7, units = "in")
}