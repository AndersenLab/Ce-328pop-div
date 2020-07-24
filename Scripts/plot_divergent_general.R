library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

load("Processed_Data/divergent_classification.RData")
load("Processed_Data/strain_geo.RData")

cluster_threshold = 9000

### statistics for divergent regions

plot_mask_stat <- df_divergent_final %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=div_id) +
  theme_bw() +
  theme(text = element_text(size=12, color='black'), axis.text.x = element_text(angle = 90, size = 12, color = 'black', vjust=0.5, hjust = 1), axis.title.x = element_blank()) +
  labs(y="Bin (1kb) counts")

plot_mask_stat

ggsave("Plots/plot_mask_stat.pdf", plot = plot_mask_stat, width = 7.5, height = 5, units = "in")  

## Frequency plot

df_divergent_final$freq_bin <- factor(df_divergent_final$freq_bin, levels = c("Rare", "Intermediate", "Common"))

plot_window_freq_bar <- df_divergent_final %>%
  dplyr::distinct(., window_ID, .keep_all = T) %>%
  dplyr::group_by(window_ID) %>%
  dplyr::ungroup() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=freq_bin) +
  theme_bw() +
  theme(axis.text=element_text(size=10, color='black'), axis.title=element_text(size=11, color='black')) +
  labs(x="Frequency", y="Size (kb)")

plot_window_freq_bar

ggsave("Plots/etc/plot_div_freq.pdf", plot_window_freq_bar, width =  4, height = 4, unit='in')

####### Fig. 2 Divergent regions #############

### Fig. 2a Divergent regions of 327 isotypes (no N2) - grouped by admixture groups

div_summary <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  na.omit() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(sum = sum(cluster_size)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-sum) %>%
  dplyr::rename(isotype=STRAIN) %>%
  dplyr::mutate(plotpoint=row_number())

median((na.omit(div_summary$sum)))/1e6  ## 1.639 Mb (median) 
mean((na.omit(div_summary$sum)))/1e6  ## 1.856 Mb (mean) 

df_join_div <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  na.omit() %>%
  dplyr::rename(isotype=STRAIN) %>%
  dplyr::left_join(., div_summary, by="isotype")

#save(df_join_div, file="Processed_Data/divergent_joined.RData")
load("Processed_Data/divergent_joined.RData")

df_chr_length <- data.table::fread("Processed_Data/chr_info.tsv") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::rename(CHROM = V1, stop= V2) %>%
  dplyr::mutate(start=0)

plot_div_pop <- df_join_div %>%
  na.omit() %>%
  ggplot(.) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=1.01), color='transparent', fill='transparent', size =0.1) +
  #geom_segment(data=dplyr::distinct(df_pals_div,div_start, .keep_all=T), aes(x=(start+stop)/2e6, xend=(start+stop)/2e6, y=0.4, yend=327.6), color='red', alpha=0.75, size=0.1) +
  #geom_segment(data=data.frame(CHROM="I", start=2342215, stop=2356238), aes(x=start/1e6, xend=stop/1e6, y=0.4, yend=327.6), color='blue', alpha=0.75, size=0.1) + ## peel-1, zeel-1
  #geom_segment(data=data.frame(CHROM="III", start=11116993, stop=11136380), aes(x=start/1e6, xend=stop/1e6, y=0.4, yend=327.6), color='purple', alpha=0.75, size=0.1) + ## sup-35, pha-1
  geom_rect(aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = plotpoint-0.5 , ymax = plotpoint+0.5), fill = 'black',color='black', size = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y  = element_blank(), strip.text = element_text(size=10, color = "black"), legend.position = 'none',axis.title.y = element_text(size=11, color = "black"), axis.title.x=element_blank(), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0.00, 0.00), limits=c(0.4,327.6)) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  facet_grid(~CHROM, scales="free",space = 'free') +
  labs(x="Genomic position (Mb)",y="327 wild isotypes")

plot_div_pop

ggsave("Plots/plot_div_pop.pdf", plot = plot_div_pop, width = 7.5, height = 5, units = "in")

### Fig. 2b Divergent region landscape (merged and joined) ### 

df_all_cluster_size_freq <- df_divergent_final %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, freq) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  join_masks(.) %>%
  na.omit() %>%
  dplyr::group_by(CHROM, cluster_start, cluster_end, cluster_size) %>%
  dplyr::summarise(Frequency = mean(freq)) %>%
  dplyr::ungroup()

mean(df_all_cluster_size_freq$Frequency)
sd(df_all_cluster_size_freq$Frequency)
mean(df_all_cluster_size_freq$cluster_size)
sd(df_all_cluster_size_freq$cluster_size)

plot_div_pop_freq <- df_all_cluster_size_freq %>%
  na.omit() %>%
  ggplot(.) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 0.2, ymax=0.3), color='transparent', fill='transparent', size =0.1) +
  geom_rect(aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = Frequency-0.01, ymax = Frequency+0.01), color='black', fill='black', size = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=10, color = "black"), legend.position = 'none', axis.title.x = element_text(size=11, color = "black"), axis.title.y = element_text(size=11, color = "black"), panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in"), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_blank()) +
  labs(x="Genomic position (Mb)", y="Frequency") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.02,0.7)) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  facet_grid(~CHROM, scales="free", space = 'free')

plot_div_pop_freq

### Fig. 2c Divergent regions size for each strain / fraiction variants in divergent regions

#### Fraction variants in divergent regions for every strain

population_all <- read.table(file = "Processed_Data/divergent_output_all_strains_all_bins.tsv")

df_count <- population_all %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::select(STRAIN, window_ID, COUNT)

mean(df_count$COUNT)

df_divergent_count <- df_div_cluster_gapjoined_size %>%
  dplyr::mutate(is_count = ifelse(div_id %in% c("Count", "CountLowCoverage"), 1, 0)) %>%
  dplyr::group_by(STRAIN, CHROM, cluster_start, cluster_end) %>%
  dplyr::mutate(flag_del = ifelse(max(is_count)==0, "Del", "Pass")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(divergent_cluster = ifelse(is.na(cluster_size), "Non-divergent", 
                                           ifelse(cluster_size >= cluster_threshold & flag_del == "Pass", "Divergent", "Non-divergent"))) %>%
  dplyr::left_join(., df_count, by=c('STRAIN', 'window_ID'))


df_divergent_count_summary <- df_divergent_count %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(total_bin_count = n(), total_variant_count = sum(COUNT)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN, divergent_cluster) %>%
  dplyr::summarise(bin_count = n(), variant_sum = sum(COUNT), variant_mean = mean(COUNT), 
                   total_bin_count = mean(total_bin_count), total_variant_count = mean(total_variant_count)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(div_fraction = bin_count/total_bin_count, variant_fraction = variant_sum/total_variant_count) %>%
  dplyr::filter(divergent_cluster == "Divergent")

save(df_divergent_count_summary, file="Processed_Data/divergent_count_summary.RData")

DataS6 <- df_divergent_count_summary %>%
  dplyr::mutate(Variant_density_genomewide=total_variant_count/total_bin_count) %>%
  dplyr::select(STRAIN, Divergent_size_kb=bin_count, Fraction_divergent=div_fraction, 
                Variant_count_all=total_variant_count, Varaint_count_divergent=variant_sum, Fraction_variant_divergent=variant_fraction,
                Variant_density_genomewide, Variant_density_divergent=variant_mean) %>%
  dplyr::mutate(density_fold_divergent = Variant_density_divergent/Variant_density_genomewide)

plot_divergent_strain_size_varfrac <- df_divergent_count_summary %>%
  ggplot(.) +
  geom_point(aes(x=div_fraction*100, y=variant_fraction*100), size = 0.6, alpha =0.8) +
  geom_text_repel(data=dplyr::filter(df_divergent_count_summary, variant_fraction>0.4), aes(x=div_fraction*100, y=variant_fraction*100, label = STRAIN), size = 2.5, box.padding = 0.2, point.padding = 0.1, segment.size = 0.02, nudge_x = -0.1, nudge_y= 0, segment.alpha = 0) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color='black'), axis.text = element_text(size=10, color='black'), panel.grid=element_blank()) +
  labs(x="Fraction of hyper-divergent regions (%)", y="Variant fraction (%)",plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in")) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))

plot_divergent_strain_size_varfrac

df_divergent_count_mean <- df_divergent_count %>%
  dplyr::group_by(divergent_cluster) %>%
  dplyr::summarise(mean_count = mean(COUNT)) %>%
  dplyr::ungroup() ## at least 19.1/1.28=14.92 times higher variant density

df_divergent_count_mean

### Fig. 2d Divergent regions x geography

df_divergent_strain <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  na.omit() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(sum_div = sum(cluster_size)) %>%
  dplyr::ungroup()

df_divergent_strain_geo <- df_divergent_strain %>%
  dplyr::rename(isotype = STRAIN) %>%
  dplyr::left_join(., indep_strain_info_geo, by ='isotype') %>%
  dplyr::filter(geo != "Unknown" & reference_strain == 1)

df_divergent_strain_geo$geo <- factor(df_divergent_strain_geo$geo, levels = c("Africa", "Asia", "Atlantic", "Austrailia", "Europe",
                                  "Hawaii", "New Zealand", "N. America", "S. America"))

plot_divergent_strain_geo_size <- df_divergent_strain_geo %>%
  ggplot(.) +
  geom_point(alpha = 0.6, position=position_jitterdodge(jitter.width = 3), size = 0.4, 
             aes(x=geo, y=sum_div/1e6, fill = geo)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.8, aes(x=geo, y=sum_div/1e6, fill = geo)) +
  #geom_text_repel(data=dplyr::filter(df_divergent_strain_geo, sum_div>62e5), aes(x=geo, y=sum_div/1e6, label = isotype), size = 2.5, box.padding = 0.2, point.padding = 0.15, segment.size = 0.02, segment.alpha = 0) +
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11, color='black'), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'), legend.position = 'none', 
        axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.6, hjust=1),
        plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in"), panel.grid = element_blank()) +
  labs(x="Origin", y="Hyper-divergent\nregions (Mb)") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12))

plot_divergent_strain_geo_size

Fig2ab <- cowplot::plot_grid(plot_div_pop, plot_div_pop_freq, 
                             labels = c("a","b"), rel_heights = c(2,1), ncol=1, label_size = 12,
                           align="v", axis='lr')

Fig2cd <- cowplot::plot_grid(plot_divergent_strain_size_varfrac, 
                             plot_divergent_strain_geo_size, 
                             labels = c("c","d"), 
                             rel_widths = c(1,1), ncol=2, label_size = 12)

Fig2 <- cowplot::plot_grid(Fig2ab, Fig2cd, 
                           rel_heights = c(2,1), ncol=1,
                             align="hv", axis='tblr')

#Fig2

ggsave("Plots/Main/Fig2.png", plot = Fig2,width = 7.5, height = 7.5, units = "in")


### Common divergent regions

df_div_common_fraction <- df_divergent_final %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(total_div_strain = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN, freq_bin) %>%
  dplyr::mutate(n_bin_strain=n(), frac_freq_strain=n()/total_div_strain) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN, freq_bin, total_div_strain, frac_freq_strain) %>%
  dplyr::filter(freq_bin == "Common")

df_divergent_strain_common_geo <- df_div_common_fraction %>%
  dplyr::rename(isotype = STRAIN) %>%
  dplyr::left_join(., indep_strain_info_geo, by ='isotype') %>%
  dplyr::filter(geo != "Unknown" & reference_strain == 1)

df_divergent_strain_common_geo$geo <- factor(df_divergent_strain_common_geo$geo, levels = c("Africa", "Asia", "Atlantic", "Austrailia", "Europe",
                                                                              "Hawaii", "New Zealand", "N. America", "S. America"))

plot_divergent_strain_geo_common <- df_divergent_strain_common_geo %>%
  ggplot(.) +
  geom_point(alpha = 0.6, position=position_jitterdodge(jitter.width = 3), size = 0.4, 
             aes(x=geo, y=frac_freq_strain*100, fill = geo)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.8, aes(x=geo, y=frac_freq_strain*100, fill = geo)) +
  #geom_text_repel(data=dplyr::filter(df_divergent_strain_geo, sum_div>62e5), aes(x=geo, y=sum_div/1e6, label = isotype), size = 2.5, box.padding = 0.2, point.padding = 0.15, segment.size = 0.02, segment.alpha = 0) +
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11, color='black'), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'), legend.position = 'none', 
        axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.6, hjust=1),
        plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in"), panel.grid = element_blank()) +
  labs(x="Origin", y="Common hyper-divergent regions (%)")

plot_divergent_strain_geo_common
