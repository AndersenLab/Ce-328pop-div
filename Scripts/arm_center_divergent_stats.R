library(dplyr)
library(tidyr)

load("Processed_Data/divergent_classification.RData")

### all, arms, centers

df_tip_bins <- read.table(file="Processed_Data/df_bins_tips.bed") %>%
  dplyr::rename(CHROM=V1, START_BIN=V2, END_BIN=V3) %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN))

df_arm_bins <- read.table(file="Processed_Data/df_bins_arms.bed") %>%
  dplyr::rename(CHROM=V1, START_BIN=V2, END_BIN=V3) %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN))

df_center_bins <- read.table(file="Processed_Data/df_bins_centers.bed") %>%
  dplyr::rename(CHROM=V1, START_BIN=V2, END_BIN=V3) %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN))

### summary stats

df_div_stats_summary <- df_div_stats %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::group_by(window_ID, div_class) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(div_class) %>%
  dplyr::summarise(mean_count=mean(mean_count_bin), mean_cov = mean(mean_cov_bin)) %>%
  dplyr::ungroup()

df_div_stats_summary_arms <- df_div_stats %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_arm_bins$window_ID & CHROM != "X") %>%
  dplyr::group_by(window_ID, div_class) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(div_class) %>%
  dplyr::summarise(mean_count=mean(mean_count_bin), mean_cov = mean(mean_cov_bin)) %>%
  dplyr::ungroup()

#divergent vs non-divergent auto arms
21.348534/2.080217
#divergent vs non-divergent genome
21.348534/1.288661



### Non-divergent regions in arms

df_div_arms_nondiv <- df_div_stats %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_arm_bins$window_ID) %>%
  dplyr::mutate(class_num = ifelse(div_class=="Pass",0,1)) %>%
  dplyr::group_by(window_ID) %>%
  dplyr::filter(sum(class_num) == 0) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN)

readr::write_tsv(df_div_arms_nondiv, "Processed_Data/divergent_regions/df_arm_nondiv_coords.tsv")

## stats for each chromosome, arms and centers for each

df_div_stats_tips <- df_div_stats %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_tip_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Tip", class = "All") %>%
  dplyr::ungroup()

df_div_stats_arms <- df_div_stats %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_arm_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Arm", class = "All") %>%
  dplyr::ungroup()

df_div_stats_center <- df_div_stats %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_center_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Center", class = "All") %>%
  dplyr::ungroup()

df_div_stats_all <- df_div_stats %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "All", class = "All") %>%
  dplyr::ungroup()

## divergent stats

df_div_stats_tips_div <- df_div_stats %>%
  dplyr::filter(div_class == "Divergent") %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_tip_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Tip", class = "Hyper-divergent") %>%
  dplyr::ungroup()

df_div_stats_arms_div <- df_div_stats %>%
  dplyr::filter(div_class == "Divergent") %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_arm_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Arm", class = "Hyper-divergent") %>%
  dplyr::ungroup()

df_div_stats_center_div <- df_div_stats %>%
  dplyr::filter(div_class == "Divergent") %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_center_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Center", class = "Hyper-divergent") %>%
  dplyr::ungroup() 

df_div_stats_all_div <- df_div_stats %>%
  dplyr::filter(div_class == "Divergent") %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "All", class = "Hyper-divergent") %>%
  dplyr::ungroup()

## Non-divergent stats

df_div_stats_tips_nondiv <- df_div_stats %>%
  dplyr::filter(div_class != "Divergent") %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_tip_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Tip", class = "Non-divergent") %>%
  dplyr::ungroup()

df_div_stats_arms_nondiv <- df_div_stats %>%
  dplyr::filter(div_class != "Divergent") %>%  
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_arm_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Arm", class = "Non-divergent") %>%
  dplyr::ungroup()

df_div_stats_center_nondiv <- df_div_stats %>%
  dplyr::filter(div_class != "Divergent") %>%  
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_center_bins$window_ID) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "Center", class = "Non-divergent") %>%
  dplyr::ungroup() 

df_div_stats_all_nondiv <- df_div_stats %>%
  dplyr::filter(div_class != "Divergent") %>%  
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::group_by(CHROM, window_ID) %>%
  dplyr::summarise(mean_count_bin = mean(COUNT), mean_cov_bin = mean(fraction_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mean_count = mean(mean_count_bin), mean_cov = mean(mean_cov_bin), size=n_distinct(window_ID)) %>%
  dplyr::mutate(location = "All", class = "Non-divergent") %>%
  dplyr::ungroup()

## bind data.frame

df_div_stats_genome <- rbind(df_div_stats_tips, df_div_stats_arms,df_div_stats_center,df_div_stats_all, 
                             df_div_stats_tips_nondiv, df_div_stats_arms_nondiv,df_div_stats_center_nondiv,df_div_stats_all_nondiv,
                             df_div_stats_tips_div, df_div_stats_arms_div,df_div_stats_center_div,df_div_stats_all_div)

df_div_stats_genome_summary <- df_div_stats_genome %>%
  dplyr::group_by(CHROM, class, location) %>%
  dplyr::summarise(Variant_density = round(mean(mean_count),2), Mean_coverage = round(mean(mean_cov*100),2),
                   Size = round(mean(size)/1e3,2)) %>%
  dplyr::ungroup()

I_tip_div_NA <- c("I", "Hyper-divergent", "Tip", NA, NA, NA)

df_div_stats_genome_summary <- rbind(df_div_stats_genome_summary, I_tip_div_NA) %>%
  dplyr::arrange(CHROM, location)

write.csv(df_div_stats_genome_summary, file="Manuscript/Supplementary_Table/div_stats.csv")

df_div_stats_genome_summary_gather <- df_div_stats_genome_summary %>%
  tidyr::gather(-CHROM, -location, -class, key="stats", value="value")

df_div_stats_genome_summary_gather$stats <- factor(df_div_stats_genome_summary_gather$stats, 
                                                   levels = c("Variant_density", "Mean_coverage", "Size"), 
                                                    labels = c("Variant density (per kb)", "Mean coverage (%)", "Size"))
  
plot_div_chrom_comp <- df_div_stats_genome_summary_gather %>%
  dplyr::filter(stats != "Size") %>%
  ggplot(.) +
  geom_col(position = position_dodge(width = 1)) +
  aes(x=location, y=as.numeric(value), fill=class) +
  theme_bw() +
  scale_fill_brewer(palette = 'Set1') +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1, color='black'),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=10, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.position = c(0.87,0.41),
        legend.margin = margin(0,0,0,0))+
  facet_grid(stats~CHROM, scales = 'free') +
  labs(y="Statistics")

plot_div_chrom_comp

df_div_stats_genome_summary_gather_spread <- df_div_stats_genome_summary_gather %>%
  dplyr::mutate(class=ifelse(class=="Hyper-divergent", "Divergent", class)) %>%
  dplyr::mutate(class=ifelse(class=="Non-divergent", "Non_divergent", class)) %>%
  tidyr::spread(key=class, value=value) %>%
  dplyr::mutate(ratio=as.numeric(Divergent)/as.numeric(Non_divergent))

plot_fold_diff_div <- df_div_stats_genome_summary_gather_spread %>%
  dplyr::filter(stats != "Size") %>%
  ggplot(.) +
  geom_col() +
  aes(x=location, y=as.numeric(ratio)) +
  theme_bw() +
  #scale_fill_brewer(palette = 'Set1') +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1, color='black'),
        axis.title.y = element_text(size=11, color='black'), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=10, color='black'),
        legend.text = element_text(size=10, color='black'))+
  labs(y="Fold difference") +
  facet_grid(stats~CHROM, scales = 'free')

FigS9 <- cowplot::plot_grid(plot_div_chrom_comp, plot_fold_diff_div, 
                            labels = c("a","b"), ncol=1, 
                            align = 'v', rel_heights = c(1,1), axis = 'l', label_size = 12)

FigS9

ggsave(FigS9, file = "Plots/supplementary_figures/FigS9.png", width=7.5, height=9, unit = 'in')

save(df_div_stats_genome, df_div_stats_genome_summary_gather, df_div_stats_genome_summary_gather_spread, file="Processed_Data/Divergent_chromosome_stats.RData")

### size summary

df_div_stats_genome_summary_gather %>%
  na.omit() %>%
  #dplyr::filter(CHROM != "X") %>% ## for autosomal arms only
  dplyr::filter(class == "Hyper-divergent" & stats == "Size") %>%
  dplyr::group_by(location) %>%
  dplyr::summarise(sum_size=sum(as.numeric(value))) %>%
  dplyr::ungroup()

14.2/20.5

as.numeric(df_div_stats_genome_summary_gather$value)

###############

df_div_stats_arms_auto <- df_div_stats_arms %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::filter(window_ID %in% df_arm_bins$window_ID & CHROM != "X") %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(mean_count = mean(COUNT), mean_cov = mean(fraction_cov)) %>%
  dplyr::ungroup()

mean(df_div_stats_arms_auto$mean_count)

df_div_stats_arms_auto_div <- df_div_stats_arms_auto %>%
  dplyr::filter(div_class == "Divergent") %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(mean_count = mean(COUNT), mean_cov = mean(fraction_cov)) %>%
  dplyr::ungroup()

mean(df_div_stats_iso_div$mean_count)
mean(df_div_stats_iso_div$mean_cov)

