library(stringr)

cbr_variant_counts <- NULL

for(strains in list.files("Processed_Data/cbriggsae/cbr_pi/variant_counts/")){
  
  strain_name <- str_split(strains, pattern = "\\.")[[1]][1]
  
  cts <- data.table::fread(glue::glue("Processed_Data/cbriggsae/cbr_pi/variant_counts/{strains}")) %>%
    dplyr::mutate(STRAIN = strain_name) %>%
    dplyr::rename(CHROM=V1, START_BIN = V2, END_BIN = V3, COUNT = V4)
  
  if(!exists("cbr_variant_counts")){
    cbr_variant_counts <- cts
   } else {
    cbr_variant_counts <- dplyr::bind_rows(cbr_variant_counts, cts)
  }
  
}

mean(cbr_variant_counts$COUNT)
median(cbr_variant_counts$COUNT)

### Mean adjusted threshold (11 (C. elegans) x mean ratio)

 ## cel mean count = 1.60856 cbr mean count = 2.496798, ratio = 2.496798/1.60856

ct_cbr <- 15*1.552195

cbr_nstrain <- length(unique(cbr_variant_counts$STRAIN))

cbr_df_chr_length <- cbr_variant_counts %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(stop=max(END_BIN)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(start=0)

###

cbr_population_all_div <- cbr_variant_counts  %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::mutate(div_id = ifelse(COUNT>ct_cbr, "Count", "Pass")) %>%
  dplyr::mutate(div_id = ifelse((grepl("Pass", dplyr::lag(div_id)) | grepl("Pass", dplyr::lead(div_id))) | !grepl("Pass", div_id), div_id, "TwoFlank")) %>%
  dplyr::mutate(div_id = ifelse(is.na(div_id), "Pass", div_id)) %>%
  dplyr::mutate(div_class = ifelse(div_id == "Pass", "Pass", "Divergent"))
  
  
###### filtering bins in >=xkb cluster, join xkb gaps ########

### joining ###

cluster_start <- NA
cluster_end <- NA

cbr_df_div_cluster <- cbr_population_all_div %>%
  dplyr::group_by(STRAIN, CHROM) %>%
  dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
  dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
  dplyr::ungroup() %>%
  dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)

cbr_df_div_cluster_filtered <- cbr_df_div_cluster %>%
  dplyr::filter(cluster == T)

cbr_cluster_start_vec <- cbr_df_div_cluster_filtered$cluster_start
cbr_start_bin_vec <- cbr_df_div_cluster_filtered$START_BIN

for (i in 1:length(cbr_cluster_start_vec)) {
  cbr_cluster_start_vec[i] <- ifelse(cbr_start_bin_vec[i]==0, 0, 
                                     ifelse(!is.na(cbr_cluster_start_vec[i]), cbr_cluster_start_vec[i], cbr_cluster_start_vec[i-1]))
  
}

cbr_df_div_cluster_filtered <- cbr_df_div_cluster_filtered %>%
  dplyr::select(-cluster_start) %>%
  dplyr::mutate(cluster_start=cbr_cluster_start_vec)

cbr_df_div_cluster_size <- cbr_df_div_cluster %>%
  dplyr::select(-cluster_start) %>%
  dplyr::left_join(., dplyr::select(cbr_df_div_cluster_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
  dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
  dplyr::mutate(cluster_size=ifelse(div_class =="Pass", NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cluster_end=cluster_start+cluster_size)

cbr_df_divergent <- cbr_df_div_cluster_size %>%
  #na.omit() %>%
  #dplyr::filter(div_class == "Divergent") %>%
  dplyr::filter(cluster_size >= 2e3) %>%
  dplyr::mutate(cluster_start=ifelse(is.na(cluster_start), START_BIN, cluster_start), 
                cluster_end =ifelse(is.na(cluster_end), END_BIN, cluster_end)) %>%
  dplyr::group_by(window_ID) %>%
  dplyr::mutate(freq=n()/cbr_nstrain, mwf = ifelse(freq > 0.5, 1-freq, freq)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq_bin=ifelse(freq < 0.01, "Rare", ifelse(freq < 0.05, "Intermediate", "Common")))

#write.table(cbr_df_divergent, file = "Processed_Data/cbr_df_divergent.tsv")

cbr_df_divergent <- read.table(file = "Processed_Data/cbr_df_divergent.tsv")

others <- c("NIC19","NIC20","QR24","QR25","BW287","EG4181","WG4181","QX1547","JU516","VX0034","ED3101","JU1341","JU1348")

cbr_df_divergent %>%
  dplyr::distinct(window_ID, freq,freq_bin) %>%
  dplyr::group_by(freq_bin) %>%
  dplyr::summarise(size=n()/1e3)

cbr_df_divergent_summary <- cbr_df_divergent %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  dplyr::mutate(group=ifelse(STRAIN %in% others,"Others","Tropical"))

cbr_div_summary <- cbr_df_divergent_summary %>%
  na.omit() %>%
  dplyr::group_by(STRAIN, group) %>%
  dplyr::summarise(sum = sum(cluster_size)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-sum) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(plotpoint = 1:n()) %>%
  dplyr::ungroup()

cbr_df_join_div <- cbr_df_divergent_summary %>%
  dplyr::select(-group) %>%
  dplyr::left_join(., cbr_div_summary, by = "STRAIN") 

cbr_df_all_cluster_size_freq <- cbr_df_divergent %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, freq) %>%
  dplyr::arrange(CHROM, START_BIN)

cbr_df_common_div <- cbr_df_all_cluster_size_freq %>%
  dplyr::filter(freq >= 0.1) %>%
  join_masks(.) %>%
  dplyr::mutate(cluster_start=ifelse(is.na(cluster_start), START_BIN, cluster_start), 
                cluster_end =ifelse(is.na(cluster_end), END_BIN, cluster_end)) %>%
  dplyr::distinct(CHROM, cluster_start, .keep_all=T) %>%
  dplyr::mutate(group = "Common")

cbr_df_join_div$group <- factor(cbr_df_join_div$group, levels = c("Tropical", "Others"))

cbr_plot_div_pop <- cbr_df_join_div %>%
  na.omit() %>%
  ggplot(.) +
  geom_rect(aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = plotpoint-0.5 , ymax = plotpoint+0.5), fill = 'black',color='black', size = 0.2) +
  #geom_segment(data=cbr_df_common_div, aes(x=(cluster_start+cluster_end)/2e6, xend=(cluster_start+cluster_end)/2e6, y=0.48, yend=36.52), color='brown', alpha=0.7, size=0.1) +
  #geom_rug(data=cbr_df_common_div, sides="b", color='brown', aes(x=(cluster_start+cluster_end)/2e6), length = unit(0.03, "npc"), size=0.3) +
  geom_hline(yintercept = 0, color='black', size=0.4) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y  = element_blank(), axis.title.y = element_text(size=11, color = "black"), strip.text = element_text(size=10, color = "black"), 
        legend.position = 'none', axis.title.x = element_blank(), axis.ticks.y = element_blank(), 
        strip.text.y=element_text(size=9, color = "black", angle = 180),  
        panel.spacing = unit(0.1, "lines"), panel.grid = element_blank()) +
  labs(x="Genomic position (Mb)", y="35 wild C. briggsae strains") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  #facet_grid(cluster~CHROM, scales = 'free', space = 'free') +
  facet_grid(group~CHROM, scales="free", switch = "y", space = 'free')

cbr_plot_div_pop

cbr_plot_div_pop_rug <- cbr_df_common_div %>%
  ggplot(.) +
  geom_rect(data=dplyr::mutate(cbr_df_join_div, group="Common"), aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin=0.2, ymax=0.5), 
            fill = 'transparent',color='transparent', size = 0.2) +
  geom_rect(color='brown', aes(xmin=cluster_start/1e6, xmax=cluster_end/1e6, ymin=0.2, ymax=0.5), size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y  = element_blank(), axis.title.y = element_text(size=11, color = "black"), 
        legend.position = 'none', axis.title.x = element_text(size=11, color = "black"), axis.ticks.y = element_blank(), 
        panel.spacing = unit(0.1, "lines"), panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank()) +
  labs(x="Genomic position (Mb)", y="") +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(group~CHROM, scales="free", switch = "y", space = 'free')
  
  cbr_plot_div_pop_rug
  
  ExFig6 <- cowplot::plot_grid(cbr_plot_div_pop, cbr_plot_div_pop_rug, ncol=1, rel_heights = c(10,1), align = "v")
  
  ExFig6
  
ggsave("Plots/Extended/cbr_div.png", ExFig6, width =  7.5, height = 6)

save(cbr_df_divergent, cbr_df_join_div, cbr_df_all_cluster_size_freq, cbr_df_common_div, file = "Processed_Data/Cbr_divergent.RData")
