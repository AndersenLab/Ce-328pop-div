library(dplyr)
library(ggplot2)
library(valr)
library(readr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

join_masks <- function(mask_file = NULL){
  cluster <- NA
  cluster_start <- NA
  
  for (i in c(1:nrow(mask_file)-1)) {
    
    if(i==1){
      
      cluster[i] <- ifelse(mask_file$END_BIN[i] == mask_file$START_BIN[i+1], 'yes', 'no')
      cluster_start[i] <- mask_file$START_BIN[i]
      
    } else {
      
      cluster[i] <- ifelse(mask_file$END_BIN[i-1] == mask_file$START_BIN[i] | mask_file$END_BIN[i] == mask_file$START_BIN[i+1], 'yes', 'no')
      cluster_start[i] <- ifelse(mask_file$CHROM[i] == mask_file$CHROM[i-1] & cluster[i] == 'yes',
                                 ifelse(mask_file$END_BIN[i-1] == mask_file$START_BIN[i] | mask_file$END_BIN[i] == mask_file$START_BIN[i+1],
                                        ifelse(cluster[i-1] == "no", mask_file$START_BIN[i],
                                               ifelse(mask_file$END_BIN[i-1] == mask_file$START_BIN[i], cluster_start[i-1], mask_file$START_BIN[i])), mask_file$START_BIN[i]), mask_file$START_BIN[i])
      
    }
  }
  
  cluster[nrow(mask_file)] <- NA
  cluster_start[nrow(mask_file)] <- NA
  
  mask_file_cluster <- data.frame(mask_file, cluster, cluster_start) %>%
    dplyr::group_by(CHROM, cluster_start) %>%
    dplyr::mutate(cluster_size=ifelse(is.na(cluster_start), 1000, n()*1000)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_end=cluster_start+cluster_size)
  
  
}
divergent_gap <- function(df=df_joined_masks_freq, threshold = cluster_threshold, gap_threshold=gap) {
  
  df_gap <- df %>%
    dplyr::filter(cluster_size >= threshold) %>%
    dplyr::distinct(STRAIN, CHROM, cluster_start,cluster_size, cluster_end)
  
  gap <- NA
  gap_size <- NA
  gap_join <- NA
  extended_cluster_end <- NA
  extended_cluster_size <-NA
  
  for (i in c(1:nrow(df_gap)-1)) {
    
    gap[i] <- ifelse(df_gap$STRAIN[i] == df_gap$STRAIN[i+1] & df_gap$CHROM[i] == df_gap$CHROM[i+1], 'yes', NA)
    gap_size[i] <- ifelse(gap[i]=="yes", df_gap$cluster_start[i+1]-df_gap$cluster_end[i], NA)
    gap_join[i] <- ifelse(gap[i]=="yes" & gap_size[i] <= gap_threshold & gap_size[i] > 0, 'join', NA)
    extended_cluster_end[i] <- ifelse(gap[i]=="yes" & gap_size[i] <= gap_threshold & gap_size[i] > 0, df_gap$cluster_end[i+1],df_gap$cluster_end[i])
    extended_cluster_size[i] <- ifelse(gap[i]=="yes" & gap_size[i] <= gap_threshold & gap_size[i] > 0, extended_cluster_end[i]-df_gap$cluster_start[i], df_gap$cluster_size[i])
    
  }
  
  gap[nrow(df_gap)] <- NA
  gap_size[nrow(df_gap)] <- NA
  gap_join[nrow(df_gap)] <- NA
  extended_cluster_end[nrow(df_gap)] <- NA
  extended_cluster_size[nrow(df_gap)] <- NA
  
  df_gap_size <- data.frame(df_gap, gap, gap_size, gap_join, extended_cluster_end, extended_cluster_size) %>%
    dplyr::mutate(threshold = threshold, threshold_pass = ifelse(cluster_size >= threshold, T, F)) %>%
    dplyr::group_by(STRAIN, CHROM, extended_cluster_end) %>%
    dplyr::mutate(extended_cluster_start = min(cluster_start)) %>%
    dplyr::mutate(extended_cluster_size = extended_cluster_end - extended_cluster_start) %>%
    dplyr::ungroup()
  
  df_gap_size_join <- df %>%
    dplyr::left_join(., df_gap_size, by=c('STRAIN', 'CHROM', 'cluster_start', 'cluster_size', 'cluster_end'))
  
  
  return(list(df_gap_size, df_gap_size_join))
  
}


longread_strains <- c("CB4856", "DL238", "ECA36", "ECA396", "EG4725", "JU310", "JU1400", "JU2526", "JU2600", "MY2147", "MY2693", "NIC2", "NIC526","QX1794", "XZ1516")

df_accuracy <- NULL

for (strain in longread_strains) {
  
  Shortread_raw <- data.table::fread(glue::glue("Processed_Data/Processed_Masks_nohet_DELINS20kb/{strain}_Mask_DF.tsv")) %>%
    dplyr::select(CHROM, START_BIN, END_BIN, STRAIN, COUNT, COVERAGE, fraction_SV_bases)
  
  Shortread_raw_coverage <- Shortread_raw %>%
    dplyr::group_by(STRAIN) %>%
    dplyr::summarise(mean_cov_genome = mean(COVERAGE)) %>%
    dplyr::ungroup()
  
  Shortread_coverage_fraction <- Shortread_raw %>%
    dplyr::left_join(., Shortread_raw_coverage, by = "STRAIN") %>%
    dplyr::mutate(fraction_cov = COVERAGE/mean_cov_genome)
  
   Longread_raw <- read_tsv(file=glue::glue("Processed_Data/longread_align/{strain}_vs_N2.coords")) %>%
    dplyr::rename(Ref_start = '[S1]', Ref_end = '[E1]', Wild_start = '[S2]', Wild_end = '[E2]', 
                  Ref_length = '[LEN 1]', Wild_length = '[LEN 2]', Identity = "[% IDY]", 
                  Ref_chrom_length = '[LEN R]', Wild_chrom_length='[LEN Q]', CHROM='[TAGS]')
  
   Longread_raw_ID <- Longread_raw %>%
     dplyr::select(chrom=CHROM, start=Ref_start, end=Ref_end, Identity)
   
   df_bins <- read_tsv(file="Processed_Data/df_bins.tsv") %>%
     dplyr::rename(chrom=chr, end=stop)
   
   df_bins_longread <- valr::bed_intersect(df_bins,Longread_raw_ID) %>%
     dplyr::rename(CHROM=chrom, START_BIN=start.x, END_BIN=end.x, 
                   start=start.y, stop=end.y, coverage=.overlap, Identity=Identity.y) %>%
     dplyr::mutate(full_cov = ifelse(coverage == 1e3, T, F)) %>% 
     dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
     dplyr::mutate(max_coverage = max(coverage), max_ID=max(Identity)) %>%
     dplyr::filter(coverage == max_coverage) %>%
     dplyr::mutate(bin_ID = max(Identity)) %>%
     dplyr::ungroup() %>%
     dplyr::mutate(bin_ID_normalized = max(Identity)*coverage/1e3)
   
   df_divergent_comp <- Shortread_coverage_fraction %>%
     dplyr::left_join(., df_bins_longread, by=c("CHROM", "START_BIN","END_BIN")) %>%
     dplyr::mutate(bin_ID_normalized = ifelse(is.na(bin_ID), 0, bin_ID), 
                   max_coverage = ifelse(is.na(bin_ID), 0, max_coverage)) %>%
     dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
     dplyr::mutate(n=n()) %>%
     dplyr::ungroup() %>%
     dplyr::distinct(CHROM, START_BIN, END_BIN, bin_ID, .keep_all=T)
  
   for (id in c(95)) {
     
     for (long_size in c(9)) {
       
       for (long_cov in c(60)) {
         
         df_divergent_comp_div <- df_divergent_comp %>%
           dplyr::mutate(div_class = ifelse(max_ID < id | is.na(max_ID) | max_coverage <= long_cov*10, "div", "Pass")) %>%
           dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
           dplyr::filter(div_class == "div" | div_twoflank == "TwoFlank") %>%
           dplyr::select(CHROM, START_BIN, END_BIN, bin_ID, max_coverage, bin_ID_normalized, COUNT, fraction_cov, div_class, div_twoflank) %>%
           dplyr::rename(longread_coverage=max_coverage, shortread_depth = fraction_cov) %>%
           join_masks(.) %>%
           dplyr::mutate(div_cluster=ifelse(cluster_size >= long_size*1e3, 'Divergent', 'Non-divergent'))
         
   for (cov in c(0.25,0.3,0.35,0.4,0.45)) {
     
     for (ct in c(9,12,15,18,21)) {  
      
   shortread_div <- Shortread_coverage_fraction  %>%
     dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
     dplyr::mutate(div_count = ifelse(COUNT>ct, "Count", "Pass"), div_lowcov = ifelse(fraction_cov<cov, "LowCoverage", "Pass")) %>%
     dplyr::mutate(div_class = ifelse(paste(div_count, div_lowcov, sep="_") == "Pass_Pass", "Pass", "Divergent")) %>%
     dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
     dplyr::mutate(div_id = paste(ifelse(div_count=="Pass","",div_count), ifelse(div_lowcov=="Pass","",div_lowcov), ifelse(div_twoflank=="Pass","",div_twoflank), sep=""), div_class = ifelse(paste(div_count, div_lowcov, div_twoflank, sep="_") == "Pass_Pass_Pass", "Pass", "Divergent"))
   
   cluster_start <- NA
   cluster_end <- NA
   
   df_div_cluster <- shortread_div %>%
     dplyr::group_by(STRAIN, CHROM) %>%
     dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
     dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
     dplyr::ungroup() %>%
     dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)
   
   df_div_cluster_filtered <- df_div_cluster %>%
     dplyr::filter(cluster == T)
   
   cluster_start_vec <- df_div_cluster_filtered$cluster_start
   
   for (i in 1:length(cluster_start_vec)) {
     cluster_start_vec[i] <- ifelse(!is.na(cluster_start_vec[i]), cluster_start_vec[i], cluster_start_vec[i-1])
   }
   
   df_div_cluster_filtered <- df_div_cluster_filtered %>%
     dplyr::select(-cluster_start) %>%
     dplyr::mutate(cluster_start=cluster_start_vec)
   
   df_div_cluster_size <- df_div_cluster %>%
     dplyr::select(-cluster_start) %>%
     dplyr::left_join(., dplyr::select(df_div_cluster_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
     dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
     dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
     dplyr::ungroup() %>%
     dplyr::mutate(cluster_end=cluster_start+cluster_size)
   
   for (x in c(6,9,12,15)) {
     
     clu_threshold = 1000*x
      
   df_gap_bind_list <- divergent_gap(df=df_div_cluster_size, threshold = clu_threshold, gap_threshold=3000)
   
   df_gap_bind <- df_gap_bind_list[[1]]
   df_gap_bind_all_bins <- df_gap_bind_list[[2]] %>%
     dplyr::mutate(div_id = ifelse(div_class == "Divergent", div_id, ifelse(dplyr::lag(END_BIN, 1) == dplyr::lag(cluster_end, 1) & dplyr::lag(gap_join, 1) == 'join' | dplyr::lag(END_BIN, 2) == dplyr::lag(cluster_end, 2) & dplyr::lag(gap_join, 2) == "join" | dplyr::lag(END_BIN, 3) == dplyr::lag(cluster_end, 3) & dplyr::lag(gap_join, 3) == "join", "gap_3kb", div_id))) %>%
     dplyr::mutate(div_id = ifelse(is.na(div_id),"", div_id)) %>%
     dplyr::mutate(div_class = ifelse(div_id == "gap_3kb", "Divergent", div_class)) %>%
     dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id)
   
   df_div_cluster_gapjoined <- df_gap_bind_all_bins %>%
     dplyr::group_by(STRAIN, CHROM) %>%
     dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
     dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
     dplyr::ungroup() %>%
     dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)
   
   df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined %>%
     dplyr::filter(cluster == T)
   
   cluster_start_gapjoined_vec <- df_div_cluster_gapjoined_filtered$cluster_start
   
   for (i in 1:length(cluster_start_gapjoined_vec)) {
     cluster_start_gapjoined_vec[i] <- ifelse(!is.na(cluster_start_gapjoined_vec[i]), cluster_start_gapjoined_vec[i], cluster_start_gapjoined_vec[i-1])
   }
   
   df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined_filtered %>%
     dplyr::select(-cluster_start) %>%
     dplyr::mutate(cluster_start=cluster_start_gapjoined_vec)
   
   df_div_cluster_gapjoined_size <- df_div_cluster_gapjoined %>%
     dplyr::select(-cluster_start) %>%
     dplyr::left_join(., dplyr::select(df_div_cluster_gapjoined_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
     dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
     dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
     dplyr::ungroup() %>%
     dplyr::mutate(cluster_end=cluster_start+cluster_size)
   
   df_divergent_final_select <- df_div_cluster_gapjoined_size %>%
     dplyr::mutate(is_count = ifelse(div_id %in% c("Count", "CountLowCoverage"), 1, 0)) %>%
     dplyr::group_by(STRAIN, CHROM, cluster_start, cluster_end) %>%
     dplyr::mutate(flag_del = ifelse(max(is_count)==0, "Del", "Pass")) %>%
     dplyr::ungroup() %>%
     dplyr::select(CHROM, START_BIN, END_BIN, div_class, div_id, cluster_start, cluster_end, cluster_size,flag_del)
   
   df_final <- df_divergent_final_select %>%
     dplyr::mutate(class = ifelse(cluster_size >=clu_threshold & flag_del == "Pass", "Divergent", "Non-divergent")) %>%
     dplyr::mutate(class=ifelse(is.na(class), "Non-divergent", class)) %>%
     dplyr::filter(class == "Divergent")
   
   
   ## compare with long-read divergrent regions
   
   df_final_select_short <- df_final %>%
     dplyr::select(CHROM, START_BIN, END_BIN, cluster_start, cluster_end) %>%
     dplyr::mutate(source="shortread")
   
   df_final_select_long <- df_divergent_comp_div %>%
     dplyr::filter(div_cluster == "Divergent") %>%
     dplyr::select(CHROM, START_BIN, END_BIN, cluster_start, cluster_end) %>%
     dplyr::mutate(source="longread")
   
   df_final_select_comp <- rbind(df_final_select_short, df_final_select_long) %>%
     dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
     dplyr::mutate(n=n()) %>%
     dplyr::ungroup() %>%
     dplyr::mutate(class=ifelse(n==2, "both", ifelse(source == "shortread", "short_only", "long_only"))) %>%
     dplyr::distinct(CHROM, START_BIN, END_BIN, class) %>%
     dplyr::arrange(CHROM, START_BIN, END_BIN) %>%
     join_masks(.)
   
   accuracy <- nrow(dplyr::filter(df_final_select_comp, class=="both"))/nrow(df_final_select_comp)   
   
   accuracy
  
   div_comp_parameter <- data.frame(isotype=strain, 
                                          long_id_threshold = id, long_cluster_threshold = long_size,
                                          long_coverage_threshold = long_cov,
                                          cov_threshold = cov, count_threshold =ct, cluster_threshold = clu_threshold,
                                          both=nrow(dplyr::filter(df_final_select_comp, class=="both")), 
                                          short_only=nrow(dplyr::filter(df_final_select_comp, class=="short_only")),
                                          long_only=nrow(dplyr::filter(df_final_select_comp, class=="long_only")),
                                          accuracy=nrow(dplyr::filter(df_final_select_comp, class=="both"))/nrow(df_final_select_comp))
   
   df_accuracy <- rbind(df_accuracy, div_comp_parameter)
   
   }
     }
   }
       }
     }
   }
}

df_accuracy_ <- df_accuracy %>%
  dplyr::mutate(accuracy_short=both/(both+short_only))

write.table(na.omit(df_accuracy_), file = "Processed_Data/df_accuracy.tsv", quote=F, row.names=F)

df_accuracy_ <- read.table(file = "Processed_Data/df_accuracy.tsv", header=T)

df_accuracy_summarise <- df_accuracy_ %>%
  dplyr::distinct(isotype, long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  dplyr::group_by(long_id_threshold, long_coverage_threshold, long_cluster_threshold, count_threshold, cov_threshold,cluster_threshold) %>%
  dplyr::summarise(sum_both = sum(both), sum_short_only = sum(short_only), sum_long_only = sum(long_only), 
                   mean_accuracy=mean(accuracy), sd_accuracy=sd(accuracy),
                   mean_accuracy_short=mean(accuracy_short)) %>%
  dplyr::ungroup()
  
plot_accuracy_count <- df_accuracy_ %>%
  na.omit() %>%
  dplyr::distinct(isotype, long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  dplyr::filter(cluster_threshold==9e3) %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.5) +
  aes(x=(both+short_only+long_only)/1e3, y=accuracy*100, fill=count_threshold) +
  scale_fill_continuous(breaks=c(3,6,9,12,15,18)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'), 
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.y= element_text(size=10, color='black'), 
        axis.text.x = element_text(size=10, color='black'),
        legend.position='right')+
  facet_wrap(~isotype, ncol=3) +
  labs(x="Total hyper-divergent regions detected (Mb)", 
       y="Overlap of hyper-divergent regions between \nlong-read and short-read based approach (%)", 
       fill = "Count\nthreshold")

plot_accuracy_count

plot_accuracy_coverage <- df_accuracy_ %>%
  na.omit() %>%
  dplyr::distinct(isotype, long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  dplyr::filter(cluster_threshold==9e3) %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.5) +
  aes(x=(both+short_only+long_only)/1e3, y=accuracy*100, fill=cov_threshold) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'), 
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.y= element_text(size=10, color='black'), 
        axis.text.x = element_text(size=10, color='black'),
        legend.position='right')+
  facet_wrap(~isotype, ncol=3) +
  labs(x="Total hyper-divergent regions detected (Mb)", 
       y="Overlap of hyper-divergent regions between \nlong-read and short-read based approach (%)", 
       fill = "Coverage\nthreshold")

plot_accuracy_coverage



#ggsave(file="Long_strains_long_short_comp_twoflank_fine1.png", width=7.5, height=5, unit='in')

df_accuracy_ %>%
  na.omit() %>%
  dplyr::distinct(isotype, long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.5) +
  aes(x=both/1e3, y=accuracy*100, fill=isotype) +
  theme_bw() +
  labs(x="Long+/short+ divergent regions (Mb)", y="Long+/short+ divergent regions (%)")

#ggsave(file="Long_strains_long_short_comp_twoflank_fine2.png", width=7.5, height=5, unit='in')

df_accuracy_ %>%
  na.omit() %>%
  dplyr::distinct(isotype, long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  ggplot(.) +
  geom_point(shape=21, alpha =0.5) +
  aes(x=long_only/1e3, y=short_only/1e3, fill=isotype) +
  theme_bw() +
  labs(x="long-read only (Mb)", y="short-read only (Mb)", fill="accuracy")

#ggsave(file="Long_strains_long_short_comp_twoflank_fine3.png", width=7.5, height=5, unit='in')

df_accuracy_ %>%
  na.omit() %>%
  dplyr::distinct(isotype, long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.5) +
  aes(x=both/1e3, y=accuracy_short*100, fill=isotype) +
  theme_bw() +
  labs(x="Short+ divergent regions (Mb)", y="Long+/short+ divergent regions (%)") 

#ggsave(file="Long_strains_long_short_comp_twoflank_fine4.png", width=7.5, height=5, unit='in')

df_accuracy_summarise %>%
  na.omit() %>%
  dplyr::distinct(long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  ggplot(.) +
  geom_point(alpha=0.5) +
  aes(x=(sum_both+sum_short_only+sum_long_only)/1e3, y=mean_accuracy*100) +
  theme_bw() +
  labs(x="Total divergent regions detected (Mb)", y="Average long+/short+ divergent regions (%)",
       fill = "Long+/short+\ndivergent\nregions (Mb)")

ggsave(file="Long_strains_long_short_comp_twoflank_fine5.png", width=7.5, height=5, unit='in')

df_accuracy_summarise %>%
  na.omit() %>%
  dplyr::distinct(long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, .keep_all=T) %>%
  ggplot(.) +
  geom_point(alpha=0.5) +
  aes(x=(sum_both+sum_short_only)/1e3, y=mean_accuracy_short*100) +
  theme_bw() +
  labs(x="Short+ divergent regions (Mb)", y="Average long+/short+ divergent regions (%)")  +
  geom_vline(xintercept = 37259/1e3, size=0.2, color='red') +
  geom_hline(yintercept = 0.8445994*100, size = 0.2, color='red')+
  geom_vline(xintercept = 38890/1e3, size=0.2, color='blue') +
  geom_hline(yintercept = 0.8759818*100, size = 0.2, color='blue')

ggsave(file="Long_strains_long_short_comp_twoflank_fine6.png", width=7.5, height=5, unit='in')


### apply new gap filling to all strains

df_accuracy_regap <- NULL

for (strain in longread_strains) {
  
  Shortread_raw <- data.table::fread(glue::glue("Processed_Data/Processed_Masks_nohet_DELINS20kb/{strain}_Mask_DF.tsv")) %>%
    dplyr::select(CHROM, START_BIN, END_BIN, STRAIN, COUNT, COVERAGE, fraction_SV_bases)
  
  Shortread_raw_coverage <- Shortread_raw %>%
    dplyr::group_by(STRAIN) %>%
    dplyr::summarise(mean_cov_genome = mean(COVERAGE)) %>%
    dplyr::ungroup()
  
  Shortread_coverage_fraction <- Shortread_raw %>%
    dplyr::left_join(., Shortread_raw_coverage, by = "STRAIN") %>%
    dplyr::mutate(fraction_cov = COVERAGE/mean_cov_genome)
  
  Longread_raw <- read_tsv(file=glue::glue("Processed_Data/longread_align/{strain}_vs_N2.coords")) %>%
    dplyr::rename(Ref_start = '[S1]', Ref_end = '[E1]', Wild_start = '[S2]', Wild_end = '[E2]', 
                  Ref_length = '[LEN 1]', Wild_length = '[LEN 2]', Identity = "[% IDY]", 
                  Ref_chrom_length = '[LEN R]', Wild_chrom_length='[LEN Q]', CHROM='[TAGS]')
  
  Longread_raw_ID <- Longread_raw %>%
    dplyr::select(chrom=CHROM, start=Ref_start, end=Ref_end, Identity)
  
  df_bins <- read_tsv(file="Processed_Data/df_bins.tsv") %>%
    dplyr::rename(chrom=chr, end=stop)
  
  df_bins_longread <- valr::bed_intersect(df_bins,Longread_raw_ID) %>%
    dplyr::rename(CHROM=chrom, START_BIN=start.x, END_BIN=end.x, 
                  start=start.y, stop=end.y, coverage=.overlap, Identity=Identity.y) %>%
    dplyr::mutate(full_cov = ifelse(coverage == 1e3, T, F)) %>% 
    dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
    dplyr::mutate(max_coverage = max(coverage), max_ID=max(Identity)) %>%
    dplyr::filter(coverage == max_coverage) %>%
    dplyr::mutate(bin_ID = max(Identity)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bin_ID_normalized = max(Identity)*coverage/1e3)
  
  df_divergent_comp <- Shortread_coverage_fraction %>%
    dplyr::left_join(., df_bins_longread, by=c("CHROM", "START_BIN","END_BIN")) %>%
    dplyr::mutate(bin_ID_normalized = ifelse(is.na(bin_ID), 0, bin_ID), 
                  max_coverage = ifelse(is.na(bin_ID), 0, max_coverage)) %>%
    dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
    dplyr::mutate(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(CHROM, START_BIN, END_BIN, bin_ID, .keep_all=T)
  
  for (id in c(95)) {
    
    for (long_size in c(9)) {
      
      for (long_cov in c(60)) {
        
        df_divergent_comp_div <- df_divergent_comp %>%
          dplyr::mutate(div_class = ifelse(max_ID < id | is.na(max_ID) | max_coverage <= long_cov*10, "div", "Pass")) %>%
          dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
          dplyr::filter(div_class == "div" | div_twoflank == "TwoFlank") %>%
          dplyr::select(CHROM, START_BIN, END_BIN, bin_ID, max_coverage, bin_ID_normalized, COUNT, fraction_cov, div_class, div_twoflank) %>%
          dplyr::rename(longread_coverage=max_coverage, shortread_depth = fraction_cov) %>%
          join_masks(.) %>%
          dplyr::mutate(div_cluster=ifelse(cluster_size >= long_size*1e3, 'Divergent', 'Non-divergent'))
        
        for (cov in c(0.35)) {
          
          for (ct in c(15)) {  
            
            shortread_div <- Shortread_coverage_fraction  %>%
              dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
              dplyr::mutate(div_count = ifelse(COUNT>ct, "Count", "Pass"), div_lowcov = ifelse(fraction_cov<cov, "LowCoverage", "Pass")) %>%
              dplyr::mutate(div_class = ifelse(paste(div_count, div_lowcov, sep="_") == "Pass_Pass", "Pass", "Divergent")) %>%
              dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
              dplyr::mutate(div_id = paste(ifelse(div_count=="Pass","",div_count), ifelse(div_lowcov=="Pass","",div_lowcov), ifelse(div_twoflank=="Pass","",div_twoflank), sep=""), div_class = ifelse(paste(div_count, div_lowcov, div_twoflank, sep="_") == "Pass_Pass_Pass", "Pass", "Divergent"))
            
            cluster_start <- NA
            cluster_end <- NA
            
            df_div_cluster <- shortread_div %>%
              dplyr::group_by(STRAIN, CHROM) %>%
              dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
              dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
              dplyr::ungroup() %>%
              dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)
            
            df_div_cluster_filtered <- df_div_cluster %>%
              dplyr::filter(cluster == T)
            
            cluster_start_vec <- df_div_cluster_filtered$cluster_start
            
            for (i in 1:length(cluster_start_vec)) {
              cluster_start_vec[i] <- ifelse(!is.na(cluster_start_vec[i]), cluster_start_vec[i], cluster_start_vec[i-1])
            }
            
            df_div_cluster_filtered <- df_div_cluster_filtered %>%
              dplyr::select(-cluster_start) %>%
              dplyr::mutate(cluster_start=cluster_start_vec)
            
            df_div_cluster_size <- df_div_cluster %>%
              dplyr::select(-cluster_start) %>%
              dplyr::left_join(., dplyr::select(df_div_cluster_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
              dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
              dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
              dplyr::ungroup() %>%
              dplyr::mutate(cluster_end=cluster_start+cluster_size)
            
            for (clu_threshold in c(9e3)) {
            for (gap_cluster in c(9)) {
              for (gap_size in c(2,3,4,5,6,7,8,9)) {
                
                gap_cluster_threshold = gap_cluster*1e3
                gap_size_threshold = gap_size*1e3
                
                df_gap_bind_list <- divergent_gap(df=df_div_cluster_size, threshold = gap_cluster_threshold, gap_threshold=gap_size_threshold)
                
                df_gap_bind <- df_gap_bind_list[[1]]
                df_gap_bind_all_bins <- df_gap_bind_list[[2]] %>%
                  dplyr::mutate(div_id = ifelse(div_class == "Divergent", div_id, 
                                                ifelse(dplyr::lag(END_BIN, 1) == dplyr::lag(cluster_end, 1) & dplyr::lag(gap_join, 1) == 'join' | 
                                                         dplyr::lag(END_BIN, 2) == dplyr::lag(cluster_end, 2) & dplyr::lag(gap_join, 2) == "join" | 
                                                         dplyr::lag(END_BIN, 3) == dplyr::lag(cluster_end, 3) & dplyr::lag(gap_join, 3) == "join" |
                                                         dplyr::lag(END_BIN, 4) == dplyr::lag(cluster_end, 4) & dplyr::lag(gap_join, 4) == "join" |
                                                         dplyr::lag(END_BIN, 5) == dplyr::lag(cluster_end, 5) & dplyr::lag(gap_join, 5) == "join" | 
                                                         dplyr::lag(END_BIN, 6) == dplyr::lag(cluster_end, 6) & dplyr::lag(gap_join, 6) == "join" |
                                                         dplyr::lag(END_BIN, 7) == dplyr::lag(cluster_end, 7) & dplyr::lag(gap_join, 7) == "join" |
                                                         dplyr::lag(END_BIN, 8) == dplyr::lag(cluster_end, 8) & dplyr::lag(gap_join, 8) == "join", "gap", div_id))) %>%
                  dplyr::mutate(div_id = ifelse(is.na(div_id),"", div_id)) %>%
                  dplyr::mutate(div_class = ifelse(div_id == "gap", "Divergent", div_class)) %>%
                  dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id)
                
              df_div_cluster_gapjoined <- df_gap_bind_all_bins %>%
                dplyr::group_by(STRAIN, CHROM) %>%
                dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
                dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
                dplyr::ungroup() %>%
                dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)
              
              df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined %>%
                dplyr::filter(cluster == T)
              
              cluster_start_gapjoined_vec <- df_div_cluster_gapjoined_filtered$cluster_start
              
              for (i in 1:length(cluster_start_gapjoined_vec)) {
                cluster_start_gapjoined_vec[i] <- ifelse(!is.na(cluster_start_gapjoined_vec[i]), cluster_start_gapjoined_vec[i], cluster_start_gapjoined_vec[i-1])
              }
              
              df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined_filtered %>%
                dplyr::select(-cluster_start) %>%
                dplyr::mutate(cluster_start=cluster_start_gapjoined_vec)
              
              df_div_cluster_gapjoined_size <- df_div_cluster_gapjoined %>%
                dplyr::select(-cluster_start) %>%
                dplyr::left_join(., dplyr::select(df_div_cluster_gapjoined_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
                dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
                dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(cluster_end=cluster_start+cluster_size)
              
              df_divergent_final_select <- df_div_cluster_gapjoined_size %>%
                dplyr::mutate(is_count = ifelse(div_id %in% c("Count", "CountLowCoverage"), 1, 0)) %>%
                dplyr::group_by(STRAIN, CHROM, cluster_start, cluster_end) %>%
                dplyr::mutate(flag_del = ifelse(max(is_count)==0, "Del", "Pass")) %>%
                dplyr::ungroup() %>%
                dplyr::select(CHROM, START_BIN, END_BIN, div_class, div_id, cluster_start, cluster_end, cluster_size,flag_del)
              
              df_final <- df_divergent_final_select %>%
                dplyr::mutate(class = ifelse(cluster_size >=clu_threshold & flag_del == "Pass", "Divergent", "Non-divergent")) %>%
                dplyr::mutate(class=ifelse(is.na(class), "Non-divergent", class)) %>%
                dplyr::filter(class == "Divergent")
              
              
              ## compare with long-read divergrent regions
              
              df_final_select_short <- df_final %>%
                dplyr::select(CHROM, START_BIN, END_BIN, cluster_start, cluster_end) %>%
                dplyr::mutate(source="shortread")
              
              df_final_select_long <- df_divergent_comp_div %>%
                dplyr::filter(div_cluster == "Divergent") %>%
                dplyr::select(CHROM, START_BIN, END_BIN, cluster_start, cluster_end) %>%
                dplyr::mutate(source="longread")
              
              df_final_select_comp <- rbind(df_final_select_short, df_final_select_long) %>%
                dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
                dplyr::mutate(n=n()) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(class=ifelse(n==2, "both", ifelse(source == "shortread", "short_only", "long_only"))) %>%
                dplyr::distinct(CHROM, START_BIN, END_BIN, class) %>%
                dplyr::arrange(CHROM, START_BIN, END_BIN) %>%
                join_masks(.)
              
              accuracy <- nrow(dplyr::filter(df_final_select_comp, class=="both"))/nrow(df_final_select_comp)   
              
              accuracy
              
              div_comp_parameter <- data.frame(isotype=strain, 
                                               long_id_threshold = id, long_cluster_threshold = long_size,
                                               long_coverage_threshold = long_cov,
                                               cov_threshold = cov, count_threshold =ct, cluster_threshold = clu_threshold,
                                               gap_cluster=gap_cluster, gap_size=gap_size,
                                               both=nrow(dplyr::filter(df_final_select_comp, class=="both")), 
                                               short_only=nrow(dplyr::filter(df_final_select_comp, class=="short_only")),
                                               long_only=nrow(dplyr::filter(df_final_select_comp, class=="long_only")),
                                               accuracy=nrow(dplyr::filter(df_final_select_comp, class=="both"))/nrow(df_final_select_comp))
              
              df_accuracy_regap <- rbind(df_accuracy_regap, div_comp_parameter)
              
              }
            }
          }
        }
      }
    }
  }
  }
}


df_accuracy_regap_ <- df_accuracy_regap %>%
  dplyr::mutate(accuracy_short=both/(both+short_only))

write.table(na.omit(df_accuracy_regap_), file = "Processed_Data/df_accuracy_regap.tsv", quote=F, row.names=F)

df_accuracy_regap_summarise <- df_accuracy_regap_ %>%
  dplyr::distinct(isotype, long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, gap_cluster, gap_size, .keep_all=T) %>%
  dplyr::group_by(long_id_threshold, long_coverage_threshold, long_cluster_threshold, count_threshold, cov_threshold,cluster_threshold, gap_cluster, gap_size) %>%
  dplyr::summarise(sum_both = sum(both), sum_short_only = sum(short_only), sum_long_only = sum(long_only), 
                   mean_accuracy=mean(accuracy), sd_accuracy=sd(accuracy),
                   mean_accuracy_short=mean(accuracy_short)) %>%
  dplyr::ungroup()

df_accuracy_regap_summarise %>%
  na.omit() %>%
  dplyr::distinct(long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, gap_cluster, gap_size, .keep_all=T) %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.7, size=2) +
  aes(x=(sum_both+sum_short_only+sum_long_only)/1e3, y=mean_accuracy*100, fill=gap_size) +
  theme_bw() +
  labs(x="Total divergent regions detected (Mb)", y="Long+/short+ divergent regions (%)")

ggsave(file="Long_strains_long_short_comp_twoflank_gap1.png", width=5, height=5, unit='in')


df_accuracy_regap_summarise %>%
  na.omit() %>%
  dplyr::distinct(long_id_threshold, long_cluster_threshold,
                  long_coverage_threshold,
                  cov_threshold, count_threshold, cluster_threshold, gap_cluster, gap_size, .keep_all=T) %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.7) +
  aes(x=(sum_both+sum_short_only)/1e3, y=mean_accuracy_short*100, fill=gap_size) +
  theme_bw() +
  labs(x="Short+ divergent regions (Mb)", y="Long+/short+ divergent regions (%)")

ggsave(file="Long_strains_long_short_comp_twoflank_gap2.png", width=7.5, height=5, unit='in')



######## final application

longread_strains <- c("CB4856", "DL238", "ECA36", "ECA396", "EG4725", "JU310", "JU1400", "JU2526", "JU2600", "MY2147", "MY2693", "NIC2", "NIC526","QX1794", "XZ1516")

id=95
long_size=9
long_cov=60
cov=0.35
ct=15
clu_threshold=9000
gap_cluster_threshold=9000
gap_size_threshold=5000

df_final_select_comp_all <- NULL
df_accuracy_final <- NULL
df_bins_stats <- NULL

for (strain in longread_strains) {
  
  Shortread_raw <- data.table::fread(glue::glue("Processed_Data/Processed_Masks_nohet_DELINS20kb/{strain}_Mask_DF.tsv")) %>%
    dplyr::select(CHROM, START_BIN, END_BIN, STRAIN, COUNT, COVERAGE, fraction_SV_bases)
  
  Shortread_raw_coverage <- Shortread_raw %>%
    dplyr::group_by(STRAIN) %>%
    dplyr::summarise(mean_cov_genome = mean(COVERAGE)) %>%
    dplyr::ungroup()
  
  Shortread_coverage_fraction <- Shortread_raw %>%
    dplyr::left_join(., Shortread_raw_coverage, by = "STRAIN") %>%
    dplyr::mutate(fraction_cov = COVERAGE/mean_cov_genome)
  
  Longread_raw <- read_tsv(file=glue::glue("Processed_Data/longread_align/{strain}_vs_N2.coords")) %>%
    dplyr::rename(Ref_start = '[S1]', Ref_end = '[E1]', Wild_start = '[S2]', Wild_end = '[E2]', 
                  Ref_length = '[LEN 1]', Wild_length = '[LEN 2]', Identity = "[% IDY]", 
                  Ref_chrom_length = '[LEN R]', Wild_chrom_length='[LEN Q]', CHROM='[TAGS]')
  
  Longread_raw_ID <- Longread_raw %>%
    dplyr::select(chrom=CHROM, start=Ref_start, end=Ref_end, Identity)
  
  df_bins <- read_tsv(file="Processed_Data/df_bins.tsv") %>%
    dplyr::rename(chrom=chr, end=stop)
  
  df_bins_longread <- valr::bed_intersect(df_bins,Longread_raw_ID) %>%
    dplyr::rename(CHROM=chrom, START_BIN=start.x, END_BIN=end.x, 
                  start=start.y, stop=end.y, coverage=.overlap, Identity=Identity.y) %>%
    dplyr::mutate(full_cov = ifelse(coverage == 1e3, T, F)) %>% 
    dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
    dplyr::mutate(max_coverage = max(coverage), max_ID=max(Identity)) %>%
    dplyr::filter(coverage == max_coverage) %>%
    dplyr::mutate(bin_ID = max(Identity)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bin_ID_normalized = max(Identity)*coverage/1e3)
  
  df_divergent_comp <- Shortread_coverage_fraction %>%
    dplyr::left_join(., df_bins_longread, by=c("CHROM", "START_BIN","END_BIN")) %>%
    dplyr::mutate(bin_ID_normalized = ifelse(is.na(bin_ID), 0, bin_ID), 
                  max_coverage = ifelse(is.na(bin_ID), 0, max_coverage)) %>%
    dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
    dplyr::mutate(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(CHROM, START_BIN, END_BIN, bin_ID, .keep_all=T)
  
        df_divergent_comp_div <- df_divergent_comp %>%
          dplyr::mutate(div_class = ifelse(max_ID < id | is.na(max_ID) | max_coverage <= long_cov*10, "div", "Pass")) %>%
          dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
          dplyr::filter(div_class == "div" | div_twoflank == "TwoFlank") %>%
          dplyr::select(CHROM, START_BIN, END_BIN, bin_ID, max_coverage, bin_ID_normalized, COUNT, fraction_cov, div_class, div_twoflank) %>%
          dplyr::rename(longread_coverage=max_coverage, shortread_depth = fraction_cov) %>%
          join_masks(.) %>%
          dplyr::mutate(div_cluster=ifelse(cluster_size >= long_size*1e3, 'Divergent', 'Non-divergent'))
        
         
            shortread_div <- Shortread_coverage_fraction  %>%
              dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
              dplyr::mutate(div_count = ifelse(COUNT>ct, "Count", "Pass"), div_lowcov = ifelse(fraction_cov<cov, "LowCoverage", "Pass")) %>%
              dplyr::mutate(div_class = ifelse(paste(div_count, div_lowcov, sep="_") == "Pass_Pass", "Pass", "Divergent")) %>%
              dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
              dplyr::mutate(div_id = paste(ifelse(div_count=="Pass","",div_count), ifelse(div_lowcov=="Pass","",div_lowcov), ifelse(div_twoflank=="Pass","",div_twoflank), sep=""), div_class = ifelse(paste(div_count, div_lowcov, div_twoflank, sep="_") == "Pass_Pass_Pass", "Pass", "Divergent"))
            
            cluster_start <- NA
            cluster_end <- NA
            
            df_div_cluster <- shortread_div %>%
              dplyr::group_by(STRAIN, CHROM) %>%
              dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
              dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
              dplyr::ungroup() %>%
              dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)
            
            df_div_cluster_filtered <- df_div_cluster %>%
              dplyr::filter(cluster == T)
            
            cluster_start_vec <- df_div_cluster_filtered$cluster_start
            
            for (i in 1:length(cluster_start_vec)) {
              cluster_start_vec[i] <- ifelse(!is.na(cluster_start_vec[i]), cluster_start_vec[i], cluster_start_vec[i-1])
            }
            
            df_div_cluster_filtered <- df_div_cluster_filtered %>%
              dplyr::select(-cluster_start) %>%
              dplyr::mutate(cluster_start=cluster_start_vec)
            
            df_div_cluster_size <- df_div_cluster %>%
              dplyr::select(-cluster_start) %>%
              dplyr::left_join(., dplyr::select(df_div_cluster_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
              dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
              dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
              dplyr::ungroup() %>%
              dplyr::mutate(cluster_end=cluster_start+cluster_size)
      
                df_gap_bind_list <- divergent_gap(df=df_div_cluster_size, threshold = gap_cluster_threshold, gap_threshold=gap_size_threshold)
                
                df_gap_bind <- df_gap_bind_list[[1]]
                df_gap_bind_all_bins <- df_gap_bind_list[[2]] %>%
                  dplyr::mutate(div_id = ifelse(div_class == "Divergent", div_id, 
                                                ifelse(dplyr::lag(END_BIN, 1) == dplyr::lag(cluster_end, 1) & dplyr::lag(gap_join, 1) == 'join' | 
                                                         dplyr::lag(END_BIN, 2) == dplyr::lag(cluster_end, 2) & dplyr::lag(gap_join, 2) == "join" | 
                                                         dplyr::lag(END_BIN, 3) == dplyr::lag(cluster_end, 3) & dplyr::lag(gap_join, 3) == "join" |
                                                         dplyr::lag(END_BIN, 4) == dplyr::lag(cluster_end, 4) & dplyr::lag(gap_join, 4) == "join" |
                                                         dplyr::lag(END_BIN, 5) == dplyr::lag(cluster_end, 5) & dplyr::lag(gap_join, 5) == "join" | 
                                                         dplyr::lag(END_BIN, 6) == dplyr::lag(cluster_end, 6) & dplyr::lag(gap_join, 6) == "join" |
                                                         dplyr::lag(END_BIN, 7) == dplyr::lag(cluster_end, 7) & dplyr::lag(gap_join, 7) == "join" |
                                                         dplyr::lag(END_BIN, 8) == dplyr::lag(cluster_end, 8) & dplyr::lag(gap_join, 8) == "join", "gap", div_id))) %>%
                  dplyr::mutate(div_id = ifelse(is.na(div_id),"", div_id)) %>%
                  dplyr::mutate(div_class = ifelse(div_id == "gap", "Divergent", div_class)) %>%
                  dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id)
                
                df_div_cluster_gapjoined <- df_gap_bind_all_bins %>%
                  dplyr::group_by(STRAIN, CHROM) %>%
                  dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
                  dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
                  dplyr::ungroup() %>%
                  dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)
                
                df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined %>%
                  dplyr::filter(cluster == T)
                
                cluster_start_gapjoined_vec <- df_div_cluster_gapjoined_filtered$cluster_start
                
                for (i in 1:length(cluster_start_gapjoined_vec)) {
                  cluster_start_gapjoined_vec[i] <- ifelse(!is.na(cluster_start_gapjoined_vec[i]), cluster_start_gapjoined_vec[i], cluster_start_gapjoined_vec[i-1])
                }
                
                df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined_filtered %>%
                  dplyr::select(-cluster_start) %>%
                  dplyr::mutate(cluster_start=cluster_start_gapjoined_vec)
                
                df_div_cluster_gapjoined_size <- df_div_cluster_gapjoined %>%
                  dplyr::select(-cluster_start) %>%
                  dplyr::left_join(., dplyr::select(df_div_cluster_gapjoined_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
                  dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
                  dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(cluster_end=cluster_start+cluster_size)
                
                df_divergent_final_select <- df_div_cluster_gapjoined_size %>%
                  dplyr::mutate(is_count = ifelse(div_id %in% c("Count", "CountLowCoverage"), 1, 0)) %>%
                  dplyr::group_by(STRAIN, CHROM, cluster_start, cluster_end) %>%
                  dplyr::mutate(flag_del = ifelse(max(is_count)==0, "Del", "Pass")) %>%
                  dplyr::ungroup() %>%
                  dplyr::select(CHROM, START_BIN, END_BIN, div_class, div_id, cluster_start, cluster_end, cluster_size,flag_del)
                
                df_final <- df_divergent_final_select %>%
                  dplyr::mutate(class = ifelse(cluster_size >=clu_threshold & flag_del == "Pass", "Divergent", "Non-divergent")) %>%
                  dplyr::mutate(class=ifelse(is.na(class), "Non-divergent", class)) %>%
                  dplyr::filter(class == "Divergent")
                
                ## compare with long-read divergrent regions
                
                df_final_select_short <- df_final %>%
                  dplyr::select(CHROM, START_BIN, END_BIN, cluster_start, cluster_end) %>%
                  dplyr::mutate(source="shortread")
                
                df_final_select_long <- df_divergent_comp_div %>%
                  dplyr::filter(div_cluster == "Divergent") %>%
                  dplyr::select(CHROM, START_BIN, END_BIN, cluster_start, cluster_end) %>%
                  dplyr::mutate(source="longread")
                
                df_final_select_comp <- rbind(df_final_select_short, df_final_select_long) %>%
                  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
                  dplyr::mutate(n=n()) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(class=ifelse(n==2, "both", ifelse(source == "shortread", "short_only", "long_only"))) %>%
                  dplyr::distinct(CHROM, START_BIN, END_BIN, class) %>%
                  dplyr::arrange(CHROM, START_BIN, END_BIN) %>%
                  join_masks(.) %>%
                  dplyr::mutate(isotype=strain)
                
                df_final_select_comp_all <- rbind(df_final_select_comp_all, df_final_select_comp)
                
                div_comp_parameter <- data.frame(isotype=strain, 
                                                 long_id_threshold = id, long_cluster_threshold = long_size,
                                                 long_coverage_threshold = long_cov,
                                                 cov_threshold = cov, count_threshold =ct, cluster_threshold = clu_threshold,
                                                 gap_size=gap_size_threshold,
                                                 both=nrow(dplyr::filter(df_final_select_comp, class=="both")), 
                                                 short_only=nrow(dplyr::filter(df_final_select_comp, class=="short_only")),
                                                 long_only=nrow(dplyr::filter(df_final_select_comp, class=="long_only")),
                                                 accuracy=nrow(dplyr::filter(df_final_select_comp, class=="both"))/nrow(df_final_select_comp))
                
                df_accuracy_final <- rbind(df_accuracy_final, div_comp_parameter)

}     

write.table(na.omit(df_final_select_comp_all), file = "Processed_Data/df_final_select_comp_all.tsv", quote=F, row.names=F)

df_accuracy_final_ <- df_accuracy_final %>%
  dplyr::mutate(accuracy_short=both/(both+short_only))

write.table(na.omit(df_accuracy_final_), file = "Processed_Data/df_accuracy_final.tsv", quote=F, row.names=F)

nrow(dplyr::distinct(dplyr::filter(df_final_select_comp_all, class != "long_only"), CHROM, START_BIN))

df_accuracy_final_ <- read.table(file = "Processed_Data/df_accuracy_final.tsv", header=T)

df_accuracy_final_ %>%
  na.omit() %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8, size = 3) +
  aes(x=(both+short_only+long_only)/1e3, y=accuracy*100, fill=isotype) +
  theme_bw() +
  labs(x="Total divergent regions detected (Mb)", y="Long+/short+ divergent regions (%)")

ggsave(file="final_long_short_comp1.png", width=7.5, height=5, unit='in')

df_accuracy_final_ %>%
  na.omit() %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8, size=3) +
  aes(x=both/1e3, y=accuracy_short*100, fill=isotype) +
  theme_bw() +
  labs(x="Short+ divergent regions (Mb)", y="Long+/short+ divergent regions (%)") 

ggsave(file="final_long_short_comp2.png", width=7.5, height=5, unit='in')


df_final_select_comp_all <- read.table(file = "Processed_Data/df_final_select_comp_all.tsv", header=T)

                
df_final_select_comp_all %>%
            ggplot(.) +
            geom_bar() +
            aes(x=class) +
            theme_bw()
          
df_final_select_comp_all2 <- df_final_select_comp_all %>%
  dplyr::group_by(CHROM, cluster_start,isotype) %>%
  dplyr::mutate(nclass=ifelse(n_distinct(class)==1, "single", "mixed")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(class == "short_only") %>%
  dplyr::rename(cluster_start_all=cluster_start, cluster_size_all=cluster_size, cluster_end_all=cluster_end) %>%
  join_masks(.) %>%
  dplyr::select(-cluster, -cluster.1) %>%
  dplyr::mutate(edge = ifelse(nclass == "single", NA, ifelse(cluster_start_all==cluster_start | cluster_end_all==cluster_end, T, F)))

df_final_select_comp_all2 %>%
  ggplot(.) +
  geom_bar() +
  aes(x=factor(nclass, levels = c("single","mixed"), labels = c("in_short_only_cluster", "in_mixed_cluster")), fill=edge) +
  theme_bw() +
  xlab(" ")


## longread only divergent bins

df_final_select_comp_all3 <- df_final_select_comp_all %>%
  dplyr::group_by(CHROM, cluster_start, isotype) %>%
  dplyr::mutate(nclass=ifelse(n_distinct(class)==1, "single", "mixed")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(class == "long_only") %>%
  dplyr::rename(cluster_start_all=cluster_start, cluster_size_all=cluster_size, cluster_end_all=cluster_end) %>%
  join_masks(.) %>%
  dplyr::select(-cluster, -cluster.1) %>%
  dplyr::mutate(edge = ifelse(nclass == "single", NA, ifelse(cluster_start_all==cluster_start | cluster_end_all==cluster_end, T, F)))

df_final_select_comp_all3 %>%
  ggplot(.) +
  geom_bar() +
  aes(x=factor(nclass, levels = c("single","mixed"), labels = c("in_long_only_cluster", "in_mixed_cluster")), fill=edge) +
  theme_bw() +
  xlab(" ")

