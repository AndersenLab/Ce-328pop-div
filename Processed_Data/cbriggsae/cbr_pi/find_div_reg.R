# joining adjacent masks
join_masks <- function(mask_file = NULL, jump_tolerance = 2001){
  common_cluster <- NA
  common_cluster_start <- NA
  
  for (i in c(1:nrow(mask_file)-1)) {
    
    if(i==1){
      
      common_cluster[i] <- ifelse(mask_file$startbin[i] == mask_file$startbin[i+1], 'yes', 'no')
      common_cluster_start[i] <- mask_file$endbin[i]
      
    } else {
      
      common_cluster[i] <- ifelse(mask_file$endbin[i-1] > mask_file$startbin[i] - jump_tolerance | mask_file$endbin[i] > mask_file$startbin[i+1] - jump_tolerance, 'yes', 'no')
      common_cluster_start[i] <- ifelse(mask_file$chrom[i] == mask_file$chrom[i-1] & common_cluster[i] == 'yes', 
                                        ifelse(mask_file$endbin[i-1] > mask_file$startbin[i] - jump_tolerance | mask_file$endbin[i] > mask_file$startbin[i+1] - jump_tolerance, 
                                               ifelse(common_cluster[i-1] == "no", mask_file$startbin[i], 
                                                      ifelse(mask_file$endbin[i-1] == mask_file$startbin[i], common_cluster_start[i-1], mask_file$startbin[i])), mask_file$startbin[i]), mask_file$startbin[i])
      
    }
  }
  
  common_cluster[nrow(mask_file)] <- NA
  common_cluster_start[nrow(mask_file)] <- NA
  
  mask_file_cluster <- data.frame(mask_file, common_cluster, common_cluster_start) %>%
    dplyr::group_by(chrom, common_cluster_start) %>%
    dplyr::mutate(common_cluster_size=ifelse(is.na(common_cluster_start), 1000, n()*1000)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(common_cluster_end=common_cluster_start+common_cluster_size)
  
}


for(strains in list.files("variant_counts/")){
  
  strain_name <- str_split(strains, pattern = "\\.")[[1]][1]
  
  cts <- data.table::fread(glue::glue("variant_counts/{strains}")) %>%
    dplyr::mutate(strain = strain_name) %>%
    dplyr::rename(chrom=V1, startbin = V2, endbin = V3, vct = V4)
  
  # 19 comes from the 99% of variants in bin across the population
  if(nrow(cts %>% dplyr::filter(vct >= 19)) > 0){
    temp_join <- join_masks(cts %>% dplyr::filter(vct >= 19))
  } else {
    print(glue::glue("{strain_name} has no outlier regions"))
  }
  
  
  if(!exists("cbr_variant_counts")){
    cbr_variant_counts <- cts
    cbr_div_counts <- temp_join
  } else {
    cbr_variant_counts <- dplyr::bind_rows(cbr_variant_counts, cts)
    cbr_div_counts <- dplyr::bind_rows(cbr_div_counts, temp_join)
  }
  
}

div_bins <- cbr_variant_counts %>%
  dplyr::filter(vct>=19) %>%
  dplyr::arrange(chrom, startbin,endbin, strain) %>%
  na.omit() %>%
  dplyr::select( chrom, start = startbin, end = endbin, strain)  

# 2)
div_bins_merged <- valr::bed_merge(dplyr::group_by(div_bins, strain),max_dist = 5001 ) %>%
  dplyr::mutate(size = end-start)%>%
  dplyr::mutate(large_bin = ifelse(size >= 9e3, "large", "small"))


strain_plot_id <- div_bins_merged %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(sum = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-sum) %>%
  dplyr::mutate(ID = 1:n())

div_bins_merged %>%
  dplyr::left_join(.,strain_plot_id,by = "strain") %>%
  ggplot(.) +
  geom_rect(aes(xmin =  start/1e6, xmax = end/1e6, ymin = ID-0.5 , ymax = ID+0.5), color = 'black', size = 0.5) +
  geom_rect(aes(xmin =  start/1e6, xmax = end/1e6, ymin = ID-0.5 , ymax = ID+0.5), color = 'red', fill="red", size = 0.5, data = div_bins_merged %>% dplyr::filter(large_bin == "large") %>% dplyr::left_join(.,strain_plot_id,by = "strain") ) +
  theme_bw() +
  facet_grid(~chrom, scales = 'free') +
  theme(axis.text.y = element_blank(), 
        text = element_text(size=12, color='black'), 
        panel.grid = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = 'none',
        axis.title.x=element_text(size=12, color = 'black'), 
        axis.text.x=element_text(size=11, color = 'black')) +
  labs(x="Genomic position (Mb)") +
  scale_y_continuous(expand = c(0.01, 0.01))
  

test_df <- cbr_variant_counts %>%
  dplyr::filter(vct>=quantile(cbr_variant_counts$vct, probs = 0.99)) %>%
  dplyr::filter(strain == "JU1341")


quantile(cbr_variant_counts$vct, probs = 0.01)

cbr_variant_counts %>%
  dplyr::filter(vct>quantile(cbr_variant_counts$vct, probs = 0.99)) %>%
  ggplot()+
  aes(x = vct) +
  geom_histogram(size = 0.25) +
  geom_vline(aes(xintercept = quantile(cbr_variant_counts$vct, probs = 0.99)))
