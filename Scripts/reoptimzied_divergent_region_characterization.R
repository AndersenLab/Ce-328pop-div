library(dplyr)
library(tidyr)
library(glue)
library(valr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

###################### Divergent region ###################


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

#options(scipen=999)

masked_directory <- "Processed_Data/Processed_Masks_nohet_DELINS20kb/"

mask_files <- list.files(masked_directory)

n_sm <- length(mask_files)

# mask frequency thresholds
mask_frq_thresh <- c(0.01, 0.05)

# Extract all regions, select variant counts and coverage
for(sm in 1:length(mask_files)){
  temp_sm_all <- data.table::fread(glue::glue("{masked_directory}{mask_files[sm]}")) %>%
    dplyr::select(CHROM, START_BIN, END_BIN, STRAIN, COUNT, COVERAGE, fraction_SV_bases, window_mask)
  
  if(!exists("population_all")){
    population_all <- temp_sm_all
  } else {
    population_all <- dplyr::bind_rows(population_all, temp_sm_all)
  }
}

#write.table(population_all, file = "Processed_Data/divergent_output_all_strains_all_bins.tsv")
population_all <- read.table(file = "Processed_Data/divergent_output_all_strains_all_bins.tsv")

## Coverage

df_coverage <- population_all %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(mean_cov_genome = mean(COVERAGE)) %>%
  dplyr::ungroup()

df_coverage_fraction <- population_all %>%
  dplyr::left_join(., df_coverage, by = "STRAIN") %>%
  dplyr::mutate(fraction_cov = COVERAGE/mean_cov_genome)

## long-read based optimized paramters (2020/05/11)

cov=0.35
ct=15
cluster_threshold=9000
gap_cluster_threshold=9000
gap=5000

population_all_div <- df_coverage_fraction  %>%
  dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
  dplyr::mutate(div_count = ifelse(COUNT>ct, "Count", "Pass"), div_lowcov = ifelse(fraction_cov<cov, "LowCoverage", "Pass")) %>%
  dplyr::mutate(div_class = ifelse(paste(div_count, div_lowcov, sep="_") == "Pass_Pass", "Pass", "Divergent")) %>%
  dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
  dplyr::mutate(div_id = paste(ifelse(div_count=="Pass","",div_count), ifelse(div_lowcov=="Pass","",div_lowcov), ifelse(div_twoflank=="Pass","",div_twoflank), sep=""), div_class = ifelse(paste(div_count, div_lowcov, div_twoflank, sep="_") == "Pass_Pass_Pass", "Pass", "Divergent"))

###### filtering bins in >=xkb cluster, join xkb gaps ########

### joining ###

cluster_start <- NA
cluster_end <- NA

df_div_cluster <- population_all_div %>%
  dplyr::group_by(STRAIN, CHROM) %>%
  dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
  dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
  dplyr::ungroup() %>%
  dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)

df_div_cluster_filtered <- df_div_cluster %>%
  dplyr::filter(cluster == T)

cluster_start_vec <- df_div_cluster_filtered$cluster_start
start_bin_vec <- df_div_cluster_filtered$START_BIN

for (i in 1:length(cluster_start_vec)) {
  cluster_start_vec[i] <- ifelse(start_bin_vec[i]==0, 0, 
                                 ifelse(!is.na(cluster_start_vec[i]), cluster_start_vec[i], cluster_start_vec[i-1]))

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

## filtering and gap filling

divergent_gap <- function(df=df_joined_masks_freq, threshold = gap_cluster_threshold, gap_threshold=gap) {
  
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

df_gap_bind_list <- divergent_gap(df=df_div_cluster_size, threshold = gap_cluster_threshold, gap_threshold=gap)

df_gap_bind <- df_gap_bind_list[[1]]
df_gap_bind_all_bins <- df_gap_bind_list[[2]] %>%
  dplyr::mutate(div_id = ifelse(div_class == "Divergent", div_id, 
                                ifelse(dplyr::lag(END_BIN, 1) == dplyr::lag(cluster_end, 1) & dplyr::lag(gap_join, 1) == 'join' | 
                                         dplyr::lag(END_BIN, 2) == dplyr::lag(cluster_end, 2) & dplyr::lag(gap_join, 2) == "join" | 
                                         dplyr::lag(END_BIN, 3) == dplyr::lag(cluster_end, 3) & dplyr::lag(gap_join, 3) == "join" |
                                         dplyr::lag(END_BIN, 4) == dplyr::lag(cluster_end, 4) & dplyr::lag(gap_join, 4) == "join", "gap", div_id))) %>%
  dplyr::mutate(div_class = ifelse(div_id == "gap", "Divergent", div_class)) %>%
  dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id)

### Re-clustring after gap joing (to resovle divergent regions with multiple 3kb gaps)

df_div_cluster_gapjoined <- df_gap_bind_all_bins %>%
  dplyr::group_by(STRAIN, CHROM) %>%
  dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
  dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
  dplyr::ungroup() %>%
  dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)

df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined %>%
  dplyr::filter(cluster == T)

cluster_start_gapjoined_vec <- df_div_cluster_gapjoined_filtered$cluster_start
start_bin_gapjoined_vec <- df_div_cluster_gapjoined_filtered$START_BIN

for (i in 1:length(cluster_start_gapjoined_vec)) {
  cluster_start_gapjoined_vec[i] <- ifelse(start_bin_gapjoined_vec[i]==0, 0, 
                                           ifelse(!is.na(cluster_start_gapjoined_vec[i]), cluster_start_gapjoined_vec[i], cluster_start_gapjoined_vec[i-1]))
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

df_divergent_final <- df_div_cluster_gapjoined_size %>%
  dplyr::filter(cluster_size >= cluster_threshold) %>%
  dplyr::group_by(window_ID) %>%
  dplyr::mutate(freq=n()/328, mwf = ifelse(freq > 0.5, 1-freq, freq)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq_bin=ifelse(freq < 0.01, "Rare", ifelse(freq < 0.05, "Intermediate", "Common"))) %>%
  dplyr::mutate(is_count = ifelse(div_id %in% c("Count", "CountLowCoverage"), 1, 0)) %>%
  dplyr::group_by(STRAIN, CHROM, cluster_start, cluster_end) %>%
  dplyr::mutate(flag_del = ifelse(max(is_count)==0, "Del", "Pass")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(flag_del == "Pass")

write.table(df_divergent_final, file = "Processed_Data/df_divergent_final.tsv")

df_divergent_final <- read.table(file = "Processed_Data/df_divergent_final.tsv")

### summary data for each isotype ###

df_divergent_final_isotype <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end)

write.csv(df_divergent_final_isotype, file="Manuscript/Supplementary_Data/DataS3_divergent_regions_isotypes.csv")

## total size

nrow(dplyr::distinct(df_divergent_final, window_ID))/1e3

nrow(dplyr::distinct(df_divergent_final, window_ID))/sum(df_chr_length$stop)*1e3

df_divergent_final %>%
  dplyr::distinct(window_ID, mwf,freq_bin) %>%
  dplyr::group_by(freq_bin) %>%
  dplyr::summarise(size=n()/1e3)

## divergent vs non-divergent stats

df_div_stats <- df_coverage_fraction %>%
  dplyr::select(CHROM, START_BIN, END_BIN, STRAIN, COUNT, fraction_cov) %>%
  dplyr::left_join(., dplyr::select(df_divergent_final,CHROM, START_BIN,END_BIN,STRAIN,div_class), by=c('CHROM','START_BIN','END_BIN','STRAIN')) %>%
  dplyr::mutate(div_class = ifelse(is.na(div_class), "Pass", "Divergent"))

save(join_masks, df_div_cluster_gapjoined_size, df_divergent_final, df_div_stats, file="Processed_Data/divergent_classification.RData")
#load("Processed_Data/divergent_classification.RData")

df_div_stats %>%
  dplyr::group_by(div_class) %>%
  dplyr::summarise(mean_count = mean(COUNT), mean_cov = mean(fraction_cov)) %>%
  dplyr::ungroup()

df_div_stats_iso <- df_div_stats %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(mean_count = mean(COUNT), mean_cov = mean(fraction_cov)) %>%
  dplyr::ungroup()

mean(df_div_stats_iso$mean_count)

df_div_stats_iso_div <- df_div_stats %>%
  dplyr::filter(div_class == "Divergent") %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(mean_count = mean(COUNT), mean_cov = mean(fraction_cov)) %>%
  dplyr::ungroup()
  
mean(df_div_stats_iso_div$mean_count)
mean(df_div_stats_iso_div$mean_cov)

## total size of example strains

df_divergent_final %>%
  dplyr::filter(STRAIN == "CB4856") %>%
  dplyr::summarise(size=n()/1e3)

df_divergent_final %>%
  dplyr::filter(STRAIN == "XZ1516") %>%
  dplyr::summarise(size=n()/1e3)

### total divergent region size with different threshold

df_div_common <- df_divergent_final %>%
  dplyr::filter(freq_bin == "Common") %>%
  dplyr::distinct(chr=CHROM, start=START_BIN, stop = END_BIN) %>%
  dplyr::arrange(., chr, start)

readr::write_tsv(df_div_common, "Processed_Data/divergent_regions/Common_divergent_regions.tsv")

df_div_intermediate <- df_divergent_final %>%
  dplyr::filter(freq_bin == "Intermediate") %>%
  dplyr::distinct(chr=CHROM, start=START_BIN, stop = END_BIN) %>%
  dplyr::arrange(., chr, start)

readr::write_tsv(df_div_intermediate, "Processed_Data/divergent_regions/Intermediate_divergent_regions.tsv")

df_div_rare <- df_divergent_final %>%
  dplyr::filter(freq_bin == "Rare") %>%
  dplyr::distinct(chr=CHROM, start=START_BIN, stop = END_BIN) %>%
  dplyr::arrange(., chr, start)

readr::write_tsv(df_div_rare, "Processed_Data/divergent_regions/Rare_divergent_regions.tsv")

df_div_all <- df_divergent_final %>%
  dplyr::distinct(chr=CHROM, start=START_BIN, stop = END_BIN) %>%
  dplyr::arrange(., chr, start)

readr::write_tsv(df_div_all, "Processed_Data/divergent_regions/All_divergent_regions.tsv")

## clustered

df_div_common_clustered <- df_divergent_final %>%
  dplyr::filter(freq_bin == "Common") %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN) %>%
  dplyr::arrange(., CHROM, START_BIN) %>%
  join_masks(.) %>%
  dplyr::distinct(chr=CHROM, start=cluster_start, stop = cluster_end) %>%
  dplyr::arrange(., chr, start) %>%
  na.omit()

readr::write_tsv(df_div_common_clustered, "Processed_Data/divergent_regions/Common_divergent_regions_clustered.tsv")

df_div_intermediate_clustered <- df_divergent_final %>%
  dplyr::filter(freq_bin == "Intermediate") %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN) %>%
  dplyr::arrange(., CHROM, START_BIN) %>%
  join_masks(.) %>%
  dplyr::distinct(chr=CHROM, start=cluster_start, stop = cluster_end) %>%
  dplyr::arrange(., chr, start) %>%
  na.omit()

readr::write_tsv(df_div_intermediate_clustered, "Processed_Data/divergent_regions/Intermediate_divergent_regions_clustered.tsv")

df_div_rare_clustered <- df_divergent_final %>%
  dplyr::filter(freq_bin == "Rare") %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN) %>%
  dplyr::arrange(., CHROM, START_BIN) %>%
  join_masks(.) %>%
  dplyr::distinct(chr=CHROM, start=cluster_start, stop = cluster_end) %>%
  dplyr::arrange(., chr, start) %>%
  na.omit()

readr::write_tsv(df_div_rare_clustered, "Processed_Data/divergent_regions/Rare_divergent_regions_clustered.tsv")

df_div_all_clustered <- df_divergent_final %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN) %>%
  dplyr::arrange(., CHROM, START_BIN) %>%
  join_masks(.) %>%
  dplyr::distinct(chrom=CHROM, start=cluster_start, end = cluster_end) %>%
  dplyr::arrange(., chrom, start) %>%
  na.omit()

readr::write_tsv(df_div_all_clustered, "Processed_Data/divergent_regions/All_divergent_regions_clustered.tsv")
#readr::write_tsv(df_div_all_clustered, "Processed_Data/divergent_regions/All_divergent_regions_clustered.bed", col_names=F)

#df_div_all_clustered <- readr::read_tsv(file="Processed_Data/divergent_regions/All_divergent_regions_clustered.tsv")

WS270_gene <- read.table("Processed_Data/gene.protein_coding_WS270.bed") %>%
  dplyr::select(V1,V2,V3,V4) %>%
  tidyr::separate(V4, into=c("V4","V5")) %>%
  dplyr::select(chrom=V1, start=V2, end=V3, gene=V5)

### Genes in species-wide divergent regions

Divergent_genes <- valr::bed_intersect(df_div_all_clustered, WS270_gene)

Divergent_genes_only <- unique(Divergent_genes$gene.y)

write.table(Divergent_genes_only, "Processed_Data/divergent_regions/genes/Divergent_genes_only.tsv", row.names=F, col.names=F, quote=F)

### Genes in divergent regions for each strains

Divergent_genes_strains <- NA

for (st in unique(df_divergent_final$STRAIN)) {
  
  df_strain <- df_divergent_final %>%
    dplyr::filter(STRAIN == st) %>%
    dplyr::distinct(CHROM, START_BIN, END_BIN) %>%
    dplyr::arrange(., CHROM, START_BIN) %>%
    join_masks(.) %>%
    dplyr::distinct(chrom=CHROM, start=cluster_start, end = cluster_end) %>%
    dplyr::arrange(., chrom, start) %>%
    na.omit()
    
  genes_temp <- valr::bed_intersect(df_strain, WS270_gene) %>%
    dplyr::mutate(strain = st) %>%
    dplyr::distinct(gene.y, strain) %>%
    dplyr::rename(gene=gene.y)
  
  Divergent_genes_strains <- rbind(Divergent_genes_strains, genes_temp)
  
}

write.table(na.omit(Divergent_genes_strains), "Processed_Data/divergent_regions/genes/Divergent_genes_strains.tsv", row.names=F, col.names=F, quote=F)


### Divergent genes at least with two strain calls

Divergent_genes_filtered <- Divergent_genes_strains %>%
  na.omit() %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n_strain = n_distinct(strain)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_strain > 1)

write.table(Divergent_genes_filtered$gene, "Processed_Data/divergent_regions/genes/Divergent_genes_only_private_filtered.tsv", row.names=F, col.names=F, quote=F)


### Genes in non-divergent arms

df_non_div_arms <- readr::read_tsv(file="Processed_Data/divergent_regions/df_arm_nondiv_coords.tsv") %>%
  dplyr::rename(chrom=CHROM, start=START_BIN, end=END_BIN)

Non_div_arm_genes <- valr::bed_intersect(df_non_div_arms, WS270_gene)

Non_div_arm_genes_only <- unique(Non_div_arm_genes$gene.y)

write.table(Non_div_arm_genes_only, "Processed_Data/divergent_regions/genes/Non_div_arm_genes_only.tsv", row.names=F, col.names=F, quote=F)

