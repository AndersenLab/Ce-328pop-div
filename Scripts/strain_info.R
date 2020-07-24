library(dplyr)
library(tidyr)
library(geosphere)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load(file="Raw/WI_328_list.RData")

################ Figure 1. Geographic structure of global C. elegans populations ####################

### Geographic grouping

df_CeNDR <- read.csv(file = "Raw/20190710_CeNDR_strain_data.csv", fill = TRUE, header=TRUE, sep = ",")

isotype_CeNDR <- as.character(unique(df_CeNDR$isotype))
isotype_CeNDR[!isotype_CeNDR %in% strain_vector]

df_CeNDR$isotype <- gsub("CB4851", "ECA243", df_CeNDR$isotype)
df_CeNDR$isotype <- gsub("CB4853", "ECA246", df_CeNDR$isotype)
df_CeNDR$isotype <- gsub("CB4855", "ECA248", df_CeNDR$isotype)
df_CeNDR$isotype <- gsub("CB4857", "ECA250", df_CeNDR$isotype)
df_CeNDR$isotype <- gsub("CB4858", "ECA251", df_CeNDR$isotype)
df_CeNDR$isotype <- gsub("PB306", "ECA259", df_CeNDR$isotype)

df_seq_strain <- read.table("Raw/strain_input.tsv") %>%
  dplyr::distinct(strain=V1) %>%
  dplyr::left_join(., df_CeNDR, by = 'strain') %>%
  dplyr::filter(!is.na(reference_strain)) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(n_strain = n()) %>%
  dplyr::ungroup()

write.csv(df_seq_strain, "Manuscript/Supplementary_table/TableS1_609_wild_strains.csv")

unique(df_seq_strain$isotype) %in% strain_vector

input_strain <- as.character(df_seq_strain$strain) ### 609 whole-genome sequenced strains / 328 isotypes

nrow(distinct(filter(df_seq_strain, n_strain > 1), isotype)) ### 87 isotypes have more than 1 strain

# overlap_geo

df_seq_strain_multi <- df_seq_strain %>%
  dplyr::filter(n_strain > 1) %>%
  dplyr::distinct(isotype, country, state, .keep_all=T) %>%
  dplyr::select(strain, isotype, country, state)

strain_info <- df_seq_strain %>%
  dplyr::filter(isotype %in% strain_vector) %>%
  dplyr::rename(long = longitude, lat = latitude) %>%
  dplyr::filter(lat != "None")%>%
  dplyr::filter(lat != "") %>%
  dplyr::filter(!is.na(lat))

strain_info$lat <- as.numeric(as.character(strain_info$lat))
strain_info$long <- as.numeric(as.character(strain_info$long))

length(unique(strain_info$isotype))

## 605 strains / 324 isotypes have sampling location info | 4 strains/isotypes without long/lat

strain_info_nogeo <- df_seq_strain %>%
  dplyr::rename(long = longitude, lat = latitude) %>%
  dplyr::filter(is.na(lat) | lat == "None" | lat == "")

strain_info_nogeo$lat <- as.numeric(as.character(strain_info_nogeo$lat))
strain_info_nogeo$long <- as.numeric(as.character(strain_info_nogeo$long))

## calculate distance (by lat and long) between any strains to identify independent isolation (>50km) of isotypes

df_dist_all <- NULL
df_prox_10 <- NULL
df_dist_iso <- NULL
df_prox_10_iso <- NULL

for (i in 1:nrow(strain_info)) {
  
  strain_comp <- strain_info$strain[i]
  isotype_comp <- strain_info$isotype[i]
  long_comp <- strain_info$long[i]
  lat_comp <- strain_info$lat[i]
  
  df_dist <- strain_info %>%
    dplyr::mutate(long_comp=long_comp, lat_comp=lat_comp, strain_comp=strain_comp, isotype_comp=isotype_comp) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(dist = geosphere::distm(cbind(long, lat), cbind(long_comp, lat_comp), fun=distGeo)/1e3) %>%
    dplyr::select(strain, strain_comp, isotype,isotype_comp, long, long_comp, lat, lat_comp, dist)
  
  df_dist_all <- rbind(df_dist_all, df_dist)
  
  df1 <- df_dist %>%
    dplyr::filter(dist < 1e1, strain != strain_comp) %>%
    dplyr::filter(isotype != isotype_comp)
  
  df_prox_10 <- rbind(df_prox_10, df1)
  
  df2 <- df_dist %>%
    dplyr::filter(isotype == isotype_comp & dist > 5e1) 
  
  df_dist_iso <- rbind(df_dist_iso, df2)
  
  df3 <- df_dist %>%
    dplyr::filter(dist < 1e1, isotype == isotype_comp)
  
  df_prox_10_iso <- rbind(df_prox_10_iso, df3)
  
}

length(unique(df_dist_iso$strain))
length(unique(df_dist_iso$isotype))

length(unique(df_prox_10$strain))
length(unique(df_prox_10$isotype))

## distribuition of  in-isotype and between-isotype pairwise distances

df_dist_iniso <- df_dist_all %>%
  dplyr::filter(isotype_comp == isotype) %>%
  dplyr::mutate(class = 'in_isotype')

df_dist_bwiso <- df_dist_all %>%
  dplyr::filter(strain == isotype, strain_comp == isotype_comp) %>%
  dplyr::filter(!isotype_comp == isotype) %>%
  dplyr::mutate(class = 'between_isotype')

df_dist_in_bw_iso <- rbind(df_dist_iniso, df_dist_bwiso)

df_dist_in_bw_iso %>%
  ggplot(.) +
  geom_histogram(binwidth=200, origin = 0, position="identity", alpha = 0.6) +
  aes(x=dist, fill=class) +
  scale_fill_manual(values=c("#377EB8", "#E41A1C")) +
  theme_bw() +
  labs(x="Distance (km)")

df_dist_in_bw_iso %>%
  ggplot(.) +
  geom_histogram(binwidth=5, origin = 0, position="identity", alpha = 0.6) +
  aes(x=dist, fill=class) +
  scale_fill_manual(values=c("#377EB8", "#E41A1C")) +
  theme_bw() +
  xlim(0,200) +
  labs(x="Distance (km)")

df_dist_in_bw_iso %>%
  ggplot(.) +
  geom_histogram(binwidth=2, origin = 0, position="identity", alpha = 0.6) +
  aes(x=dist, fill=class) +
  scale_fill_manual(values=c("#377EB8", "#E41A1C")) +
  theme_bw() +
  xlim(0,100) +
  labs(x="Distance (km)")

## group local strains with same isotypes -- 50km as a distance threshold for same local strains

df_dist_indep <- df_dist_all %>%
  dplyr::filter(isotype %in% unique(df_dist_iso$isotype)) %>%
  dplyr::filter(isotype == isotype_comp & strain != strain_comp) %>%
  dplyr::filter(!strain == isotype)

for (st in unique(df_dist_indep$isotype)) {
  
  df1 <- df_dist_indep  %>%
    dplyr::filter(strain_comp == st) %>%
    dplyr::filter(dist < 50)
  
  rm_strains <- unique(df1$strain)
  
  df_dist_indep <- df_dist_indep %>%
    dplyr::filter(!strain %in% rm_strains) %>%
    dplyr::filter(!strain_comp %in% rm_strains)
  
}

for (st in unique(df_dist_indep$strain)) {
  
  df1 <- df_dist_indep  %>%
    dplyr::filter(strain != st & strain_comp == st) %>%
    dplyr::filter(dist < 50)
  
  rm_strains <- unique(df1$strain)
  
  df_dist_indep <- df_dist_indep %>%
    dplyr::filter(!strain %in% rm_strains) %>%
    dplyr::filter(!strain_comp %in% rm_strains)
  
}

indep_strains <- unique(df_dist_indep$strain)
length(unique(df_dist_indep$strain)) ### 328 isotype + 19 additional independent samples = 347 / 343 with available geo info 

indep_strain_info <- strain_info %>%
  dplyr::filter(strain == isotype | strain %in% indep_strains) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(n_isotype = n()) %>%
  dplyr::ungroup() %>%
  dplyr::bind_rows(strain_info_nogeo)

nrow(distinct(filter(indep_strain_info, n_isotype >1), isotype))

write.csv(indep_strain_info, file = "347_indep_strains.csv")
save(df_dist_in_bw_iso, indep_strain_info, file="Processed_Data/indep_strain_info.RData")
