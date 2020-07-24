library(tidyverse)

tsv_files <- fs::dir_ls("./", regexp = "\\.tsv$")
interspec_aa <- tsv_files %>% 
  map_dfr(data.table::fread) %>%
  dplyr::rename(ref_trans = V1, comp_trans = V2, intra_aa = V3, brig_trans = V4, nig_trans = V5, inter_aa = V6, divergent_status = V7, comp_strain = V8)

interspec_aa %>%
  dplyr::group_by(ref_trans, brig_trans, nig_trans, divergent_status) %>%
  dplyr::summarise(mean_intra = mean(intra_aa, na.rm = T),
                   mean_inter = mean(inter_aa, na.rm = T)) %>%
  tidyr::gather(stra_comp, aa_id, -(ref_trans:divergent_status)) %>%
  dplyr::mutate(stra_comp = factor(stra_comp, levels = c("mean_intra","mean_inter"))) %>%
  ggplot()+
  aes(x = stra_comp, y = aa_id, fill = divergent_status)+
  geom_boxplot(outlier.colour = NA) + scale_x_discrete(labels=c("mean_intra" = "Comparison between\n16 C. elegans isotypes",
                                                                "mean_inter" = "Comparison between\nC. nigoni and C. briggsae")) +
  ylim(90, 100) +
  theme_bw(18) +
  theme(axis.title.x = element_blank())+
  labs(y = "Amino acid identity (%)")