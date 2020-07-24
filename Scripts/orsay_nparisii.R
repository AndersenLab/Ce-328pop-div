library(dplyr)
library(tidyr)
library(ggplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/divergent_classification.RData")

# parse gff
wbgff = ape::read.gff("Processed_Data/c_elegans.PRJNA13758.WS270.annotations.gff3")

wbgff_gene_info <- wbgff %>% 
  dplyr::filter(type == "gene") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(parse_atr = strsplit(attributes, split = ";")) %>% 
  dplyr::mutate(transcript_name = ifelse(length(grep("sequence_name=", parse_atr, value = T))>0, 
                                         gsub("sequence_name=", "", grep("sequence_name=", parse_atr, value = T)), NA),
                locus = ifelse(length(grep("locus=", parse_atr, value = T))>0, 
                               gsub("locus=", "", grep("locus=", parse_atr, value = T)), NA),
                wbgene = ifelse(length(grep("sequence_name=", parse_atr, value = T))>0, 
                                gsub("Name=", "", grep("Name=", parse_atr, value = T)), NA),
                map_pos = ifelse(length(grep("interpolated_map_position=", parse_atr, value = T))>0, 
                                 gsub("interpolated_map_position=", "", grep("interpolated_map_position=", parse_atr, value = T)), NA)) 


# load divergent
final_div_regions <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  dplyr::select(chrom = CHROM, start = cluster_start, end = cluster_end, isotype = STRAIN) %>%
  dplyr::group_by(chrom,   start,     end) %>%
  dplyr::mutate(divergent_strains = paste(isotype, collapse = ":")) %>%
  dplyr::distinct(chrom,   start,     end, .keep_all = T) %>%
  dplyr::select(chrom, start, end, divergent_strains)

final_gene_info  <- wbgff_gene_info %>%
  dplyr::rename(chrom = seqid) %>%
  valr::bed_intersect(wbgff_gene_info, final_div_regions) %>%
  dplyr::select(wbgene = wbgene.x, divergent_strains = divergent_strains.y, .overlap) %>%
  dplyr::group_by(wbgene) %>%
  dplyr::mutate(divergent_strains = paste(divergent_strains, collapse = ";")) %>%
  dplyr::distinct(wbgene, .keep_all=T) %>%
  dplyr::left_join(wbgff_gene_info,., by = "wbgene") %>%
  dplyr::mutate(is_divergent = ifelse(is.na(divergent_strains), F,T)) %>%
  dplyr::select(chrom = seqid, start, end, transcript_name:is_divergent)

#orsay
orsay <- data.table::fread("Processed_data/Chen2017_osray_DE_N2.csv") %>%
  na.omit()

gene_haps_orsay <- orsay %>%
  dplyr::rename(transcript_name = 'Sequence Name') %>%
  dplyr::left_join(., final_gene_info,by = "transcript_name") %>%
  dplyr::select(wbgene = `WormBase Gene ID`, log2_fold_change = `Fold-change (log2)`, is_divergent, FDR) %>%
  dplyr::mutate(pathogen="Orsay virus")

#Nparisii
npari <- data.table::fread("Processed_data/Chen2017_Nparisii_DE_N2.csv") %>%
  na.omit()

gene_haps_npari <- npari %>%
  dplyr::rename(transcript_name = 'Sequence Name') %>%
  dplyr::left_join(., final_gene_info,by = "transcript_name") %>%
  dplyr::select(wbgene = `WormBase Gene ID`, log2_fold_change = `Fold-change (log2)`, is_divergent, FDR) %>%
  dplyr::mutate(pathogen="N. parisii")

patho_genes <- rbind(gene_haps_orsay, gene_haps_npari)

patho_genes %>%
  dplyr::group_by(pathogen) %>%
  #dplyr::group_by(pathogen, is_divergent) %>%
  dplyr::summarise(n=n())

131/196 #nparisii
85/130 #orsay

plot_patho1 <- patho_genes %>%
  ggplot(.) +
  geom_jitter(width=0.3, alpha= 0.7) +
  geom_boxplot(outlier.alpha = 0, alpha=0.7) +
  aes(x=is_divergent, y=log2_fold_change, fill=is_divergent) +
  theme_bw() +
  scale_fill_brewer(palette = 'Set1') +
  theme(axis.text = element_text(size=10, color = "black"),
        axis.title = element_text(size=11, color = "black"),
        panel.grid = element_blank(),
        legend.position = 'none') +
  facet_grid(~pathogen) +
  labs(x="Hyper-divergent", y="log2(fold change)")

plot_patho2 <- patho_genes %>%
  ggplot(.) +
  geom_jitter(width=0.3, alpha= 0.7) +
  geom_boxplot(outlier.alpha = 0, alpha=0.7) +
  aes(x=is_divergent, y=-log10(FDR), fill=is_divergent) +
  theme_bw() +
  scale_fill_brewer(palette = 'Set1') +
  theme(axis.text = element_text(size=10, color = "black"),
        axis.title = element_text(size=11, color = "black"),
        panel.grid = element_blank(),
        legend.position = 'none') +
  facet_grid(~pathogen) +
  labs(x="Hyper-divergent", y="-log10(FDR)")

plot_patho3 <- patho_genes %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.6, size =2) +
  aes(x=log2_fold_change, y=-log10(FDR), fill=is_divergent) +
  theme_bw() +
  scale_fill_brewer(palette = 'Set1') +
  theme(axis.text = element_text(size=10, color = "black"),
        axis.title = element_text(size=11, color = "black"),
        panel.grid = element_blank(),
        legend.position = c(0.7,0.8)) +
  facet_grid(~pathogen) +
  labs(x="log2(fold change)", y="-log10(FDR)", fill = "Hyper-divergent")

plot_patho <- cowplot::plot_grid(plot_patho1, plot_patho2, plot_patho3,
                              labels = c("a","b","c"), ncol=1,
                              label_size = 12, align= 'hv')
plot_patho

save(patho_genes, file="Processed_Data/pathogen_exp.RData")
