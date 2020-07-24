library(tidyverse)
library(fs)
library(valr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# parse gff
wbgff = ape::read.gff("data/c_elegans.PRJNA13758.WS270.annotations.gff3")

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
                                 gsub("interpolated_map_position=", "", grep("interpolated_map_position=", parse_atr, value = T)), NA)) %>%
  dplyr::filter(grepl("protein_coding", attributes))

# genomic 1kb bins
n_1kb_bins <- data.table::fread("data/df_bins.tsv") %>%
  dplyr::rename(chrom=chr,end=stop)

# # load mapping results results - takes a while
# csv_files <- fs::dir_ls("data/emma_results/", regexp = "\\.tsv$")
# 
# csv_files <- csv_files[csv_files!=c("data/emma_results/QTL_peaks.tsv")]
# 
# #64,053 snps for BF
# 
# td_df <- csv_files %>% 
#   map_dfr(read_tsv) %>%
#   na.omit() 
# 
# write.table(td_df, file = "data/concatenated_emma_results.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

td_df <- data.table::fread("data/concatenated_emma_results.tsv")

sig_burden <- td_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(b_t = trait) %>%
  dplyr::mutate(bacteria = strsplit(trait, split = "_")[[1]][1],
                trait = strsplit(trait, split = "_")[[1]][2]) %>% 
  dplyr::filter(!grepl("green|cv|iqr", trait)) %>% 
  dplyr::filter(grepl("JUb", bacteria)) %>%
  dplyr::arrange(desc(log10p)) %>%
  dplyr::distinct(bacteria, peakPOS, .keep_all = T) %>%
  dplyr::select(chrom = CHROM, start = startPOS, end = endPOS, log10p, bacteria, trait, b_t) 

phenotyped_strains <- unique(td_df$strain)

# load divergent
final_div_regions <- data.table::fread("data/DataS3_divergent_regions_isotypes.csv") %>%
  dplyr::filter(STRAIN %in% phenotyped_strains) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::select(chrom = CHROM, start = cluster_start, end = cluster_end, isotype = STRAIN) %>%
  valr::bed_merge() %>%
  dplyr::group_by(chrom,   start,     end) %>%
  dplyr::mutate(divergent_strains = paste(isotype, collapse = ":")) %>%
  dplyr::distinct(chrom,   start,     end, .keep_all = T) %>%
  dplyr::select(chrom, start, end, divergent_strains) %>%
  dplyr::rowwise() %>%
  dplyr::filter(end-start>=9000)

d_regions <- valr::bed_merge(final_div_regions %>% dplyr::ungroup()) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(div_marker = paste(chrom,start,end, collapse = "_")) %>%
  valr::bed_intersect(.,n_1kb_bins) %>%
  dplyr::filter(.overlap!=0) %>%
  dplyr::select(chrom, start = start.y, end = end.y, div_marker = div_marker.x)

fraction_d_region <- valr::bed_intersect(d_regions,final_div_regions %>% dplyr::ungroup()) %>%
  dplyr::filter(.overlap!=0) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(frac_strains = length((strsplit(divergent_strains.y, split = ":")[[1]]))/length(phenotyped_strains)) %>%
  dplyr::filter(frac_strains > 0.05) %>%
  dplyr::distinct(chrom, start.x, end.x, div_marker.x) %>%
  dplyr::select(chrom, start = start.x, end = end.x, div_marker = div_marker.x)

sig_burden_1kb <- sig_burden %>%
  dplyr::filter(grepl("JUb", bacteria)) %>%
  valr::bed_intersect(.,n_1kb_bins) %>%
  dplyr::filter(.overlap!=0) %>%
  dplyr::distinct(chrom, start.y, end.y, .keep_all = T) %>%
  dplyr::select(chrom, start = start.y, end = end.y, log10p.x:b_t.x)

div_intersect_qtl <- valr::bed_intersect(sig_burden_1kb,fraction_d_region) %>%
  dplyr::filter(.overlap!=0)

div_not_qtl <- valr::bed_intersect(sig_burden_1kb, fraction_d_region, invert = T) %>%
  dplyr::distinct(chrom, start, end, .keep_all = T)

# over representation:
phyper(nrow(div_intersect_qtl), 
       nrow(fraction_d_region), 
       nrow(n_1kb_bins)-nrow(fraction_d_region), 
       nrow(div_not_qtl)+nrow(div_intersect_qtl), 
       lower.tail= FALSE)

length(unique(div_intersect_qtl$div_marker.y))
length(unique(fraction_d_region$div_marker))

# [1] 2.588783e-38

ggplot()+
  geom_rect(data = fraction_d_region %>% dplyr::ungroup() %>% dplyr::select(div_marker) %>% tidyr::separate(div_marker, into = c("chrom", "start", "end"), convert = T) %>% dplyr::distinct(),
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = 0 , ymax = 30), size = 0.2, alpha = 1, color = "gray70") + 
  geom_rect(data = sig_burden %>% dplyr::filter(grepl("JUb", bacteria)) %>% dplyr::ungroup() %>% dplyr::mutate(nplot = 1:n()),
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = nplot-0.5 , ymax = nplot+0.5, fill = bacteria), size = 0.2, alpha = 1, color = "black")+
  geom_rect(data = n_1kb_bins,
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = 0 , ymax = 0), size = 0.2, alpha = 1, color = NA, fill = NA)+
  facet_grid(.~chrom, scales = "free", space = "free") +
  theme_bw(18) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top") +
  labs(x = "Genomic Position (Mb)", y = "Trait", fill = "Bacteria")

ggsave("plots/natural_bacteria_qtl_divergent.pdf", height = 6, width = 12)
ggsave("plots/natural_bacteria_qtl_divergent.png", height = 6, width = 12)

ggplot()+
  geom_rect(data = fraction_d_region %>% dplyr::ungroup() %>% dplyr::select(div_marker) %>% tidyr::separate(div_marker, into = c("chrom", "start", "end"), convert = T) %>% dplyr::distinct(),
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = 0 , ymax = 30), size = 0.2, alpha = 1, color = "gray70") + 
  geom_rect(data = sig_burden %>% dplyr::filter(grepl("JUb", bacteria)) %>% dplyr::ungroup() %>% dplyr::mutate(nplot = 1:n()),
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = nplot-0.5 , ymax = nplot+0.5, fill = bacteria), size = 0.2, alpha = 1, color = "black")+
  geom_rect(data = div_intersect_qtl %>% dplyr::select(chrom, start = start.x , end = end.x, b_t.x.x) %>% dplyr::group_by(b_t.x.x) %>% valr::bed_merge(),
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = 30 , ymax = 32), size = 0.2, alpha = 1, color = "black", fill = "red")+
  geom_rect(data = div_not_qtl %>% dplyr::group_by(b_t.x) %>% valr::bed_merge(),
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = 32 , ymax = 34), size = 0.2, alpha = 1, color = "black", fill = "blue")+
  geom_rect(data = n_1kb_bins,
            aes(xmin =  start/1e6, xmax = end/1e6, ymin = 0 , ymax = 0), size = 0.2, alpha = 1, color = NA, fill = NA)+
  facet_grid(.~chrom, scales = "free", space = "free") +
  theme_bw(18) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Trait", fill = "Bacteria")




####### old analysis
# 

# sig_burden <- td_df %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(b_t = trait) %>%
#   dplyr::mutate(bacteria = strsplit(trait, split = "_")[[1]][1],
#                 trait = strsplit(trait, split = "_")[[1]][2]) %>%
#   dplyr::distinct(bacteria, peakPOS, .keep_all = T) %>%
#   dplyr::filter(!grepl("green|cv|iqr", trait)) %>%
#   dplyr::select(chrom = CHROM, start = startPOS, end = endPOS, log10p, bacteria, trait, b_t) 


# 
# 
# sig_burden_1kb <- sig_burden %>%
#   valr::bed_intersect(.,n_1kb_bins) %>%
#   dplyr::filter(.overlap!=0) %>%
#   dplyr::distinct(start.y, end.y, .keep_all = T) %>%
#   dplyr::select(chrom, start = start.y, end = end.y, log10p.x:b_t.x)
# 
# div_intersect_qtl <- valr::bed_intersect(sig_burden_1kb,d_regions) %>%
#   dplyr::filter(.overlap!=0)
# 
# div_not_qtl <- valr::bed_intersect(sig_burden_1kb, d_regions, invert = T) %>%
#   dplyr::distinct(chrom, start, end, .keep_all = T)
# 
# 
# 
# # nrow(div_intersect_qtl) = n 1kb bins that intersect divergent bins
# # nrow(d_regions) = n of 1kb divergent bins
# # nrow(n_1kb_bins)-nrow(d_regions) = n 1kb non-divergent bins in genome (n bins in genome - n divergent genes)
# # nrow(div_not_qtl) +nrow(div_intersect_qtl) = number of 1kb bins with signidicant qtl
# 
# # over representation:
# phyper(nrow(div_intersect_qtl), 
#        nrow(d_regions), 
#        nrow(n_1kb_bins)-nrow(d_regions), 
#        nrow(div_not_qtl)+nrow(div_intersect_qtl), 
#        lower.tail= FALSE)
# 
# # [1] 1
# 
# # depletion:
# phyper(nrow(div_intersect_qtl), 
#        nrow(d_regions), 
#        nrow(n_1kb_bins)-nrow(d_regions), 
#        nrow(div_not_qtl)+nrow(div_intersect_qtl), 
#        lower.tail= T)
# 
# # [1] 2.285546e-16
# 
# # over representation:
# phyper(nrow(div_intersect_qtl), 
#        nrow(d_regions), 
#        nrow(n_1kb_bins)-nrow(d_regions), 
#        nrow(div_not_qtl)+nrow(div_intersect_qtl), 
#        lower.tail= FALSE)
# 
# # [1] 1.135497e-26
# 
# # depletion:
# phyper(nrow(div_intersect_qtl), 
#        nrow(d_regions), 
#        nrow(n_1kb_bins)-nrow(d_regions), 
#        nrow(div_not_qtl)+nrow(div_intersect_qtl), 
#        lower.tail= T)
# 
# # [1] 1