library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/haplotype_plot_df.Rda")
load("Processed_Data/indep_strain_info.RData")
load("Processed_Data/PCA.RData")
load("Processed_Data/umap.RData")
load("Processed_Data/strain_geo.RData")

df_chr_length <- data.table::fread("Processed_Data/chr_info.tsv") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::rename(CHROM = V1, stop= V2) %>%
  dplyr::mutate(start=0)

### Defining swept chromosomes and swept haplotype

plot_df_sweep_summary <- plot_df %>%
  dplyr::distinct(isotype, chromosome, filtered_sweep_ratio) %>%
  dplyr::group_by(chromosome) %>%
  dplyr::summarise(mean_sweep_ratio=mean(filtered_sweep_ratio)) %>%
  dplyr::ungroup()

### Defining swept strains

df_swept_isotype_gw <- plot_df %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::distinct(isotype, chromosome, filtered_sweep_ratio) %>%
  dplyr::group_by(isotype, chromosome) %>%
  dplyr::summarise(max_swept_ratio = max(filtered_sweep_ratio), 
                   mean_swept_ratio = mean(filtered_sweep_ratio),
                   min_swept_ratio = min(filtered_sweep_ratio)) %>%
  dplyr::arrange(max_swept_ratio) %>%
  dplyr::mutate(is_swept = ifelse(max_swept_ratio > 0.3, "swept", "unswept"))


chrom_length <- sum(dplyr::filter(df_chr_length, !CHROM %in% c("II","III"))$stop)

df_swept_isotype <- plot_df %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::distinct(isotype, chromosome, filtered_sweep_len) %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(filtered_swept_fraction = sum(filtered_sweep_len)/chrom_length) %>%
  dplyr::ungroup() %>%
  dplyr::rename(strain=isotype)

save(df_swept_isotype, file = "Processed_Data/df_swept_isotype.RData")


### overlay with pop structure

### PCA with complete site vcf

plot_PC12_all_sweep <- df_pcs_noremoval_complete %>%
  dplyr::filter(geo != "Unknown") %>%
  dplyr::left_join(., df_swept_isotype, by='strain') %>%
  ggplot(.) +
  #geom_text_repel(data=dplyr::filter(df_pcs_noremoval_clean, PC1 < -0.2 | PC2 >0.07),  
  #                label=dplyr::filter(df_pcs_noremoval_clean, PC1 < -0.2 | PC2 >0.07)$isotype, 
  #                aes(x=PC1,y=PC2)) +  
  geom_point(shape=21, aes(x=PC1, y=PC2, fill = geo, size=filtered_swept_fraction, alpha=-filtered_swept_fraction))+
  scale_fill_manual(values=geo.colours) +
  scale_size_continuous(range = c(1,3), trans = 'exp', guide=FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position='none',
        panel.grid = element_blank())+
  labs(x="PC1", y="PC2", fill ="")  +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  guides(fill= guide_legend(nrow=2))

plot_PC12_all_sweep

##### outlier removed

plot_PC12_pruned_sweep <- df_pcs_out_removal_complete %>%
  dplyr::filter(geo != "Unknown") %>%
  dplyr::left_join(., df_swept_isotype, by='strain') %>%
  ggplot(.) +
  geom_point(shape=21, aes(x=PC1, y=PC2, fill = geo, size=filtered_swept_fraction, alpha =-log(filtered_swept_fraction)))+
  scale_fill_manual(values=geo.colours) +
  scale_size_continuous(range = c(0.5,3), trans = 'exp', guide=FALSE) +  
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.title = element_blank(),
        legend.position ='none',
        panel.grid = element_blank())+
  labs(x="PC1", y="PC2", fill ="")

plot_PC12_pruned_sweep


## UMAP

plot_umap_sweep <- temp_umap %>%
  dplyr::filter(geo != "Unknown") %>%
  dplyr::left_join(., df_swept_isotype, by='strain') %>%
  ggplot(.)+
  aes(x=X1,y=X2)+
  geom_point(shape=21, aes(fill = geo, size=filtered_swept_fraction, alpha=-filtered_swept_fraction))+
  scale_fill_manual(values=geo.colours) +
  scale_size_continuous(range = c(0.5,3), trans = 'exp', guide=FALSE) +  
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title=element_blank(),
        axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position ='none') +
  labs(x = "UMAP1", y = "UMAP2")

plot_umap_sweep

plot_PC1_all_sweep <- df_pcs_noremoval_complete %>%
  dplyr::left_join(., df_swept_isotype, by='strain') %>%
  ggplot(.) +
  geom_point(shape=21, size = 1.5, alpha = 0.8)+
  aes(x=PC1, y=filtered_swept_fraction, fill=geo) +
  scale_fill_manual(values=geo.colours) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title=element_blank(),
        axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position ='none') +
  labs(y="Fraction swept")+
  scale_x_continuous(breaks = c(0, 0.1, 0.2))+
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))

plot_PC1_pruned_sweep <- df_pcs_out_removal_complete %>%
  dplyr::left_join(., df_swept_isotype, by='strain') %>%
  ggplot(.) +
  geom_point(shape=21, size = 1.5, alpha = 0.8)+
  aes(x=PC1, y=filtered_swept_fraction, fill=geo) +
  scale_fill_manual(values=geo.colours) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title=element_blank(),
        axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position ='none') +
  labs(y="Fraction swept")+
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))

plot_umap1_sweep <- temp_umap %>%
  dplyr::left_join(., df_swept_isotype, by='strain') %>%
  ggplot(.) +
  geom_point(shape=21, aes(fill = geo), size = 1.5, alpha = 0.8)+
  aes(x=X1, y=filtered_swept_fraction) +
  scale_fill_manual(values=geo.colours) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title=element_blank(),
        axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position ='none') +
  labs(y="Fraction swept", x="UMAP1")+
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))


PCA_sweep_merged <- cowplot::plot_grid(plot_PC12_all_sweep, plot_PC12_pruned_sweep, plot_umap_sweep,
                                       plot_PC1_all_sweep, plot_PC1_pruned_sweep, plot_umap1_sweep,
                                       labels = c("a","b","c","d","e","f"), 
                                       label_size = 12, nrow=2,
                                       align = 'h', axis='tblr')  

PCA_sweep_merged

ggsave(PCA_sweep_merged, file="PCA_sweep_merged.pdf", width = 7.5, height=2.5, unit='in')

############# Supple Fig. xx Selective sweeps #####################################

#==============================#
# Plot ~ Swept haplotypes only #
#==============================#

strain_labels <- plot_df %>%
  dplyr::ungroup() %>%
  dplyr::select(plotpoint, isotype) %>%
  dplyr::distinct() %>%
  dplyr::arrange(plotpoint)

df_sweep_all <- indep_strain_info_geo %>%
  dplyr::filter(reference_strain == 1) %>%
  dplyr::distinct(isotype, geo) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::distinct(isotype, geo, chromosome, filtered_sweep_ratio)

df_sweep_all_summary <- indep_strain_info_geo %>%
  dplyr::filter(reference_strain == 1) %>%
  dplyr::distinct(isotype, geo) %>%
  dplyr::arrange(as.character(geo)) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::distinct(isotype, geo, chromosome, isotype_swept_haplotype_length) %>%
  dplyr::group_by(isotype, geo) %>%
  dplyr::summarise(sum_sweep = sum(isotype_swept_haplotype_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(isotype, geo, sum_sweep) %>%
  dplyr::group_by(geo) %>%
  dplyr::mutate(mean_sweep = mean(sum_sweep)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(mean_sweep, sum_sweep)) %>%
  dplyr::mutate(plotpoint=row_number())

plot_df_sweep_all <- plot_df %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::left_join(., df_sweep_all_summary, by = 'isotype')

## sweep in all strains

plot_df_sweep_all$geo <- factor(plot_df_sweep_all$geo, levels = c("Asia", "S. America", "Austrailia", "Africa","Europe",
                                                                  "N. America","New Zealand","Atlantic","Hawaii","Unknown"))

plot_sweep_all_hap_hor <- plot_df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III"), geo != "Unknown") %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$geo,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=10, color = "black"), 
        axis.text.y  = element_blank(), 
        strip.text.x = element_text(size=9, color = "black"), 
        legend.position = 'none', axis.title.x = element_text(size=11, color = "black"), 
        axis.ticks.y = element_blank(), strip.text.y=element_text(size=8, color = "black", angle = 90),  
        strip.text.y.left = element_text(angle=0),
        panel.spacing = unit(0.05, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in")) +
  facet_grid(geo~chromosome, scales="free", space="free_y", switch = 'y')

plot_sweep_all_hap_hor

plot_sweep_all_hor <- df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III"), geo != "Unknown") %>%
  ggplot(.) +
  #geom_jitter(alpha=0.8, size = 0.6) +
  #geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  geom_quasirandom(shape=21, alpha = 0.7) +
  aes(x=reorder(geo, filtered_sweep_ratio), y= filtered_sweep_ratio, fill=geo) +
  scale_fill_manual(values = geo.colours) +
  labs(x="Ancestry", y="Fraction swept")+
  scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1))+
  #scale_x_discrete(limits = rev(levels(df_unadmix_sweep$geo)))+
  coord_flip()+
  theme_bw() +
  theme(axis.title = element_text(size=11, color='black'), legend.position = 'none', 
        axis.text = element_text(size=10, color='black'),
        axis.title.y = element_blank(), strip.background = element_blank(), 
        strip.text = element_blank(),  panel.spacing = unit(0.1, "lines"), 
        plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in"),
        panel.grid = element_blank()) +
  facet_wrap(~chromosome, nrow=1)

plot_sweep_all_hor

plot_sweep_all <- cowplot::plot_grid(plot_sweep_all_hap_hor, plot_sweep_all_hor, 
                                     ncol=1, rel_heights = c(3,1), align = "hv", axis = "lr",
                                     labels = c("a", "b"), label_size = 12)

plot_sweep_all

ggsave(plot_sweep_all, file = "Plots/Supplementary/sweep_geo.pdf", 
       width = 7.5,
       height = 11, unit = "in")

e## sweep in unswept strains for various threshold (mean swept ratio = 5, 10, 20)

plot_df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III"), isotype %in% sweep10_isotype) %>%
  dplyr::left_join(., plotpoint10, by='isotype') %>%
  #dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=11, color = "black"), axis.text.y  = element_blank(), strip.text = element_text(size=11, color = "black"), legend.position = 'none', axis.title.x = element_text(size=13, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in")) +
  facet_grid(~chromosome, scales="free", space="free_y", switch = "y")

ggsave("Plots/sweep10_isotypes.png", width = 7.5, height = 5, units = "in")

plot_df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III"), isotype %in% sweep20_isotype) %>%
  dplyr::left_join(., plotpoint20, by='isotype') %>%
  #dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=11, color = "black"), axis.text.y  = element_blank(), strip.text = element_text(size=11, color = "black"), legend.position = 'none', axis.title.x = element_text(size=13, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in")) +
  facet_grid(~chromosome, scales="free", space="free_y", switch = "y")

ggsave("Plots/sweep20_isotypes.png", width = 7.5, height = 5, units = "in")

plot_df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III"), isotype %in% sweep30_isotype) %>%
  dplyr::left_join(., plotpoint30, by='isotype') %>%
  #dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=11, color = "black"), axis.text.y  = element_blank(), strip.text = element_text(size=11, color = "black"), legend.position = 'none', axis.title.x = element_text(size=13, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in")) +
  facet_grid(~chromosome, scales="free", space="free_y", switch = "y")

ggsave("Plots/sweep30_isotypes.png", width = 7.5, height = 5, units = "in")




## sweep in N2

plot_df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::filter(isotype == "N2") %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=11, color = "black"), axis.text.y  = element_blank(), strip.text = element_text(size=11, color = "black"), legend.position = 'none', axis.title.x = element_text(size=13, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in")) +
  facet_grid(cluster~chromosome, scales="free", space="free_y", switch = "y")


### Fig. S7 Pacific Rim strains, sweep ######

world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica

########### Plotting for the identification of geographic structure

plot_east_pacificrim <- ggplot()+ geom_map(data=world, map=world,
                                           aes(x=long, y=lat, map_id=region),
                                           color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat, fill = 'red'), shape =21, size =3, alpha = 0.9) +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-140, -70)) 

plot_east_pacificrim ## take QG2075, EG4349, EG4946 out

plot_west_pacificrim <- ggplot()+ geom_map(data=world, map=world,
                                           aes(x=long, y=lat, map_id=region),
                                           color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat, fill = 'red'), shape =21, size =3, alpha = 0.9) +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(120, 200)) 

plot_west_pacificrim

## Austrailian strains

plot_AU <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat, fill = 'red'), shape =21, size =3, alpha = 0.9) +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(110, 160), y=c(-50,-10)) 

plot_AU

df_au_strain <- indep_strain_info %>%
  dplyr::filter(long > 110  & long < 160) %>%
  dplyr::filter(lat > -50 & lat < -10) 

au_strain <- as.character(df_au_strain$strain)
au_isotype <- unique(as.character(df_au_strain$isotype))

plot_au_confirm <- ggplot()+ geom_map(data=world, map=world,
                                      aes(x=long, y=lat, map_id=region),
                                      color="black", fill="white", size=0.2)+
  geom_point(data = dplyr::filter(indep_strain_info, strain %in% au_strain), aes(x=long, y=lat, fill = 'red'), shape =21, size =3, alpha = 0.9) +
  #geom_label_repel(aes(long, lat, label = strain), data = dplyr::filter(indep_strain_info, strain %in% au_strain), fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12))

plot_au_confirm

## Hawaiian strains

plot_hw <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat, fill = 'red'), shape =21, size =3, alpha = 0.9) +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-170, -150)) 

plot_hw

df_hw_strain <- indep_strain_info %>%
  dplyr::filter(long > -170  & long < -150)

hw_strain <- as.character(df_hw_strain$strain)
hw_isotype <- unique(as.character(df_hw_strain$isotype))

plot_hw_confirm <- ggplot()+ geom_map(data=world, map=world,
                                      aes(x=long, y=lat, map_id=region),
                                      color="black", fill="white", size=0.2)+
  geom_point(data = dplyr::filter(indep_strain_info, strain %in% hw_strain), aes(x=long, y=lat, fill = 'red'), shape =21, size =3, alpha = 0.9) +
  #geom_label_repel(aes(long, lat, label = strain), data = dplyr::filter(indep_strain_info, strain %in% hw_strain), fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12))

plot_hw_confirm

## pacific rim strain list

df_pr_strain <- indep_strain_info %>%
  dplyr::filter(long > 120 | long < -70 & long > -140, !strain %in% c("QG2075", "EG4349", "EG4946", au_strain)) 

pr_strain <- as.character(df_pr_strain$strain)
pr_isotype <- unique(as.character(df_pr_strain$isotype))

plot_pacificrim <- ggplot()+ geom_map(data=world, map=world,
                                      aes(x=long, y=lat, map_id=region),
                                      color="gray42", fill="white", size=0.2)+
  geom_point(data = dplyr::filter(indep_strain_info, strain %in% pr_strain), aes(x=long, y=lat, fill = 'red'), shape =21, size =1.5, alpha = 0.9) +
  #geom_label_repel(aes(long, lat, label = strain), data = dplyr::filter(indep_strain_info, strain %in% pr_strain), fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  geom_point(data = dplyr::filter(df_subpop, strain %in% pr_strain), aes(x=long, y=lat, fill = cluster), shape =21, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = ancestry.colours) +
  scale_color_manual(values = ancestry.colours) +
  theme_map()+
  theme(legend.position = 'none')

plot_pacificrim

#ggsave("Plots/plot_pacificrim.pdf", plot = plot_pacificrim, width = 7.5, height = 3.6, units = "in")

### pr strain info

df_pr_info <- df_subpop %>%
  dplyr::filter(isotype %in% pr_isotype)

plot_pr_admix <- df_pr_info %>%
  ggplot(.) +
  geom_bar() +
  aes(x=cluster) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=11, color='black'), axis.text = element_text(size=10, color='black')) +
  labs(y="Strain count")

## pr strains haplotype

df_sweep_pr <- df_subpop %>%
  dplyr::filter(isotype %in% pr_isotype) %>%
  dplyr::distinct(isotype, cluster) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::distinct(isotype, cluster, chromosome, filtered_sweep_ratio)

df_sweep_pr_summary <- df_subpop %>%
  dplyr::filter(isotype %in% pr_isotype) %>%
  dplyr::distinct(isotype, cluster) %>%
  dplyr::arrange(as.character(cluster)) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::distinct(isotype, cluster, chromosome, isotype_swept_haplotype_length) %>%
  dplyr::group_by(isotype, cluster) %>%
  dplyr::summarise(sum_sweep = sum(isotype_swept_haplotype_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(isotype, cluster, sum_sweep) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(mean_sweep = mean(sum_sweep)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(mean_sweep, sum_sweep)) %>%
  dplyr::mutate(plotpoint=row_number())

plot_df_sweep_pr <- plot_df %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::filter(isotype %in% pr_isotype) %>%
  dplyr::left_join(., df_sweep_pr_summary, by = 'isotype')

plot_sweep_pr_hap_hor <- plot_df_sweep_pr %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_pr$isotype)),
                     labels = df_sweep_pr_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=11, color = "black"), axis.text.y  = element_blank(), strip.text = element_text(size=11, color = "black"), legend.position = 'none', axis.title.x = element_text(size=13, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in"), strip.background.y = element_rect(colour = "transparent", fill = "transparent")) +
  facet_grid(cluster~chromosome, scales="free", space="free_y", switch = "y")

plot_sweep_pr_hap_hor

save(plot_df_sweep_pr, file = "Processed_Data/Pacific_rim_sweep.RData")

FigS7ab <- cowplot::plot_grid(plot_pacificrim, plot_pr_admix, ncol=1, rel_heights = c(4,3), labels = c("a", "b"))
FigS7c <- cowplot::plot_grid(plot_sweep_pr_hap_hor, labels = c("c"))

FigS7 <- cowplot::plot_grid(FigS7ab, FigS7c, nrow=1, rel_widths = c(1,1))

FigS7

ggsave("Plots/FigS7_PacificRim.pdf", plot = FigS7, width = 10, height = 5, units = "in")

### pacific vs non-pacific isotypes
## pacfic & non-pacific isotypes - ECA251, MY23

pr_isotype
hw_isotype
pr_hw_isotype <- c(pr_isotype,hw_isotype)


df_isotype_admix_Pr <- df_isotype_admix %>%
  dplyr::filter(isotype %in% pr_isotype)

unique(df_isotype_admix_Pr$cluster)

nrow(dplyr::filter(df_isotype_admix_Pr, cluster %in% c("WD","PA","HW1","HW2","HW3")))/nrow(df_isotype_admix_Pr)
nrow(dplyr::filter(df_isotype_admix_Pr, cluster %in% c("WS1","WS2","CA","PT1","AT", "AU1", "AU2")))/nrow(df_isotype_admix_Pr)


df_isotype_admix_PA <- df_isotype_admix %>%
  dplyr::filter(isotype %in% pr_hw_isotype)

nrow(dplyr::filter(df_isotype_admix_PA, cluster %in% c("WD","PA","HW1","HW2","HW3")))/nrow(df_isotype_admix_PA)

df_isotype_admix_nonPA <- df_isotype_admix %>%
  dplyr::filter(!isotype %in% c(pr_hw_isotype, "CB4852", "ECA259", "PB303", "XZ2018") | isotype %in% c("MY23", "ECA251"))

nrow(dplyr::filter(df_isotype_admix_nonPA, cluster %in% c("WS1","WS2","EU1","EU2","FR","AU1","AU2","PT1","AT")))/nrow(df_isotype_admix_nonPA)
