


load("Processed_Data/admixture_full.RData")
load("Processed_Data/admix_K16.RData")
load("Processed_Data/haplotype_plot_df.Rda")

load("Processed_Data/Ce_Genome-wide_Neutrality_stats_K16.Rda")


load("Processed_Data/divergent_joined.RData")
load("Processed_Data/divergent_count_summary.RData")
load("Processed_Data/ILS_trees.RData")

load("Processed_Data/pals_div.RData")
load("Processed_Data/LD_decay.RData")
load("Processed_Data/divergent_hap_zoom.RData")
#load("Processed_Data/ILS_trees.RData")
#load("Processed_Data/divergent_freq_popgenome_example.Rdata")
load("Processed_Data/divergent_gene_freq_popgenome_example.Rdata")
load("Processed_Data/ortho_branch_length.RData")
load("Processed_Data/ortho_id.RData")
load("Processed_Data/Cbr_divergent.RData")

load("Processed_Data/admixture_full.RData")
load("Processed_Data/treemix_ve.RData")
load("Processed_Data/div_optimization.RData")
load("Processed_Data/Divergent_chromosome_stats.RData")
load("Processed_Data/Pop_div_sweep.RData")
load("Processed_Data/Divergent_haplotype_examples.Rdata")
load("Processed_Data/swept_isotype_admixutre.RData")
load("Processed_Data/div_haplo_counts.RData")


ct=11
cov=0.15
cluster_threshold = 15000


########## Main figures ############

##### Fig. 1 The pacific origin of C. elegans ######

### Fig. 1b 343 independent isolation (with labels for non-admixted strains) on the world map ###

admix_strain_info <- indep_strain_info %>%
  dplyr::filter(isotype %in% unadmix_isotype) %>%
  dplyr::left_join(., df_unadmix, by = 'isotype')

plot_map_world <- 
  ggplot()+ geom_map(data=world, map=world,
                     aes(x=long, y=lat, map_id=region),
                     color="gray51", fill="white", size=0.2, alpha = 0.8)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size = 1, alpha = 0.7, fill = "red") +
  #ggsn::scalebar(data = admix_strain_info, dist = 2000, dist_unit = "km", location = "bottomleft",
  #           transform = TRUE, model = "WGS84", height = 0.02, st.size = 2.3, st.dist	= 0.03, border.size = 0.3) +
  #geom_label_repel(aes(long, lat, label = strain, fill = cluster), data = admix_strain_info, fontface = 'bold', color = 'white', size = 0, box.padding = 0.2, point.padding = 0.3, segment.color = '#9acd32', segment.size = 0.3, segment.alpha = 0.6, nudge_y = 0.05) +
  theme_map()+
  theme(text = element_text(size=12), legend.position = 'bottom', legend.direction = "horizontal", legend.title = element_blank(), legend.key.size = unit(0.05, "in"),  legend.key.width = unit(0.05,"in")) +
  guides(fill= guide_legend(nrow=2))

plot_map_world

Fig1a <- cowplot::plot_grid(plot_map_world, labels = c("a"), label_size = 12)  
ggsave("Plots/main_figures/Fig1a.pdf", Fig1a, width = 7.5, height = 3.5)

####### Fig. 1b - Genome-wide tree of unadmixed strains (grouped by admixture ancestry)

plot_tree_NJ_admix <- ggtree(tree_ph_admix, layout="fan", open.angle=60, branch.length="rate") + 
  aes(size = label) + 
  theme(legend.position = c(0.68,0.4), legend.text = element_text(size = 13), legend.title = element_blank(), plot.margin = margin(b=-7, l=-8, r=-1.6, t=-8, unit = "in")) +
  guides(col= guide_legend(ncol=9))  +
  scale_size_manual(values=c(.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.2))

plot_tree_NJ_admix_rotate <- rotate_tree(plot_tree_NJ_admix,90) %>% flip(639, 642) %>% flip(579,585)

plot_tree_NJ_admix_rotate

Fig1a <- cowplot::plot_grid(plot_tree_NJ_admix_rotate, labels = c("b"), label_size = 12)  
ggsave("Plots/Fig1b_NJ.pdf", Fig1b, width = 7.5, height = 5)



####### Fig. 1a - Genome-wide tree of unadmixed strains (grouped by admixture ancestry)

plot_tree_NJ_admix <- ggtree(tree_ph_admix, layout="fan", open.angle=60, branch.length="rate") + 
  aes(color = label, size = label) + 
  theme(legend.position = c(0.68,0.4), legend.text = element_text(size = 13), legend.title = element_blank(), plot.margin = margin(b=-7, l=-8, r=-1.6, t=-8, unit = "in")) +
  #theme(legend.position = 'none', legend.text = element_text(size = 10), legend.title = element_blank(), plot.margin = margin(b=-8, l=-10, r=-6, t=-10, unit = "in")) +
  scale_color_manual(values = c("WD"="#ff38d1", "WS1"="#FF0009","HW1"= "darkorange1", 
                                "EU1"="lightskyblue2", "AU1"="#57412F","CA"= "burlywood3", "PA"="turquoise", 
                                "PT2"="springgreen4", "AU2"="lightpink2", "HW3"="#C24CF6", "EU2"="black", 
                                "FR"="#0047ab","AT"= "#9b351b","HW2"= "maroon","PT1"= "yellow3","WS2"= "#39FF14", "Z_Admixed" = "gray")) + 
  guides(col= guide_legend(ncol=9))  +
  scale_size_manual(values=c(.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.2))

plot_tree_NJ_admix_rotate <- rotate_tree(plot_tree_NJ_admix,90) %>% flip(639, 642) %>% flip(579,585)

plot_tree_NJ_admix_rotate

Fig1a <- cowplot::plot_grid(plot_tree_NJ_admix_rotate, labels = c("a"), label_size = 12)  
ggsave("Plots/Fig1a_NJ.pdf", Fig1a, width = 7.5, height = 5)


### Fig. 1b 343 independent isolation (with labels for non-admixted strains) on the world map ###

admix_strain_info <- indep_strain_info %>%
  dplyr::filter(isotype %in% unadmix_isotype) %>%
  dplyr::left_join(., df_unadmix, by = 'isotype')

plot_map_world_unadmix <- 
  ggplot()+ geom_map(data=world, map=world,
                     aes(x=long, y=lat, map_id=region),
                     color="gray51", fill="white", size=0.2, alpha = 0.8)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size = 1, alpha = 0.7, fill = "red") +
  ggsn::scalebar(data = admix_strain_info, dist = 2000, dist_unit = "km", location = "bottomleft",
                 transform = TRUE, model = "WGS84", height = 0.02, st.size = 2.3, st.dist	= 0.03, border.size = 0.3) +
  geom_label_repel(aes(long, lat, label = strain, fill = cluster), data = admix_strain_info, fontface = 'bold', color = 'white', size = 0, box.padding = 0.2, point.padding = 0.3, segment.color = '#9acd32', segment.size = 0.3, segment.alpha = 0.6, nudge_y = 0.05) +
  scale_fill_manual(values = ancestry.colours) +
  theme_map()+
  theme(text = element_text(size=12), legend.position = 'bottom', legend.direction = "horizontal", legend.title = element_blank(), legend.key.size = unit(0.05, "in"),  legend.key.width = unit(0.05,"in")) +
  guides(fill= guide_legend(nrow=2))

plot_map_world_unadmix 

Fig1b <- cowplot::plot_grid(plot_map_world_unadmix, labels = c("b"), label_size = 12)  
ggsave("Plots/Fig1b_scalebar.pdf", Fig1b, width = 7.5, height = 3.5)

##### Fig. 2 Discovery of hyper-divergent regions #####

### Fig. 2a Divergent region landscape ### pals, peel-1/zeel-1, sup-35/pha-1 labeled

plot_div_pop <- df_join_div %>%
  na.omit() %>%
  ggplot(.) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=1.01), color='transparent', fill='transparent', size =0.1) +
  #geom_segment(data=dplyr::distinct(df_pals_div,div_start, .keep_all=T), aes(x=(start+stop)/2e6, xend=(start+stop)/2e6, y=0.4, yend=327.6), color='red', alpha=0.75, size=0.1) +
  #geom_segment(data=data.frame(CHROM="I", start=2342215, stop=2356238), aes(x=start/1e6, xend=stop/1e6, y=0.4, yend=327.6), color='blue', alpha=0.75, size=0.1) + ## peel-1, zeel-1
  #geom_segment(data=data.frame(CHROM="III", start=11116993, stop=11136380), aes(x=start/1e6, xend=stop/1e6, y=0.4, yend=327.6), color='purple', alpha=0.75, size=0.1) + ## sup-35, pha-1
  geom_rect(aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = plotpoint-0.5 , ymax = plotpoint+0.5), fill = 'black',color='black', size = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y  = element_blank(), strip.text = element_text(size=10, color = "black"), legend.position = 'none',axis.title.y = element_text(size=11, color = "black"), axis.title.x=element_blank(), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0.00, 0.00), limits=c(0.4,327.6)) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  facet_grid(~CHROM, scales="free",space = 'free') +
  labs(x="Genomic position (Mb)",y="327 wild isotypes")

plot_div_pop

### Fig. 2b Divergent region landscape (merged and joined) + pals overlap ### 

df_all_cluster_size_freq <- df_divergent_final %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, freq) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  join_masks(.) %>%
  na.omit() %>%
  dplyr::group_by(CHROM, cluster_start, cluster_end, cluster_size) %>%
  dplyr::summarise(Frequency = mean(freq)) %>%
  dplyr::ungroup()

mean(df_all_cluster_size_freq$Frequency)
sd(df_all_cluster_size_freq$Frequency)
mean(df_all_cluster_size_freq$cluster_size)
sd(df_all_cluster_size_freq$cluster_size)

plot_div_pop_freq <- df_all_cluster_size_freq %>%
  na.omit() %>%
  ggplot(.) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 0.2, ymax=0.3), color='transparent', fill='transparent', size =0.1) +
  #geom_segment(data=dplyr::distinct(df_pals_div,div_start, .keep_all=T), aes(x=(start+stop)/2e6, xend=(start+stop)/2e6, y=-0.01, yend=0.43), color='red', alpha=0.75, size=0.1) +
  #geom_segment(data=data.frame(CHROM="I", start=2342215, stop=2356238), aes(x=start/1e6, xend=stop/1e6, y=-0.01, yend=0.43), color='blue', alpha=0.75, size=0.1) + ## peel-1, zeel-1
  #geom_segment(data=data.frame(CHROM="III", start=11116993, stop=11136380), aes(x=start/1e6, xend=stop/1e6, y=-0.01, yend=0.43), color='purple', alpha=0.75, size=0.1) + ## sup-35, pha-1
  geom_rect(aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = Frequency-0.01, ymax = Frequency+0.01), color='black', fill='black', size = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=10, color = "black"), legend.position = 'none', axis.title.x = element_text(size=11, color = "black"), axis.title.y = element_text(size=11, color = "black"), panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in"), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_blank()) +
  labs(x="Genomic position (Mb)", y="Frequency") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.02,0.57)) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  facet_grid(~CHROM, scales="free", space = 'free')

plot_div_pop_freq

### Fig. 2c Divergent fraction x variant count fraction ###

plot_divergent_strain_size_varfrac <- df_divergent_count_summary %>%
  ggplot(.) +
  geom_point(aes(x=div_fraction*100, y=variant_fraction*100), size = 0.6, alpha =0.8) +
  geom_text_repel(data=dplyr::filter(df_divergent_count_summary, variant_fraction>0.35), aes(x=div_fraction*100, y=variant_fraction*100, label = STRAIN), size = 2.7, box.padding = 0.2, point.padding = 0.1, segment.size = 0.02, nudge_x = -0.1, nudge_y= 0, segment.alpha = 0) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color='black'), axis.text = element_text(size=10, color='black'), panel.grid=element_blank()) +
  labs(x="Hyper-divergent regions (%)", y="Variants (%)",plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in")) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(expand=c(0.07,0.07))

plot_divergent_strain_size_varfrac

df_divergent_count_summary2 <- df_divergent_count_summary %>%
  dplyr::mutate(div_var_fraction = variant_fraction/div_fraction) %>%
  dplyr::arrange(-div_var_fraction)

Fig2 <- cowplot::plot_grid(plot_div_pop, plot_div_pop_freq, plot_divergent_strain_size_varfrac, labels = c("a","b", "c"), rel_heights = c(2,1,1.3), ncol=1, label_size = 12,
                           align="v", axis='lr')

ggsave("Plots/main_figures/Fig2.png", Fig2, width = 7.5, height = 8)

##### Fig. 3 Common divergent regions are shared among divergent strains #####

### Fig. 3a LDdecay ###

df_LDdecay$class <- factor(df_LDdecay$class, leves = c("All", "Divergent", "Non-divergent"), labels = c("All", "Hyper-divergent", "Non-divergent"))

plot_LDdecay <- df_LDdecay %>%
  dplyr::filter(class != "All", dist >= 10) %>%
  ggplot(.) +
  aes(x=dist/1e3, y=mean_r2) +
  geom_point(size = 0.2, alpha =0.6, aes(color = class)) +
  geom_smooth(se = F, size = 0.7, n = 5000, aes(group=class, linetype = class), color='black') + ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
  scale_color_brewer(palette = 'Set1') +
  scale_linetype(guide=FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'), 
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.y= element_text(size=10, color='black'), 
        axis.text.x = element_text(size=10, color='black'), 
        legend.position = c(0.75,0.85), panel.grid = element_blank(),
        legend.text =element_text(size=10, color='black'),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.08,"in"),
        plot.margin = unit(c(t=0.1,r=0.1,b=0.1,l=0.2), "in")) +
  labs(x = "Distance (kb)", y = "Linkage\ndisequilibrium (r2)", color = "") +
  guides(colour = guide_legend(override.aes = list(size=1.5, alpha = 1)))

plot_LDdecay

### Fig. 3b 

### Fig. 3b Divergent region haplotype zoom in example (V:391000-739000) ###

plot_div_hap_zoom <- divergent_hap_overlay_unadmix(df2=plot_df_unadmix, chr = "V", left=391000, right=739000, threshold = cluster_threshold)

plot_div_hap_zoom

### Fig. 3c Tajima's D ###

df_divergent_freq_simple <- df_divergent_final %>%
  dplyr::mutate(startWindow=START_BIN+1, endWindow=END_BIN) %>%
  dplyr::distinct(CHROM, startWindow, endWindow,freq_bin)

neutrality_df_div <- neutrality_df %>%
  dplyr::left_join(., df_divergent_freq_simple, by=c('CHROM','startWindow','endWindow')) %>%
  dplyr::mutate(freq_bin=ifelse(is.na(freq_bin), "Non-divergent", as.character(freq_bin)))

neutrality_df_div$freq_bin <- factor(neutrality_df_div$freq_bin, levels = c("Non-divergent","Rare","Intermediate", "Common"), labels=  c("Not\ndivergent","< 1%","1-5%", "≥ 5%"))
neutrality_df_div$statistic <- factor(neutrality_df_div$statistic, levels = c("theta_Watterson","Pi","Tajima.D"), labels = c("Watterson's estimator (θ)","Nucleotide diversity (π)","Tajima's D"))

tajimad_div <- neutrality_df_div %>%
  dplyr::filter(statistic=="Tajima's D") %>%
  na.omit()

plot_div_TajD <- tajimad_div %>%
  na.omit() %>%
  ggplot(.) +
  #geom_violin() +
  geom_boxplot(outlier.alpha = 0) +
  #geom_jitter(size=0.1, alpha = 0.2) +
  aes(x=freq_bin, y=value) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=11, color='black'),
        axis.text.x = element_text(size=10, color='black'), 
        axis.text.y = element_text(size=10, color='black'),
        panel.grid = element_blank()) +
  labs(y="Tajima's D") +
  coord_cartesian(ylim=c(-3,3.5))

plot_div_TajD



### Fig. 3e Wormcat - layer 2 / all vs control ###

plot_wormcat_2 <- df_cat2 %>%
  dplyr::left_join(., df_cat2_plotpoint, by=c("Category", "layer")) %>%
  ggplot(.) +
  geom_point(alpha=0.9) +
  geom_vline(xintercept = -log10(0.05), color='blue', size=0.4) +
  aes(x=-log10(Bonferroni), y=plotpoint, shape=group, size=EnrichRatio, fill=RGS) +
  theme_bw() +
  scale_y_continuous(breaks = 1:length(unique(df_cat2_plotpoint$Category)), labels = unique(df_cat2_plotpoint$Category), expand=c(0.07,0.07)) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22), guide=FALSE) +
  scale_size_continuous(range = c(1,4), breaks = c(2,3,4)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title.x = element_text(size=11, color='black'), axis.title.y = element_blank(), legend.title = element_text(size=10, color='black'), legend.text = element_text(size=9.5, color='black'), legend.margin = margin(0,0,0,0), legend.direction = "horizontal", legend.box = "vertical", legend.position = c(0.6,0.65), legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.2, "in"), panel.grid.minor.y = element_blank(), legend.spacing.y = unit(0.05, 'in'), panel.grid = element_blank()) +
  guides(shape= guide_legend(nrow=2, order = 3, override.aes = list(size=3), title = NULL, title.position = "top"), size= guide_legend(nrow=1, order = 1, title.position = "top"), fill =  guide_colourbar(nrow=1, order =2, title.position = "top")) +
  labs(x="-log10(corrected p-value)", shape = "Frequency", size = "Fold enrichment", fill = "Gene counts") 

plot_wormcat_2


Fig3ab <- cowplot::plot_grid(plot_div_common_fraction_admix, plot_LDdecay, labels = c("a","b"),ncol=1, align = 'v', label_size = 12)
Fig3c <- cowplot::plot_grid(plot_div_hap_zoom, labels = c("c"), label_size = 12)
Fig3abc <- cowplot::plot_grid(Fig3ab, Fig3c, rel_widths = c(1,1), nrow=1, label_size = 12, align = 'hv', axis='lr') 
Fig3de <- cowplot::plot_grid(plot_div_TajD, plot_wormcat_2, labels = c("d","e"), rel_widths = c(1,3), label_size = 12)
Fig3 <- cowplot::plot_grid(Fig3abc, Fig3de, rel_heights = c(2,1.2), ncol=1, label_size = 12) 

ggsave("Plots/main_figures/Fig3.png", Fig3, width = 7.5, height = 8)




## old version cowplot

Fig3ab <- cowplot::plot_grid(plot_div_common_fraction_admix, plot_div_TajD, labels = c("a","b"),nrow=1, rel_widths = c(2,1), align = 'h', label_size = 12)
Fig3c <- cowplot::plot_grid(plot_div_hap_zoom, labels = c("c"), label_size = 12)
Fig3d <- cowplot::plot_grid(plot_LDdecay, labels = c("d"), label_size = 12)
Fig3abd <- cowplot::plot_grid(Fig3ab, Fig3d, rel_heights = c(1,1), ncol=1, label_size = 12, align = 'hv', axis='lr') 
Fig3e <- cowplot::plot_grid(plot_wormcat_2, labels = c("e"), label_size = 12)
Fig3abde <- cowplot::plot_grid(Fig3abd, Fig3e, rel_heights = c(2,1), ncol=1, label_size = 12, align = 'v', axis='l') 

Fig3 <- cowplot::plot_grid(Fig3abde, Fig3c, rel_widths = c(2,1), nrow=1, label_size = 12) 

ggsave("Plots/main_figures/Fig3.png", Fig3, width = 7.5, height = 7)



##### Fig. 4 Ancient balancing selection for environmental response genes #####

### ortholog identity plot ###

plot_ortho_id <- df_ortho_id %>%
  ggplot(.) +
  aes(x=pair, y=identity, fill=class) +
  #geom_point() +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0) +
  scale_fill_brewer(palette = 'Set1') +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.y= element_text(size=10, color='black'), 
        axis.text.x = element_text(size=10, color='black'),
        legend.text = element_text(size=10, color='black'),
        panel.grid = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.position = c(0.6,0.2)) +
  labs(y="Identity (%)", fill="") +
  guides(fill= guide_legend(nrow=2)) +
  coord_cartesian(ylim=c(80,100))

plot_ortho_id

ggsave("Plots/main_figures/Fig4c.pdf", plot_ortho_id, width = 3, height = 3.5)


### Fig. 4c overlay of divergent regions with long-read alignments ###

CB4856_divergent <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  dplyr::filter(STRAIN == "CB4856")

plot_CB4856_divcomp_II <- df_pacbio %>%
  dplyr::filter(CHROM %in% c("II"), N2_length > 5e3, N2_start > 1e6, N2_end < 4e6, CB4856_start > 1e6, CB4856_end < 4e6) %>%
  ggplot(.) +
  geom_rect(aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymin=CB4856_start/1e6, ymax=CB4856_end/1e6), 
            fill = 'red', color='red', alpha=0.8) +
  geom_rect(data=dplyr::filter(CB4856_divergent,CHROM %in% c("II"), cluster_start > 1e6, cluster_end < 4e6), aes(xmin=cluster_start/1e6, xmax=cluster_end/1e6, ymin=1, ymax=4), 
            fill = 'grey', alpha=0.4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title=element_text(size=11, color='black'), 
        axis.text=element_text(size=10, color='black'),
        plot.margin = unit(c(t=0.1,r=0,b=0.1,l=0.1), "in")) +
  facet_wrap(~CHROM, nrow=1) +
  labs(x="N2 genomic position (MB)", y="CB4856 genomic position (MB)") +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(1,4)) +
  scale_y_continuous(expand = c(0.01,0.01), limits = c(1,4))

plot_CB4856_divcomp_V <- df_pacbio %>%
  dplyr::filter(CHROM %in% c("V"), N2_length > 5e3, N2_start > 15e6, N2_end < 20e6, CB4856_start > 15e6, CB4856_end < 20e6) %>%
  ggplot(.) +
  geom_rect(aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymin=CB4856_start/1e6, ymax=CB4856_end/1e6), 
            fill = 'red', color='red', alpha=0.8) +
  geom_rect(data=dplyr::filter(CB4856_divergent,CHROM %in% c("V"), cluster_start > 15e6, cluster_end < 20e6), aes(xmin=cluster_start/1e6, xmax=cluster_end/1e6, ymin=15, ymax=20), 
            fill = 'grey', alpha=0.4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x=element_text(size=11, color='black'), 
        axis.title.y=element_blank(),
        axis.text=element_text(size=10, color='black'),
        plot.margin = unit(c(0.1,0.1,0.1,0), "in")) +
  facet_wrap(~CHROM, nrow=1) +
  labs(x="N2 genomic position (MB)", y="CB4856 genomic position (MB)") +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(15,20)) +
  scale_y_continuous(expand = c(0.01,0.01), limits = c(15,20))

CB4856_divcomp <- cowplot::plot_grid(plot_CB4856_divcomp_II, plot_CB4856_divcomp_V, nrow=1, 
                                     align = 'v', rel_heights = c(1,0.95), axis = "l")

CB4856_divcomp

### Fig. 4d ortholog comparisons across species ###

plot_bl <- df_bl %>%
  ggplot(.) +
  aes(x=pair, y=branch_length, fill=class) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0) +
  scale_fill_brewer(palette = 'Set1') +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.y= element_text(size=10, color='black'), 
        axis.text.x = element_text(size=9, color='black'),
        legend.text = element_text(size=9, color='black'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.38,0.88),
        legend.background = element_rect(fill='transparent'),
        plot.margin = unit(c(0.1,0.1,0.1,0,1), "in")) +
  labs(y="Branch length") +
  guides(fill= guide_legend(nrow=2)) +
  coord_cartesian(ylim=c(0,0.18))

plot_bl

Fig4ab <- cowplot::plot_grid(plot_div_TajD, plot_wormcat_2, nrow=1, rel_widths = c(1,4), labels = c("a","b"), label_size = 12) 
Fig4cd <- cowplot::plot_grid(CB4856_divcomp, plot_bl, nrow=1, rel_widths = c(3,2), labels = c("c","d"), label_size = 12) 

Fig4 <- cowplot::plot_grid(Fig4ab, Fig4cd, ncol=1, rel_heights = c(1,1), label_size = 12) 

Fig4

ggsave("Plots/main_figures/Fig4.png", Fig4, width = 7.5, height = 5.5)


########## Extended Data Figures ############

##### Extended Fig. 1 Admixture analysis ######

### Ext Fig. 1a - Admixture composition of 328 wild isotypes at K=16 ###

## load admixture 

# establish plot order of strains based on anc pop and max fraction

plot_order <- long_admix_pops %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))

admix_plots_K16 <- long_admix_pops %>%
  dplyr::mutate(ordered_samples = factor(samples, levels = plot_order$samples)) %>%
  ggplot() +
  geom_bar(stat = "identity", 
           aes(x = ordered_samples, 
               y = frac_cluster, 
               fill = cluster)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(fill = "", x = "") +
  guides(fill= guide_legend(nrow=1)) +
  theme_bw() +
  theme(axis.text.x=element_blank(),    
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(t=0,r=0.1,b=0.2,l=0.2), units = "in"),
        legend.position = c(0.5,0.03),
        legend.key.size = unit(0.03,"in"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.text = element_text(size=10, color='black'))

admix_plots_K16

#ggsave("Plots/admix_plots_K16.pdf", admix_plots_K16, width = 7.5, height = 2.5)

### Ext Fig. 1b Histogram of major admixture fraction ###

plot_admix_max_freq <-df_admix_max_freq %>%
  ggplot(.) +
  geom_histogram(binwidth=0.05) +
  aes(x=max_frac) +
  theme_bw() +
  theme(axis.text = element_text(size=11, color='black'), axis.title = element_text(size=12, color='black'), panel.grid = element_blank()) +
  labs(x="Major ancestry fraction", y = "Isotype count")

plot_admix_max_freq

#ggsave("Plots/plot_admix_max_freq.pdf", plot = plot_admix_max_freq , width = 7.5, height = 4, units = "in")

### Ext Fig. 1c admixture fraction ###

plot_admix_fraction <- long_admix_pops %>%
  dplyr::filter(frac_cluster==max_frac) %>%
  ggplot(.) +
  aes(x=cluster, y=max_frac, fill=cluster) +
  geom_point(alpha = 0.7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.8) +
  #geom_violin(alpha = 0.5, bw =0.07) +
  scale_fill_manual(values = ancestry.colours) +
  theme_bw() +
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.text = element_text(size=10, color='black'), axis.title = element_text(size=11, color='black'), panel.grid = element_blank()) +
  labs(y= "Major ancestry\nfraction") +
  ylim(0, 1)

plot_admix_fraction 

#ggsave("Plots/plot_admix_fraction.pdf", plot_admix_fraction, width = 7.5, height = 2.5)

ExFig1  <- cowplot::plot_grid(admix_plots_K16, plot_admix_max_freq, plot_admix_fraction, ncol=1, rel_heights = c(4,4,3), labels = c("a","b","c"), label_size = 12)

ExFig1

ggsave("Plots/extended_data_figures/ExtFig1.png", ExFig1, width = 7.5, height = 7, units = "in")

##### Extended Fig. 2 Treemix analysis (m=9) ######

### Ext Fig. 2a, b Treemix tree and residuals (m=9) ###

# load functions for ploting TreeMix output
source("Scripts/PLOT_TREEMIX_DL.R")

plot_tree(stem = "Processed_Data/treemix/K-16_Outgroup=PA=TREEMIX_input_m2", cex=0.8)
tree_mix_m2 <- ggdraw(recordPlot())

plot_resid(stem = glue::glue("Processed_Data/treemix/K-16_Outgroup=PA=TREEMIX_input_m2"), pop_order = "Processed_Data/poporder", cex=0.7)
tree_mix_m2_resid <- ggdraw(recordPlot())

####### Ext Fig. 2c Haplotype map for the unadmixed isotypes

df_unadmix_sweep <- df_unadmix %>%
  dplyr::select(isotype, cluster) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::distinct(isotype, cluster, chromosome, filtered_sweep_ratio)

df_unadmix_summary <- df_unadmix %>%
  dplyr::select(isotype, cluster) %>%
  dplyr::arrange(as.character(cluster)) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
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

mcolor_grp <- plot_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
mcolor <- mcolor_grp$color
names(mcolor) <- mcolor_grp$haplotype

plot_df_unadmix <- plot_df %>%
  dplyr::filter(isotype %in% unadmix_isotype) %>%
  dplyr::left_join(., df_unadmix_summary, by = 'isotype')

plot_unadmix_hap <- plot_df_unadmix %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_unadmix$isotype)),
                     labels = df_unadmix_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Genomic position (Mb)") +
  theme_bw() +
  facet_grid(cluster~chromosome, scales="free", space="free", switch = "y") +
  theme(axis.text.x = element_blank(), axis.text.y  = element_blank(), strip.text.x = element_text(size=8, color = "black", margin = margin(0.03,0.03,0.03,0.03, "in")), legend.position = 'none', axis.title.x = element_text(size=10, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=8, color = "black", angle = 180, margin = margin(0.03,0.03,0.03,0.03, "in")),panel.spacing = unit(0.1, "lines")) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20))

plot_unadmix_hap

ExFig2a <- cowplot::plot_grid(tree_mix_m2, labels = "a", label_size = 12)
ExFig2b <- cowplot::plot_grid(tree_mix_m2_resid, labels = "b", label_size = 12)

ExFig2c <- cowplot::plot_grid(plot_unadmix_hap, labels = "c")

ExFig2  <- cowplot::plot_grid(NULL, ExFig2c, nrow=1, label_size = 12)

ggsave("Plots/ExtFig2a.pdf", ExFig2a, width = 3.75, height = 2.8, units = "in")
ggsave("Plots/ExtFig2b.pdf", ExFig2b, width = 3.75, height = 2.2, units = "in")
ggsave("Plots/ExtFig2.pdf", ExFig2, width = 7.5, height = 5, units = "in")


##### Extended Fig.3 Selective sweeps ##### 

strain_labels <- plot_df %>%
  dplyr::ungroup() %>%
  dplyr::select(plotpoint, isotype) %>%
  dplyr::distinct() %>%
  dplyr::arrange(plotpoint)

df_sweep_all <- df_subpop %>%
  dplyr::distinct(isotype, cluster) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::distinct(isotype, cluster, chromosome, filtered_sweep_ratio)

df_sweep_all_summary <- df_subpop %>%
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

df_unadmix_sweep <- df_unadmix %>%
  dplyr::select(isotype, cluster) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::distinct(isotype, cluster, chromosome, filtered_sweep_ratio)

plot_df_sweep_all <- plot_df %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  dplyr::left_join(., df_sweep_all_summary, by = 'isotype')

## sweep in all strains

plot_sweep_all_hap_hor <- plot_df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
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

plot_sweep_all_hor <- df_sweep_all %>%
  dplyr::filter(!chromosome %in% c("II","III")) %>%
  ggplot(.) +
  #geom_jitter(alpha=0.8, size = 0.6) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.6) +
  aes(x=cluster, y= filtered_sweep_ratio*100, fill=cluster) +
  scale_fill_manual(values = ancestry.colours) +
  labs(x="Ancestry", y="Proportion swept (%)")+
  scale_y_continuous(breaks = c(0,20, 40, 60, 80, 100))+
  scale_x_discrete(limits = rev(levels(df_unadmix_sweep$cluster)))+
  coord_flip()+
  theme_bw() +
  theme(text = element_text(size=11, color='black'), legend.position = 'none', axis.title.y = element_blank(), strip.background = element_blank(), strip.text = element_blank(),  panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in")) +
  facet_wrap(~chromosome, nrow=1)

plot_sweep_all_hor

plot_sweep_all <- cowplot::plot_grid(plot_sweep_all_hap_hor, plot_sweep_all_hor, 
                                     ncol=1, rel_heights = c(5,2), align = "h", axis = "lr",
                                     labels = c("a", "b"), label_size = 12)

plot_sweep_all

ggsave(plot_sweep_all, file = "Plots/extended_data_figures/ExtFig3.png", 
       width = 7.5,
       height = 9, unit = "in")



##### Extended Fig.5 Characterization of hyper-divergent regions ##### 

### Ext Fig. 5a Example shortread alignment in divergent regions ###

### Ext Fig. 5b ID-divergent-region pipeline flow chart ###

### Ext Fig. 5c Divergent region example - classification, validation with CB4856 long-read ###

plot_CB4856_div_example <- div_longread_freq(df1=df_pacbio, df2_1=div_CB4856_all_select, df2_2=df_CB4856_freq, chr="V",left=520000, right=589000, ymin=65, guide_rownum = 2)

ExFig5a <- cowplot::plot_grid(NULL, labels = c("a"), label_size = 12)

ExFig5bc <- cowplot::plot_grid(NULL, plot_CB4856_div_example, nrow=1, rel_widths = c(6,5), labels = c("b","c"), label_size = 12)

ExFig5 <- cowplot::plot_grid(ExFig5a, ExFig5bc, ncol=1, rel_heights = c(3,2), label_size = 12)

ggsave("Plots/ExtFig5.pdf", ExFig5, width = 7.5, height = 6.5, units = "in")



##### Extended Fig.6 C. briggsae divergent region ##### 

cbr_df_common_div <- cbr_df_all_cluster_size_freq %>%
  dplyr::filter(freq >= 0.3) %>%
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
  labs(x="Genomic position (Mb)", y="36 wild C. briggsae strains") +
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
        strip.text.y = element_text(size=9, color = "black", angle = 180)) +
  labs(x="Genomic position (Mb)", y="") +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(group~CHROM, scales="free", switch = "y", space = 'free')

cbr_plot_div_pop_rug

ExFig6 <- cowplot::plot_grid(cbr_plot_div_pop, cbr_plot_div_pop_rug, ncol=1, rel_heights = c(10,1), align = "v")

ExFig6

ggsave("Plots/extended_data_figures/ExtFig6.pdf", ExFig6, width =  7.5, height = 6)

```

```{r setup, warning=FALSE, message =FALSE}

########## Supplementary figures ############

##### Fig. S1 Admixture analysis with various K  ######

plot_admix_cv <- admix_cv %>%
  ggplot(.) +
  geom_jitter(size = 0.3, alpha = 0.7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
  aes(x=K, y=CV, group = K)+
  theme_bw() +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank())

admix_plots <- df_admix_full %>%
  dplyr::left_join(., na.omit(dplyr::select(df_admix_K16, samples, K16_group))) %>%
  #dplyr::filter(K %in% c(14:18)) %>%
  dplyr::mutate(ordered_samples = factor(samples, levels = plot_order_admix$samples), K = paste("K=",K, sep="")) %>%
  ggplot() +
  geom_bar(stat = "identity", 
           aes(x = ordered_samples, 
               y = frac_cluster, 
               fill = cluster)) +
  scale_fill_manual(values = K20.colours) +
  labs(fill = "", x = "", y="") +
  theme_bw() +
  theme(axis.text.x=element_blank(),    
        axis.text.y=element_blank(),
        axis.title.y = element_text(angle = 0, vjust = .5),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        legend.position = 'none',
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 9, color='black'),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")) +
  facet_grid(K~K16_group, scales = 'free', space = "free_x",  switch = "y")

admix_plots

FigS1 <- cowplot::plot_grid(plot_admix_cv, admix_plots,ncol=1, rel_heights = c(1,3), labels = c("a","b"))

ggsave(FigS1, file = "Plots/supplementary_figures/FigS1.png", width=7.5, height=9, unit = 'in')

##### Fig. S2 Geographic distribution of unadmixed isotypes  ######

### Fig. S2a Distribution of unadmixed isotypes in Europe  ####

EU_isotypes <- as.character(dplyr::filter(indep_strain_info, long > -40 & long < 25 & lat > 35 & lat <57)$isotype)

admix_strain_info_EU <- admix_strain_info %>%
  dplyr::filter(isotype %in% EU_isotypes)

plot_map_EU_unadmix <- 
  ggplot()+ geom_map(data=world, map=world,
                     aes(x=long, y=lat, map_id=region),
                     color="gray51", fill="white", size=0.2, alpha = 0.8)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size = 2, alpha = 0.7, fill = "lightgray") +
  geom_label_repel(aes(long, lat, label = strain, fill = cluster), data = admix_strain_info_EU, fontface = 'bold', color = 'white', size = 0, box.padding = 0.1, point.padding = 0.2, segment.color = '#9acd32', segment.size = 0.4, segment.alpha = 0.7) +
  #ggsn::scalebar(data = admix_strain_info, dist = 500, dist_unit = "km", location = "bottomleft", transform = TRUE, model = "WGS84", height = 0.02, st.size = 2.3, st.dist	= 0.03, border.size = 0.3) +
  ggsn::scalebar(data = admix_strain_info,  dist = 200, dist_unit = "km", location = "bottomleft",
                 transform = TRUE, model = "WGS84", height = 0.004, border.size = 0.3, st.size = 3.5, st.dist	= 0.005, anchor = c(x = -27, y = 56)) +
  scale_fill_manual(values = ancestry.colours) +
  guides(fill=guide_legend(ncol=3, nrow=3)) +
  theme_map()+
  theme(legend.text = element_text(size=11, color='black'), legend.position = c(0.05, 0.6), legend.direction = "horizontal", legend.title = element_blank(), legend.key.size = unit(0.05, "in"),  legend.key.width = unit(0.05,"in") ) +
  lims(x=c(-29, 25), y=c(36, 57)) 

plot_map_EU_unadmix

### Fig. S2b Pairwise distances among unadmixed isotypes in the same group  ####

df_dist_admix <- df_dist_in_bw_iso %>%
  dplyr::filter(isotype %in% unadmix_isotype & isotype_comp %in% unadmix_isotype) %>%
  dplyr::filter(class == "between_isotype") %>%
  dplyr::left_join(., df_isotype_admix, by='isotype') %>%
  dplyr::left_join(., dplyr::rename(df_isotype_admix, isotype_comp=isotype), by='isotype_comp') %>%
  dplyr::rename(cluster=cluster.x, cluster_comp=cluster.y) %>%
  dplyr::filter(cluster == cluster_comp)

plot_dist_unadmix <- df_dist_admix %>%
  ggplot(.) +
  aes(x=cluster, y=dist, fill = cluster) +
  geom_jitter(size = 0.5, alpha = 0.8) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.8) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.x = element_text(size=10, color='black'), 
        axis.text.y = element_text(size=10, color='black'), 
        legend.position = 'none') +
  labs(y="Distance (km)")

FigS2 <- cowplot::plot_grid(plot_map_EU_unadmix, plot_dist_unadmix,ncol=1, rel_heights = c(5,3), labels = c("a","b"))

ggsave(FigS2, file = "Plots/supplementary_figures/FigS2.pdf", width=7.5, height=7, unit = 'in')

##### Fig. S3 Treemix analysis  ######

### Fig. S3a-d Treemix outputs with m=0,3,6,9 ###

### Fig. S3e Variance explained by the population tree model with different migration numbers  ####

plot_tree_ve <- df_treesum_all %>%
  ggplot(.) +
  geom_col() +
  aes(x=migration, y=VarExplain*100) +
  theme_bw() +
  coord_cartesian(ylim=c(85,100)) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) +
  theme(axis.title.x = element_text(size=11, color='black'), 
        axis.text.x = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.y = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  labs(x="Number of migration", y="Variance explained (%)")

plot_tree_ve

FigS3d <- cowplot::plot_grid(plot_tree_ve, labels = "d", label_size = 12)

FigS3 <- cowplot::plot_grid(NULL, FigS4d, ncol=1, rel_heights = c(2.2,1))

ggsave("Plots/supplementary_figures/FigS3.pdf", FigS3, width = 7.5, height = 8, units = "in")


### Fig. S6 Haplotype plots for all 328 wild isotpes ###

plot_all_hap <- plot_df_sweep_all %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=11, color = "black"), axis.text.y  = element_blank(), strip.text = element_text(size=11, color = "black"), legend.position = 'none', axis.title.x = element_text(size=13, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in")) +
  facet_grid(cluster~chromosome, scales="free", space="free_y", switch = "y")

plot_all_hap

ggsave("Plots/supplementary_figures/plot_all_hap.pdf", plot_all_hap, width = 7.5, height = 8, units = "in")


##### Fig. S7 haplotype plot for WS2 and CA group  ######

df_Cali_summary <- df_subpop %>%
  dplyr::filter(state=="California") %>%
  #dplyr::filter(cluster %in% c("WS2", "CA")) %>%
  dplyr::arrange(as.character(cluster)) %>%
  dplyr::left_join(., plot_df, by='isotype') %>%
  dplyr::distinct(isotype, cluster, chromosome, isotype_swept_haplotype_length, state) %>%
  dplyr::group_by(isotype, cluster,state) %>%
  dplyr::summarise(sum_sweep = sum(isotype_swept_haplotype_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(isotype, cluster, sum_sweep, state) %>%
  dplyr::group_by(cluster,state) %>%
  dplyr::mutate(mean_sweep = mean(sum_sweep)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(mean_sweep, sum_sweep)) %>%
  dplyr::mutate(plotpoint=row_number())

plot_hap_Cali <- plot_df %>%
  dplyr::left_join(., df_Cali_summary, by = 'isotype') %>%
  dplyr::filter(state=="California") %>%
  #dplyr::filter(cluster %in% c("WS2", "CA")) %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  xlab("Genomic position (Mb)") +
  theme_bw() +
  facet_grid(cluster~chromosome, scales="free", space="free", switch = "y") +
  theme(axis.text.x = element_blank(), axis.text.y  = element_blank(), strip.text.x = element_text(size=10, color = "black", margin = margin(0.03,0.03,0.03,0.03, "in")), legend.position = 'none', axis.title.x = element_text(size=11, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=8, color = "black", angle = 180, margin = margin(0.03,0.03,0.03,0.03, "in")),panel.spacing = unit(0.1, "lines")) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  scale_y_continuous(expand = c(0, 0))

plot_hap_Cali

ggsave("Plots/supplementary_figures/FigS4.pdf", plot_hap_Cali, width = 7.5, height = 5.5, units = "in")



### statistics for arms vs centers

##### Fig. S7 Divergent region classification parameter search using CB4856 long-read aseemblies #####
plot_CB_div_id <- df_div_comp_parameter_gather_ID %>%
  dplyr::filter(count_threshold >= 5 & count_threshold <= 20) %>%
  ggplot(.) +
  #geom_line() +
  geom_point(aes(fill=count_threshold, shape = class), alpha = 0.7, size = 1.5) +
  aes(x=cluster_threshold/1e3, y=identity) +
  #scale_size_continuous(range = c(0.5,3), breaks = seq(0,20, by=5)) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  scale_x_continuous(expand = c(0.1, 0.1), breaks = c(10,20,30,40)) +
  scale_shape_manual(values=c(22, 21), guide=FALSE) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size=11, color='black'), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'),
        legend.text = element_text(size=10), 
        panel.grid = element_blank(),
        panel.spacing.y = unit(0.7, "lines")) +
  labs(x="Minimum cluster size (kb)", y="Long-read alignemnt identity (%)", fill = "Variant count\nthreshold") +
  facet_grid(class~cov_threshold, scales = 'free')

plot_CB_div_id

plot_CB_div_coverage <- df_div_comp_parameter_gather_coverage %>%
  dplyr::filter(count_threshold >= 5 & count_threshold <= 20) %>%
  ggplot(.) +
  #geom_line() +
  geom_point(aes(fill=count_threshold, shape = class), alpha = 0.7, size = 1.5) +
  aes(x=cluster_threshold/1e3, y=coverage*100) +
  #scale_size_continuous(range = c(0.5,3), breaks = seq(0,20, by=5)) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  scale_x_continuous(expand = c(0.1, 0.1), breaks = c(10,20,30,40)) +
  scale_shape_manual(values=c(22, 21), guide=FALSE) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size=11, color='black'), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'),
        legend.text = element_text(size=10), 
        panel.grid = element_blank(),
        panel.spacing.y = unit(0.7, "lines")) +
  labs(x="Minimum cluster size (kb)", y="Long-read alignemnt coverage (%)", fill = "Variant count\nthreshold") +
  facet_grid(class~cov_threshold, scales = 'free')

plot_CB_div_coverage

plot_CB_div_id_comp_size <- df_div_comp_parameter_gather_ID%>%
  dplyr::filter(count_threshold >= 5 & count_threshold <= 20) %>%
  dplyr::filter(class == "Hyper-divergent") %>%
  ggplot(.) +
  #geom_line() +
  geom_point(aes(fill=count_threshold), alpha = 0.7, shape =21, size = 1.5) +
  aes(x=cluster_threshold/1e3, y=Divergent_size/1e6) +
  #scale_size_continuous(range = c(0.5,3), breaks = seq(0,20, by=5)) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  scale_x_continuous(expand = c(0.1, 0.1), breaks = c(10,20,30,40)) +
  facet_grid(class~cov_threshold) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size=11, color='black'), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'),
        legend.text = element_text(size=10), 
        panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x=element_blank()) +
  labs(x="Minimum cluster size (kb)", y="Size (Mb)", fill = "Variant count\nthreshold")

plot_CB_div_id_comp_size

### N2 false positives ###

df_div_fp_parameter <- read_tsv("Processed_Data/df_div_fp_parameter.tsv")

plot_div_fp_parameter <- df_div_fp_parameter %>%
  dplyr::filter(!is.na(cov_threshold) & cov_threshold <= 0.3) %>%
  dplyr::mutate(label="False-positive in N2") %>%
  ggplot(.) +
  geom_point(aes(fill=count_threshold), alpha = 0.7, shape =21, size = 1.5) +
  aes(x=cluster_threshold/1e3, y=Divergent_clustered_size/1e6) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'), 
        axis.title.x = element_text(size=11, color='black'), 
        axis.text.y = element_text(size=10, color='black'),
        legend.text = element_text(size=10), 
        panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x=element_blank()) +
  facet_grid(label~cov_threshold, scales = 'free') +
  scale_x_continuous(expand = c(0.1, 0.1), breaks = c(10,20,30,40)) +
  labs(x="Minimum cluster size (kb)", y="False detection (Mb)", fill = "Variant count\nthreshold")

plot_div_fp_parameter

FigS7 <- cowplot::plot_grid(plot_CB_div_id, plot_CB_div_coverage, plot_CB_div_id_comp_size, plot_div_fp_parameter, 
                            labels = c("a","b","c","d"), ncol=1, 
                            align = 'v', rel_heights = c(2,2,1,1.3), label_size = 12)

FigS7 

ggsave(FigS7, file = "Plots/supplementary_figures/FigS7.pdf", width=7.5, height=11, unit = 'in')

##### Fig. S8 Validation of CB4856 divergent regions using CB4856 long-read aseemblies #####

CB4856_divergent <- df_divergent_final %>%
  dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end) %>%
  dplyr::filter(STRAIN == "CB4856")

plot_list_div_CB4856 <- list()

for (i in 1:nrow(CB4856_divergent)) {
  chr <- as.character(CB4856_divergent$CHROM[i])
  left <- as.numeric(CB4856_divergent$cluster_start[i])
  right <- as.numeric(CB4856_divergent$cluster_end[i])
  plot_list_div_CB4856[[i]] <- div_longread_freq(chr=chr, left=left, right=right, guide_rownum = 2)
}

plot_grid_div_CB4856_1 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[1:8], ncol=2)
plot_grid_div_CB4856_2 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[9:16], ncol=2)
plot_grid_div_CB4856_3 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[17:24], ncol=2)
plot_grid_div_CB4856_4 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[25:32], ncol=2)
plot_grid_div_CB4856_5 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[33:40], ncol=2)
plot_grid_div_CB4856_6 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[41:48], ncol=2)
plot_grid_div_CB4856_7 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[49:56], ncol=2)
plot_grid_div_CB4856_8 <- cowplot::plot_grid(plotlist=plot_list_div_CB4856[57:62], ncol=2)

ggsave(plot_grid_div_CB4856_1, file = "Plots/supplementary_figures/FigS8_1.pdf", width=8.5, height=11, unit = 'in')
ggsave(plot_grid_div_CB4856_2, file = "Plots/supplementary_figures/FigS8_2.pdf", width=8.5, height=11, unit = 'in')
ggsave(plot_grid_div_CB4856_3, file = "Plots/supplementary_figures/FigS8_3.pdf", width=8.5, height=11, unit = 'in')
ggsave(plot_grid_div_CB4856_4, file = "Plots/supplementary_figures/FigS8_4.pdf", width=8.5, height=11, unit = 'in')
ggsave(plot_grid_div_CB4856_5, file = "Plots/supplementary_figures/FigS8_5.pdf", width=8.5, height=11, unit = 'in')
ggsave(plot_grid_div_CB4856_6, file = "Plots/supplementary_figures/FigS8_6.pdf", width=8.5, height=11, unit = 'in')
ggsave(plot_grid_div_CB4856_7, file = "Plots/supplementary_figures/FigS8_7.pdf", width=8.5, height=11, unit = 'in')
ggsave(plot_grid_div_CB4856_8, file = "Plots/supplementary_figures/FigS8_8.pdf", width=8.5, height=8.5, unit = 'in')


##### Fig. S9 Size and frequencies of joined hyper-divergent blocks #####

plot_all_cluster_size_freq <- df_all_cluster_size_freq %>%
  ggplot(.) +
  geom_point(aes(x=cluster_size/1e3, y=Frequency), size = 0.6, alpha =0.8) +
  #geom_text_repel(data=dplyr::filter(df_div_common_fraction, frac_freq_strain<0.5), aes(x=total_div_strain/1e3, y=frac_freq_strain*100, label = STRAIN), size = 2.5, box.padding = 0.2, point.padding = 0.1, segment.size = 0.02, nudge_y = 0.1, nudge_x= 0, segment.alpha = 0) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color='black'), 
        axis.text = element_text(size=10, color='black'), 
        panel.grid=element_blank()) +
  labs(x="Size (kb)", y="Average frequencies of bins")

plot_all_cluster_size_freq

ggsave(plot_all_cluster_size_freq, file = "Plots/supplementary_figures/FigS14.pdf", width=7.5, height=5, unit = 'in')

##### Fig. S10 Summary of statistics for hyper-divergent regions across six chromosomes #####

plot_div_chrom_comp <- df_div_stats_genome_summary_gather %>%
  dplyr::filter(stats != "Size" & class != "All") %>%
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

plot_fold_diff_div <- df_div_stats_genome_summary_gather_spread %>%
  dplyr::filter(stats != "Size") %>%
  dplyr::rename(Non_divergent = "Non-divergent") %>%
  ggplot(.) +
  geom_col() +
  aes(x=location, y=as.numeric(Divergent)/as.numeric(Non_divergent)) +
  theme_bw() +
  #scale_fill_brewer(palette = 'Set1') +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1, color='black'),
        axis.title.y = element_text(size=11, color='black'), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=10, color='black'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=10, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.direction = "horizontal", legend.box = "horizontal")+
  labs(y="Fold difference") +
  facet_grid(stats~CHROM, scales = 'free') 

FigS9 <- cowplot::plot_grid(plot_div_chrom_comp, plot_fold_diff_div, 
                            labels = c("a","b"), ncol=1, 
                            align = 'v', rel_heights = c(1.2,1), axis = 'l', label_size = 12)

FigS9

ggsave(FigS9, file = "Plots/supplementary_figures/FigS9.pdf", width=7.5, height=9, unit = 'in')


##### Fig. S11 Population structure x hyper-divergent regions across the genome #####

plot_div_pop_K16 <- plot_df_sweep_all %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  ggplot(.) +
  geom_rect(data=dplyr::rename(df_join_div_sweep, chromosome=CHROM), 
            aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = plotpoint-0.5 , ymax = plotpoint+0.5), 
            fill='black', color = 'black', size = 0.1, alpha = 0.5) +
  theme_bw() +
  xlab("Genomic position (Mb)") +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),  
        axis.text.y  = element_blank(), strip.text = element_text(size=10, color = "black"),
        legend.position = 'none', 
        axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180), 
        panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in"),
        panel.grid = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  facet_grid(cluster~chromosome, scales="free", space="free_y", switch = "y")

plot_div_pop_K16

plot_div_all_hor <- df_div_fraction_chrom_fill %>%
  ggplot(.) +
  geom_jitter(alpha=0.7, size = 0.2, width = 0.2) +
  geom_boxplot(alpha = 0.8, outlier.alpha=0) +
  aes(x=cluster, y= div_fraction*100, fill=cluster) +
  scale_fill_manual(values = ancestry.colours) +
  labs(x="Ancestry", y="Proportion hyper-divergent (%)")+
  scale_x_discrete(limits = rev(levels(df_unadmix_sweep$cluster)))+
  coord_flip()+
  theme_bw() +
  theme(axis.text.x = element_text(size=10, color='black'), 
        axis.text.y = element_text(size=10, color='black'), 
        axis.title.x = element_text(size=11, color='black'), 
        legend.position = 'none', axis.title.y = element_blank(), 
        strip.background.x = element_blank(), strip.text.x = element_blank(),  
        panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in")) +
  facet_wrap(~CHROM, nrow=1, scales = 'free')

plot_div_all_hor

plot_div_all <- cowplot::plot_grid(plot_div_pop_K16, plot_div_all_hor, ncol=1, rel_heights = c(5,2.1), labels = c("a", "b","c"), label_size = 12)

ggsave(plot_div_all, file = "Plots/supplementary_figures/FigS11.pdf", 
       width = 7.5,
       height = 9.5, unit = "in")

##### Fig. S12 fraction of rare divergent regions for 16 admixture groups #####

df_div_rare_fraction <- df_divergent_final %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(total_div_strain = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN, freq_bin) %>%
  dplyr::mutate(n_bin_strain=n(), frac_freq_strain=n()/total_div_strain) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN, freq_bin, total_div_strain, frac_freq_strain) %>%
  dplyr::filter(freq_bin == "Rare")

df_div_rare_fraction_admix <- df_div_rare_fraction %>%
  dplyr::rename(isotype = STRAIN) %>%
  dplyr::left_join(., df_isotype_admix, by ='isotype') 

df_div_rare_fraction_admix$cluster <- factor(df_div_rare_fraction_admix$cluster, levels=c("WS1","WS2","EU1", "EU2", "FR","AU1","AU2","AT","PT1","PT2","WD","CA","HW1","HW2","HW3","PA"))

plot_div_rare_fraction_admix <- df_div_rare_fraction_admix %>%
  ggplot(.) +
  geom_point(alpha = 0.7, position=position_jitterdodge(jitter.width = 3), size = 1, aes(x=cluster, y=frac_freq_strain*100, fill = cluster)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.8, aes(x=cluster, y=frac_freq_strain*100, fill = cluster)) +
  geom_text_repel(data=dplyr::filter(na.omit(df_div_rare_fraction_admix), frac_freq_strain > 0.1), aes(x=cluster, y=frac_freq_strain*100, group=cluster, label = isotype), size = 2.5, box.padding = 0.2, point.padding = 0.15, segment.size = 0.02, nudge_x = 0, nudge_y= -0.02, segment.alpha = 0) +
  scale_fill_manual(values = ancestry.colours) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11, color='black'), axis.title.x = element_blank(), axis.text.y = element_text(size=10, color='black'), legend.position = 'none', axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.6, hjust=1),plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in"), panel.grid = element_blank()) +
  labs(x="Subgroup", y="Rare divergent regions (%)")

plot_div_rare_fraction_admix

ggsave(plot_div_rare_fraction_admix, file = "Plots/supplementary_figures/FigS12_rare.png", 
       width = 6,
       height = 4, unit = "in")

##### Fig. S13 Population structure x hyper-divergent regions x sweeps across the genome for swept populations #####

plot_div_sweep_pop <- plot_df_sweep_all %>%
  dplyr::rename(plotpoint=plotpoint.y) %>%
  dplyr::filter(cluster %in% swept_clusters) %>%
  ggplot(.) +
  geom_rect(aes(xmin = start/1E6, xmax = stop/1E6,ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,fill = swept_haplotype)) +
  geom_rect(data=dplyr::rename(dplyr::filter(df_join_div_sweep,cluster %in% swept_clusters), chromosome=CHROM), aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = plotpoint-0.5 , ymax = plotpoint+0.5), fill='black', color = 'black', size = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = 1:length(unique(plot_df_sweep_all$isotype)),
                     labels = df_sweep_all_summary$cluster,
                     expand = c(0, 0), position = 'left') +
  xlab("Genomic position (Mb)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),  
        axis.text.y  = element_blank(), strip.text = element_text(size=10, color = "black"),
        legend.position = 'none', 
        axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180), 
        panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0, l=0.1, r=0.1, t=0.1, unit = "in")) +
  facet_grid(cluster~chromosome, scales="free", space="free_y", switch = "y")

plot_div_sweep_pop

plot_N2_sweep <- plot_df_sweep_all %>%
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
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.title.x = element_text(size=11, color='black'),  
        legend.position = 'none', axis.title.y = element_blank(), 
        strip.background.x = element_blank(), strip.text.x = element_blank(),  
        panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in"),
        axis.ticks.y = element_blank(), 
        strip.text.y=element_text(size=9, color = "black", angle = 180)) +
  facet_grid(isotype~chromosome, scales="free", space="free_y", switch = "y")

plot_N2_sweep

df_div_fraction_chrom_fill2 <- df_div_fraction_chrom_fill %>%
  dplyr::filter(cluster %in% swept_clusters) 

df_div_fraction_chrom_fill2$cluster <- factor(df_div_fraction_chrom_fill2$cluster, levels = c("PT1", "AU2", "AU1", "FR", "EU2", "EU1", "WS2", "WS1"))

plot_div_all_hor_swept <- df_div_fraction_chrom_fill2 %>%
  ggplot(.) +
  geom_jitter(alpha=0.7, size = 0.2, width = 0.2) +
  geom_boxplot(alpha = 0.8, outlier.alpha=0) +
  aes(x=cluster, y= div_fraction*100, fill=cluster) +
  scale_fill_manual(values = ancestry.colours) +
  labs(x="Ancestry", y="Proportion hyper-divergent (%)")+
  #scale_x_discrete(limits = rev(levels(df_unadmix_sweep$cluster)))+
  coord_flip()+
  theme_bw() +
  theme(axis.text.x = element_text(size=10, color='black'), 
        axis.text.y = element_text(size=10, color='black'), 
        axis.title.x = element_text(size=11, color='black'), 
        legend.position = 'none', axis.title.y = element_blank(), 
        strip.background.x = element_blank(), strip.text.x = element_blank(),  
        panel.spacing = unit(0.1, "lines"), plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.1, unit = "in")) +
  scale_y_continuous(expand = c(0.02, 0.02), breaks = c(2,4,6,8)) +
  facet_wrap(~CHROM, nrow=1)

plot_div_all_hor_swept

plot_sweep_div_all <- cowplot::plot_grid(plot_div_sweep_pop, plot_N2_sweep, plot_div_all_hor_swept, 
                                         ncol=1, rel_heights = c(5,.45,2.1),
                                         align='v',
                                         labels = c("a", "b","c"), label_size = 12)

ggsave(plot_sweep_div_all, file = "Plots/supplementary_figures/FigS13.pdf", 
       width = 7.5,
       height = 9.5, unit = "in")


##### Fig. S12 Wormcat category 1 and 3 #####

plot_wormcat_1 <- df_cat1 %>%
  dplyr::left_join(., df_cat1_plotpoint, by=c("Category", "layer")) %>%
  ggplot(.) +
  geom_point(alpha=0.9) +
  aes(x=-log10(Bonferroni), y=plotpoint, shape=group, size=EnrichRatio, fill=RGS) +
  theme_bw() +
  #scale_x_continuous(breaks = c(2,4,6,8,10), expand = c(0.2,0.2)) +
  #scale_y_continuous(breaks = length(GO_list_MF):1, labels = c, name = "", expand = c(0.05,0.05)) +
  scale_y_continuous(breaks = 1:length(unique(df_cat1_plotpoint$Category)), labels = unique(df_cat1_plotpoint$Category)) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22)) +
  scale_size_continuous(range = c(1,4)) +
  theme(axis.text = element_text(size=10, color='black'), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.title = element_text(size=10, color='black'), 
        legend.text = element_text(size=9, color='black'), legend.margin = margin(0,0,0,0), 
        legend.direction = "horizontal", legend.box = "horizontal", legend.position = "bottom", 
        legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.2, "in"), 
        axis.title.y.right = element_text(size=10, color='black', vjust=2), panel.grid = element_blank()) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="-log10(Corrected P-value)", shape = "Frequency", size = "Fold enrichment", fill = "Gene counts") 

plot_wormcat_1

plot_wormcat_3 <- df_cat3 %>%
  dplyr::left_join(., df_cat3_plotpoint, by=c("Category", "layer")) %>%
  ggplot(.) +
  geom_point(alpha=0.9) +
  aes(x=-log10(Bonferroni), y=plotpoint, shape=group, size=EnrichRatio, fill=RGS) +
  theme_bw() +
  #scale_x_continuous(breaks = c(2,4,6,8,10), expand = c(0.2,0.2)) +
  #scale_y_continuous(breaks = length(GO_list_MF):1, labels = c, name = "", expand = c(0.05,0.05)) +
  scale_y_continuous(breaks = 1:length(unique(df_cat3_plotpoint$Category)), labels = unique(df_cat3_plotpoint$Category)) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22)) +
  scale_size_continuous(range = c(1,4)) +
  theme(axis.text = element_text(size=10, color='black'), axis.title.x = element_text(size=11, color='black'), 
        axis.title.y = element_blank(), legend.title = element_text(size=10, color='black'), 
        legend.text = element_text(size=9, color='black'), legend.margin = margin(0,0,0,0), 
        legend.direction = "horizontal", legend.box = "horizontal", legend.position = "bottom", 
        legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.2, "in"), 
        axis.title.y.right = element_text(size=10, color='black', vjust=2), panel.grid = element_blank()) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="-log10(Corrected P-value)", shape = "Frequency", size = "Fold enrichment", fill = "Gene counts") 

plot_wormcat_3

FigS12 <- cowplot::plot_grid(plot_wormcat_1, plot_wormcat_3, 
                             labels = c("a","b"), ncol=1, 
                             align = 'v', rel_heights = c(1.5,2), label_size = 12)

FigS12 

ggsave(FigS12, file = "Plots/supplementary_figures/FigS12.pdf", width=7.5, height=8, unit = 'in')


##### Fig. S13 pals loci in hyper-divergent regions #####

plot_div_pop_pals <- df_join_div %>%
  na.omit() %>%
  ggplot(.) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=1.01), color='transparent', fill='transparent', size =0.1) +
  geom_segment(data=dplyr::distinct(df_pals_div,div_start, .keep_all=T), aes(x=(start+stop)/2e6, xend=(start+stop)/2e6, y=0.4, yend=327.6), color='red', alpha=0.75, size=0.1) +
  #geom_segment(data=data.frame(CHROM="I", start=2342215, stop=2356238), aes(x=start/1e6, xend=stop/1e6, y=0.4, yend=327.6), color='blue', alpha=0.75, size=0.1) + ## peel-1, zeel-1
  #geom_segment(data=data.frame(CHROM="III", start=11116993, stop=11136380), aes(x=start/1e6, xend=stop/1e6, y=0.4, yend=327.6), color='purple', alpha=0.75, size=0.1) + ## sup-35, pha-1
  geom_rect(aes(xmin =  cluster_start/1e6, xmax = cluster_end/1e6, ymin = plotpoint-0.5 , ymax = plotpoint+0.5), fill = 'black',color='black', size = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y  = element_blank(), strip.text = element_text(size=10, color = "black"), legend.position = 'none',axis.title = element_text(size=11, color = "black"), axis.ticks.y = element_blank(), strip.text.y=element_text(size=9, color = "black", angle = 180),  panel.spacing = unit(0.1, "lines"), panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0.00, 0.00), limits=c(0.4,327.6)) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = c(5, 10, 15, 20)) +
  facet_grid(~CHROM, scales="free",space = 'free') +
  labs(x="Genomic position (Mb)",y="328 wild isotypes")

plot_div_pop_pals

ggsave(plot_div_pop_pals, file = "Plots/supplementary_figures/FigS13.png", width=7.5, height=5, unit = 'in')



##### Fig. S15 haplotype analysis at hyper-divergent regions #####

plot_div_hap_zoom_IIL <- divergent_hap_overlay_unadmix(df2=plot_df_unadmix_IIL_all, chr = "II", left=2234000, right=2438000, threshold = cluster_threshold)

plot_div_hap_zoom_IIIL <- divergent_hap_overlay_unadmix(df2=plot_df_unadmix_IIIL_all, chr = "III", left=111000, right=196000, threshold = cluster_threshold)

plot_div_hap_zoom_IVR <- divergent_hap_overlay_unadmix(df2=plot_df_unadmix_IVR_all, chr = "IV", left=17254000, right=17429000, threshold = cluster_threshold)

plot_div_hap_zoom_VR <- divergent_hap_overlay_unadmix(df2=plot_df_unadmix_VR_all, chr = "V", left=19238000, right=19514000, threshold = cluster_threshold)

FigS15 <- cowplot::plot_grid(plot_div_hap_zoom_IIL, plot_div_hap_zoom_IIIL, plot_div_hap_zoom_IVR, plot_div_hap_zoom_VR, labels = c("a", "b","c", "d"),ncol=2, align = 'hv', label_size = 12)

FigS15

ggsave("Plots/supplementary_figures/FigS15.pdf", FigS15, width = 7.5, height = 9.5)

```


``` {r}

####### Suplementary Data #######

##### DataS1. 609_wild_strains  #####

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


##### DataS2. Admixture components and classification of 328 wild isotypes #####

K=16
qfile <- pophelper::readQ(files = "Raw/LD_0.1_MAF_0.006.16.Q")[[1]]
colnames(qfile) <- LETTERS[1:K]
wide_admix_pops <- qfile %>%
  dplyr::mutate(isotype=strain_vector)
colnames(wide_admix_pops) <- c("WD", "WS1", "HW1", "EU2", "AU1", "CA", "PA", "PT2", "AU2", "HW3", "FR", "EU1", "AT", "HW2", "PT1", "WS2", "isotype")

df_admix_comp_all <- df_subpop %>%
  dplyr::distinct(isotype, cluster) %>%
  dplyr::left_join(., wide_admix_pops, by= 'isotype') %>%
  dplyr::rename(group=cluster)

write.csv(df_admix_comp_all, "Manuscript/Supplementary_Data/DataS2_admixture_components.csv")




##### DataS5. Variants count for div vs non-div in each isotype #####

DataS5 <- df_divergent_count_summary %>%
  dplyr::mutate(Variant_density_genomewide=total_variant_count/total_bin_count) %>%
  dplyr::select(STRAIN, Divergent_size_kb=bin_count, Fraction_divergent=div_fraction, 
                Variant_count_all=total_variant_count, Varaint_count_divergent=variant_sum, Fraction_variant_divergent=variant_fraction,
                Variant_density_genomewide, Variant_density_divergent=variant_mean) %>%
  dplyr::mutate(density_fold_divergent = Variant_density_divergent/Variant_density_genomewide)

DataS5 %>%
  ggplot(.) +
  geom_smooth() +
  geom_point(size = 1, color ='red', alpha = 0.5) +
  aes(x=Divergent_size_kb/1e3, y=density_fold_divergent) +
  theme_bw() +
  labs(x=)

write.csv(DataS5, file="Manuscript/Supplementary_Data/DataS5_divergent_regions_variant_counts.csv")

```

```
