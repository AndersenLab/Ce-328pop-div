library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

### Compare all vs arm (control) at layer 2 ###

all_genes <- nrow(read.csv(file="Processed_Data/WormCat/whole_genome_nov-16-2019.csv"))

for (i in c("div_all","nondiv_arms")) {
  group_genes <- nrow(read.csv(file=glue::glue("Processed_Data/WormCat/{i}/rgs_and_categories.csv")))
  for (j in 1:3) {
    assign(glue::glue("df_{i}_cat{j}"), read.csv(file=glue::glue("Processed_Data/WormCat/{i}/rgs_fisher_cat{j}_apv.csv")) %>%
             dplyr::mutate(layer=j, group = i, Total_genome = all_genes, group_total_genes = group_genes, GeneRatio = RGS/group_genes, BgRatio = AC/all_genes))
  }
}

df_all_div_cat <- rbind(df_div_all_cat1,df_div_all_cat2,df_div_all_cat3)
df_arm_cat <- rbind(df_nondiv_arms_cat1,df_nondiv_arms_cat2,df_nondiv_arms_cat3)

## For main figure (category 2)

df_cat2 <- rbind(df_all_div_cat, df_arm_cat) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::filter(Category != "Unknown") %>%
  dplyr::filter(layer == 2)

df_cat2$group <- factor(df_cat2$group, levels = c("div_all","nondiv_arms"), labels = c("Hyper-divergent regions","Non-divergent arms"))

df_cat2_plotpoint <- df_cat2 %>%
  dplyr::arrange(-Bonferroni) %>%
  dplyr::group_by(Category) %>%
  dplyr::filter(Bonferroni == min(Bonferroni)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Category, layer) %>%
  dplyr::group_by(layer) %>%
  dplyr::mutate(plotpoint=row_number()) %>%
  dplyr::ungroup()

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
  theme(axis.text = element_text(size=9, color='black'), 
        axis.title.x = element_text(size=10, color='black'), 
        axis.title.y = element_blank(), 
        legend.title = element_text(size=9.5, color='black'), 
        legend.text = element_text(size=9, color='black'), 
        legend.direction = "horizontal", 
        legend.box = "vertical", 
        legend.position = c(0.58,0.35), legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.2, "in"), 
        legend.spacing.y = unit(0.05, 'in'), 
        panel.grid = element_blank(),
        legend.margin=margin(0,0,0,0, unit= "in"),
        plot.margin=margin(0.1,0.1,0.1,0.1, unit= "in")) +
  guides(shape= guide_legend(nrow=2, order = 3, override.aes = list(size=3), title = NULL, title.position = "top"), size= guide_legend(nrow=1, order = 1, title.position = "top"), fill =  guide_colourbar(nrow=1, order =2, title.position = "top")) +
  labs(x="-log10(corrected p-value)", shape = "Frequency", size = "Fold enrichment", fill = "Gene counts") 

plot_wormcat_2

## for supplementary figure (category 1, 3)

# cat 1

df_cat1 <- rbind(df_all_div_cat, df_arm_cat) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::filter(Category != "Unknown") %>%
  dplyr::filter(layer == 1)

df_cat1$group <- factor(df_cat1$group, levels = c("div_all","nondiv_arms"), labels = c("Hyper-divergent regions","Non-divergent arms"))

df_cat1_plotpoint <- df_cat1 %>%
  dplyr::arrange(-Bonferroni) %>%
  dplyr::group_by(Category) %>%
  dplyr::filter(Bonferroni == min(Bonferroni)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Category, layer) %>%
  dplyr::group_by(layer) %>%
  dplyr::mutate(plotpoint=row_number()) %>%
  dplyr::ungroup()

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

# cat 3

df_cat3 <- rbind(df_all_div_cat, df_arm_cat) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::filter(Category != "Unknown") %>%
  dplyr::filter(layer == 3)

df_cat3$group <- factor(df_cat3$group, levels = c("div_all","nondiv_arms"), labels = c("Hyper-divergent regions","Non-divergent arms"))

df_cat3_plotpoint <- df_cat3 %>%
  dplyr::arrange(-Bonferroni) %>%
  dplyr::group_by(Category) %>%
  dplyr::filter(Bonferroni == min(Bonferroni)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Category, layer) %>%
  dplyr::group_by(layer) %>%
  dplyr::mutate(plotpoint=row_number()) %>%
  dplyr::ungroup()

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

save(df_cat1, df_cat1_plotpoint, df_cat2, df_cat2_plotpoint, df_cat3, df_cat3_plotpoint, file="Processed_Data/WormCat.RData")
