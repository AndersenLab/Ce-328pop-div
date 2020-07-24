library(dplyr)
library(ggplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

LDdecay_div <- read.table(file="Processed_Data/LD/Ce328_complete_sites-MAF05-autoarm-divergent_popLDdecay.stat") %>%
  dplyr::rename(dist=V1, mean_r2=V2, sum_r2=V4, n_pairs=V6) %>%
  dplyr::select(-V3,-V5) %>%
  dplyr::mutate(class="Hyper-divergent")

LDdecay_nondiv<- read.table(file="Processed_Data/LD/Ce328_complete_sites-MAF05-autoarm-non_divergent_popLDdecay.stat") %>%
  dplyr::rename(dist=V1, mean_r2=V2, sum_r2=V4, n_pairs=V6) %>%
  dplyr::select(-V3,-V5) %>%
  dplyr::mutate(class="Non-divergent")

df_LDdecay <- rbind(LDdecay_div,LDdecay_nondiv) %>%
  dplyr::filter(n_pairs >= 10)

#save(df_LDdecay, file="Processed_Data/LD_decay.RData")
load("Processed_Data/LD_decay.RData")

plot_LDdecay <- df_LDdecay %>%
  dplyr::filter(dist >= 10) %>%
  #dplyr::filter(n_pairs >= 20) %>%
  ggplot(.) +
  aes(x=dist/1e3, y=mean_r2) +
  geom_point(size = 0.01, alpha =0.2, aes(color = class)) +
  geom_smooth(se = T, size = 0.7, n = 5000, aes(group=class, linetype = class), color='black') + ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
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
        legend.key.size = unit(0.08,"in")) +
  labs(x = "Distance (kb)", y = "Linkage disequilibrium (r2)", color = "") +
  guides(colour = guide_legend(override.aes = list(size=1.5, alpha = 1)))

plot_LDdecay

ggsave(plot_LDdecay, file="Plots/plot_LDdecay_complete_sites_autosome_arms.png", width = 5, height=4)

plot_LDdecay_main <- df_LDdecay %>%
  dplyr::filter(dist >= 10) %>%
  ggplot(.) +
  aes(x=dist/1e3, y=mean_r2) +
  #geom_point(size = 0.1, alpha =0.4, aes(color = class)) +
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
        legend.title=element_blank(),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.08,"in")) +
  labs(x = "Distance (kb)", y = "Linkage disequilibrium (r2)", color = "") +
  guides(linetype = guide_legend(keywidth = 3, keyheight = 1))

plot_LDdecay_main

ggsave(plot_LDdecay_main, file="Plots/plot_LDdecay_main.png", width = 5, height=4)

plot_LDdecay_supple <- df_LDdecay %>%
  dplyr::filter(dist >= 10) %>%
  ggplot(.) +
  aes(x=dist/1e3, y=mean_r2) +
  geom_point(size = 0.1, alpha =0.4, aes(color = class)) +
  geom_smooth(se = F, size = 0.7, n = 5000, aes(group=class, linetype = class), color='black') + ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
  scale_color_brewer(palette = 'Set1') +
  scale_linetype(guide=FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'), 
        axis.title.y = element_text(size=11, color='black'), 
        axis.text.y= element_text(size=10, color='black'), 
        axis.text.x = element_text(size=10, color='black'), 
        legend.position = 'none', panel.grid = element_blank(),
        legend.text =element_text(size=10, color='black'),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.08,"in")) +
  labs(x = "Distance (kb)", y = "Linkage disequilibrium (r2)", color = "") +
  guides(colour = guide_legend(override.aes = list(size=1.5, alpha = 1))) +
  facet_grid(~class)

plot_LDdecay_supple

ggsave(plot_LDdecay_supple, file="Plots/plot_LDdecay_supple.png", width = 7.5, height=4)



plot_LDdecay_npairs <- df_LDdecay %>%
  dplyr::filter(dist >= 10) %>%
  ggplot(.) +
  aes(x=dist/1e3, y=n_pairs) +
  geom_point(size = 0.1, alpha =0.6, aes(color = class)) +
  geom_smooth(se = F, size = 0.7, n = 5000, aes(group=class, linetype = class), color='black') + ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
  scale_color_brewer(palette = 'Set1') +
  scale_linetype(guide=FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color='black'), 
        axis.text.y= element_text(size=10, color='black'), 
        axis.text.x = element_text(size=10, color='black'), 
        legend.position = c(0.8,0.8), panel.grid = element_blank(),
        legend.text =element_text(size=10, color='black')) +
  labs(x = "Distance (kb)", y =  "Number of pairs", color = "") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))

plot_LDdecay_npairs

ggsave(plot_LDdecay_npairs, file="Plots/plot_LDdecay_npairs.pdf", width = 7.5, height=5)
