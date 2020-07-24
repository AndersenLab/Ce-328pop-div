library(ggrepel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(uwot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/strain_geo.RData")

df_pcs_noremoval_complete <- data.table::fread("Processed_Data/eigen_results/eigenstrat_no_removal_complete.evac", skip = 1) %>%
  dplyr::select(strain=V1, V2:V16) %>% 
  dplyr::rename(PC1=V2, PC2=V3,PC3=V4, PC4=V5, PC5=V6) %>%
  dplyr::left_join(., dplyr::select(indep_strain_info_geo, strain, reference_strain, country, state, geo), by='strain') %>%
  dplyr::filter(reference_strain == 1)

df_pcs_out_removal_complete <- data.table::fread("Processed_Data/eigen_results/eigenstrat_outliers_removed_complete.evac", skip = 1) %>%
  dplyr::select(strain=V1, V2:V16) %>% 
  dplyr::rename(PC1=V2, PC2=V3,PC3=V4, PC4=V5, PC5=V6) %>%
  dplyr::left_join(., dplyr::select(indep_strain_info_geo, strain, reference_strain, country, state, geo), by='strain') %>%
  dplyr::filter(reference_strain == 1)

df_pcs_noremoval_complete_global <- data.table::fread("Processed_Data/eigen_results/eigenstrat_no_removal_complete_global.evac", skip = 1) %>%
  dplyr::select(strain=V1, V2:V16) %>% 
  dplyr::rename(PC1=V2, PC2=V3,PC3=V4, PC4=V5, PC5=V6) %>%
  dplyr::left_join(., dplyr::select(indep_strain_info_geo, strain, reference_strain, country, state, geo), by='strain') %>%
  dplyr::filter(reference_strain == 1)

non_outlier <- unique(df_pcs_out_removal_complete$strain)


### PCA with complete site vcf

Kauai_strains <- c("XZ1513", "XZ1514", "XZ1515", "XZ1516", "ECA372", "ECA701")

plot_PC12_all <- df_pcs_noremoval_complete %>%
  ggplot(.) +
  geom_text_repel(data=dplyr::filter(df_pcs_noremoval_complete, (PC1 < 0.18 & PC1>0.1) ),  
                  label=dplyr::filter(df_pcs_noremoval_complete,  (PC1 < 0.18 & PC1>0.1))$strain, 
                  aes(x=PC1,y=PC2), size = 2.5, force=10, nudge_y=0.05,
                  segment.size = 0.2, segment.alpha = 0.5, direction = 'both',
                  min.segment.length = 0) +  
  geom_point(shape=21, alpha=0.8, size=2, aes(x=PC1,y=PC2, fill=geo))+
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position='none',
        axis.title.x=element_blank(),
        panel.grid = element_blank())+
  labs(x="PC1", y="PC2", fill ="")  +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  guides(fill= guide_legend(nrow=2))

plot_PC12_all

## define 3 clusters

df_ce_clusters <- df_pcs_noremoval_complete %>%
  dplyr::select(strain, PC1, PC2) %>%
  dplyr::mutate(group = ifelse(PC1 < 0.1, "Global", ifelse(PC2>0.2, "Hawaiian", ifelse(PC1>0.2, "Pacific", "Divergent"))))

df_ce_clusters %>%
  ggplot(.) +
  geom_text_repel(data=dplyr::filter(df_pcs_noremoval_complete, (PC1 < 0.18 & PC1>0.1) ),  
                  label=dplyr::filter(df_pcs_noremoval_complete,  (PC1 < 0.18 & PC1>0.1))$strain, 
                  aes(x=PC1,y=PC2), size = 2.5, force=10, nudge_y=0.05,
                  segment.size = 0.2, segment.alpha = 0.5, direction = 'both',
                  min.segment.length = 0) +  
  geom_point(shape=21, alpha=0.8, size=2, aes(x=PC1,y=PC2, fill=group))+
  #scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position='none',
        axis.title.x=element_blank(),
        panel.grid = element_blank())+
  labs(x="PC1", y="PC2", fill ="")  +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  guides(fill= guide_legend(nrow=2))

save(df_pcs_noremoval_complete,df_pcs_out_removal_complete,df_pcs_noremoval_complete_global, df_ce_clusters, file = "Processed_Data/PCA.RData")

### zoom in to three cluster

plot_PC12_all_zoom_ww <- df_pcs_noremoval_complete %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8, size=2, aes(x=PC1,y=PC2, fill=geo))+
  geom_text_repel(data=dplyr::filter(df_pcs_noremoval_complete, strain %in% bi_strain),  
                  label=dplyr::filter(df_pcs_noremoval_complete, strain %in% bi_strain)$strain, 
                  aes(x=PC1,y=PC2), size = 2.5, force=10, nudge_y=0.02,
                  segment.size = 0.2, segment.alpha = 0.5, direction = 'both',
                  min.segment.length = 0) +  
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.position='none',
        panel.grid = element_blank())+
  labs(fill ="")  +
  scale_x_continuous(limits = c(-0.02,0.03)) +
  ylim(-0.02,0.04)+
  guides(fill= guide_legend(nrow=2))

plot_PC12_all_zoom_ww

plot_PC12_all_zoom_hw <- df_pcs_noremoval_complete %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8, size=2, aes(x=PC1,y=PC2, fill=geo))+
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        axis.title.x=element_blank(),
        legend.position='none',
        panel.grid = element_blank())+
  labs(fill ="")  +
  xlim(0.20,0.22) +
  ylim(0.26,0.29)

plot_PC12_all_zoom_hw

plot_PC12_all_zoom_pa <- df_pcs_noremoval_complete %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8, size=2, aes(x=PC1,y=PC2, fill=geo))+
  geom_text_repel(data=dplyr::filter(df_pcs_noremoval_complete, strain %in% bi_strain),  
                  label=dplyr::filter(df_pcs_noremoval_complete, strain %in% bi_strain)$strain, 
                  aes(x=PC1,y=PC2), size = 2.5, force=10, nudge_y=-0.01,
                  segment.size = 0.2, segment.alpha = 0.5, direction = 'both',
                  min.segment.length = 0) +  
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        axis.title.y=element_blank(),
        legend.position='none',
        panel.grid = element_blank())+
  labs(fill ="")  +
  xlim(0.20,0.245) +
  scale_y_continuous(limits = c(-0.23,-0.15), 
                     breaks = c(-0.16, -0.18, -0.2, -0.22))
  
plot_PC12_all_zoom_pa

PCA_allpop_merged <- cowplot::plot_grid(plot_PC12_all, plot_PC12_all_zoom_hw, plot_PC12_all_zoom_ww, 
                                        plot_PC12_all_zoom_pa, 
                                 labels = c("c","d","e","f"), 
                                 label_size = 12, nrow=2,
                                 align = 'hv', axis='tblr')  

ggsave(PCA_allpop_merged, file="PCA_allpop_merge.pdf", width = 7.5, height=5.3, unit='in')


### global only PCA

plot_PC12_global <- df_pcs_noremoval_complete_global %>%
  dplyr::filter(geo != "Unknown") %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8, size=1.5, aes(x=PC1,y=PC2, fill=geo))+
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.title = element_blank(),
        legend.position ='none',
        panel.grid = element_blank())+
  labs(x="PC1", y="PC2", fill ="")

plot_PC12_global


##### outlier removed

plot_PC12_pruned <- df_pcs_out_removal_complete %>%
  dplyr::filter(geo != "Unknown") %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8, size=1.5, aes(x=PC1,y=PC2, fill=geo))+
  scale_fill_manual(values = geo.colours) +
  theme_bw() +
  theme(axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.title = element_blank(),
        legend.position ='none',
        panel.grid = element_blank())+
  labs(x="PC1", y="PC2", fill ="")

plot_PC12_pruned

### UMAP plot

file_name = "eigenstrat_outliers_removed_complete"
take_pcs=5
n_neigh=60

pcs <- data.table::fread(glue::glue("Processed_Data/eigen_results/{file_name}.evac"), skip = 1) %>% data.frame()
cnames <- colnames(pcs)[c(1:take_pcs+1)]
pcs <- data.frame(strain = pcs$V1, pcs[,cnames])

temp_umap <- data.frame(umap(pcs, n_neighbors = n_neigh, learning_rate = 0.5, min_dist = 0.3),
                              strain = pcs$strain)%>%
        dplyr::left_join(.,dplyr::select(indep_strain_info_geo, strain, country, state, geo), by='strain') %>%
        dplyr::mutate(n_pcs = take_pcs,
                      n_nei = n_neigh,
                      vcf=file_name)

#save(temp_umap, file = "Processed_Data/umap.RData")
load("Processed_Data/umap.RData")

plot_umap <- temp_umap %>%
  dplyr::filter(geo != "Unknown") %>%
        ggplot(.)+
        aes(x=X1,y=X2)+
        geom_point(size = 1.5, shape = 21, aes(fill = geo), alpha = 0.8)+
        scale_fill_manual(values=geo.colours) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.title=element_blank(),
              axis.title = element_text(size=11, color = "black"), 
              axis.text = element_text(size=10, color = "black"),
              legend.position ='none') +
        labs(x = "UMAP1", y = "UMAP2")

plot_umap

### merged plots


PCA_merged <- cowplot::plot_grid(plot_PC12_all, plot_PC12_pruned, plot_umap, 
                                 labels = c("c","d","e"), rel_heights = c(1,1,1), 
                                 label_size = 12, nrow=1,
                                 align = 'h', axis='tblr')  

ggsave(PCA_merged, file="PCA.pdf", width = 7.5, height=3, unit='in')








###  multi PC plot

p1 <- ggpairs(df_pcs_noremoval_complete, columns = c(2:5), upper= 'blank', diag = 'blank', 
              lower = list(continuous = wrap("points", alpha = 0.8, size = 1.3, shape =21)),
              mapping = ggplot2::aes(fill = geo)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        #axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())

for(i in 1:p1$nrow) {
  for(j in 1:p1$ncol){
    p1[i,j] <- p1[i,j] + 
      scale_fill_manual(values = geo.colours)
  }
}

p1

p2 <- ggpairs(df_pcs_out_removal_complete, columns = c(2:5), lower= 'blank', diag = 'blank', 
              upper = list(continuous = wrap("points", alpha = 0.8, size = 1.3, shape =21)),
              mapping = ggplot2::aes(fill = geo)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())

for(i in 1:p2$nrow) {
  for(j in 1:p2$ncol){
    p2[i,j] <- p2[i,j] + 
      scale_fill_manual(values = geo.colours)
  }
}

p2

p12 <- p2 

for(i in 1:p12$nrow) {
  for(j in 1:p12$ncol){
    if (i>j) {
      
      p12[i,j] <- p1[i,j]
    }
    else {
      p12[i,j] <- p2[i,j] 
    }
    
  }
}

p12

ggsave("Plots/Supplementary/PCA_complete.pdf", p12, width = 7.5, height = 7.5, units = "in")



### UMAP for EU strains

file_name_eu = "eigenstrat_outliers_removed_complete"
take_pcs_eu=10
n_neigh_eu=20

pcs_eu <- data.table::fread(glue::glue("Processed_Data/eigen_results/{file_name}.evac"), skip = 1) %>% 
  data.frame() %>%
  dplyr::filter(V1 %in% eu_strain)
cnames_eu <- colnames(pcs_eu)[c(1:take_pcs_eu+1)]
pcs_eu <- data.frame(strain = pcs_eu$V1, pcs_eu[,cnames])

temp_umap_eu <- data.frame(umap(pcs_eu, n_neighbors = n_neigh_eu, learning_rate = 0.5, min_dist = 0.1),
                           strain = pcs_eu$strain)%>%
  dplyr::rename(isotype=strain) %>%
  dplyr::left_join(.,dplyr::select(indep_strain_info_geo, isotype, country, state, geo_isotype), by='isotype') %>%
  dplyr::mutate(n_pcs = take_pcs_eu,
                n_nei = n_neigh_eu,
                vcf=file_name)

temp_umap_eu %>%
  dplyr::filter(geo_isotype != "Unknown") %>%
  ggplot(.)+
  aes(x=X1,y=X2)+
  geom_point(size = 2.5, shape = 21, aes(fill = country), alpha = 0.8)+
  scale_fill_manual(values=eu.colours) +
  theme_bw(18) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=9, color = "black"),
        legend.position='top',
        legend.direction = "horizontal", 
        legend.key.size = unit(0.05, "in"), 
        legend.key.width = unit(0.05,"in")) +
  labs(x = "UMAP1", y = "UMAP2")

#save(temp_umap_eu, file = "Processed_Data/umap_eu.RData")
load("Processed_Data/umap_eu.RData")

plot_umap_eu <- temp_umap_eu %>%
  dplyr::filter(geo_isotype != "Unknown") %>%
  ggplot(.)+
  aes(x=X1,y=X2)+
  geom_point(size = 2.5, shape = 21, aes(fill = country), alpha = 0.8)+
  scale_fill_manual(values=eu.colours) +
  theme_bw(18) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11, color = "black"), 
        axis.text = element_text(size=10, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=9, color = "black"),
        legend.position='top',
        legend.direction = "horizontal", 
        legend.key.size = unit(0.05, "in"), 
        legend.key.width = unit(0.05,"in")) +
  labs(x = "UMAP1", y = "UMAP2")

plot_umap_eu



