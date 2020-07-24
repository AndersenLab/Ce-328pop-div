library(tidyverse)
library(clusterProfiler) ## BiocManager::install("clusterProfiler")
library(org.Ce.eg.db) ##BiocManager::install("org.Ce.eg.db")
library(biomaRt) ##BiocManager::install("biomaRt")
library(enrichplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

### GO enrichment

## Gene enrichment test for the most swept regions

### GO enrichment tests

enrich_go <- function(x = NULL){
  
  # GO term enrichment  
  wb_ids <- dplyr::filter(x) %>% dplyr::pull(WBGeneID)
  
  mf <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "bonferroni",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  bp <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "bonferroni",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  return(list(mf, bp))
}

### Go enrichment test for the divergent regions
# set up data base to extract gene ids

gene_files_arm <- data.table::fread("Processed_Data/divergent_regions/genes/Non_div_arm_genes_only.tsv", header=F) %>%
  dplyr::select(WBGeneID = V1)

df_GO_nondiv_arm <- enrich_go(gene_files_arm)

gene_name_Pro_BP_arm <- setReadable(df_GO_nondiv_arm[[2]], OrgDb = org.Ce.eg.db)
barplot(gene_name_Pro_BP_arm, showCategory=100) 

gene_name_Pro_MF_arm <- setReadable(df_GO_nondiv_arm[[1]], OrgDb = org.Ce.eg.db)
barplot(gene_name_Pro_MF_arm, showCategory=100)

## all 
gene_files_div_all <- data.table::fread("Processed_Data/divergent_regions/genes/Divergent_genes_only.tsv", header=F) %>%
  dplyr::select(WBGeneID = V1)

df_GO_div_all <- enrich_go(gene_files_div_all)

gene_name_Pro_BP_div_all <- setReadable(df_GO_div_all[[2]], OrgDb = org.Ce.eg.db)
barplot(gene_name_Pro_BP_div_all, showCategory=100) 

gene_name_Pro_MF_div_all <- setReadable(df_GO_div_all[[1]], OrgDb = org.Ce.eg.db)
barplot(gene_name_Pro_MF_div_all)


##### GO enrichment summary df

### Biological process (BP)

df_GO_enrich_BP <- rbind(dplyr::mutate(gene_name_Pro_BP_arm[], freq="control_arms"), dplyr::mutate(gene_name_Pro_BP_div_all[], freq="all")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(p.adjust)
  
GO_list_BP <- df_GO_enrich_BP$Description

GO_list_BP_plotpoint <- data.frame(Description=GO_list_BP, plotpoint=length(GO_list_BP):1)

df_GO_enrich_BP_sum <- df_GO_enrich_BP %>%
  dplyr::filter(Description %in% GO_list_BP) %>%
  dplyr::left_join(., GO_list_BP_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_class_total_BP <- df_GO_enrich_BP_sum %>%
  dplyr::distinct(Description, class_gene_total, plotpoint) %>%
  dplyr::arrange(-plotpoint)

df_GO_enrich_BP_sum$freq <- factor(df_GO_enrich_BP_sum$freq, levels = c("all","control_arms"), labels = c("Divergent regions","Non-divergent arms"))

plot_GO_BP <- df_GO_enrich_BP_sum %>%
  ggplot(.) +
  geom_vline(xintercept = -log10(0.05), color='blue', size=0.4) +
  geom_point(alpha=0.9) +
  aes(x=-log10(p.adjust), y=plotpoint, shape=freq, size=EnrichRatio, fill=Count, alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "",sec.axis = dup_axis(labels = df_class_total_BP$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22)) +
  scale_size_continuous(range = c(2,4), breaks = c(2,3,4)) +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=11, color='black'), 
        legend.title = element_text(size=10, color='black'), legend.text = element_text(size=9, color='black'), 
        legend.position = c(0.77,0.25), panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="-log10(Corrected P-value)", shape = "Frequency", size = "Fold enrichment", fill = "Gene counts") 

plot_GO_BP

### Molecular function (MF)

df_GO_enrich_MF <- rbind(dplyr::mutate(gene_name_Pro_MF_arm[], freq="control_arms"), dplyr::mutate(gene_name_Pro_MF_div_all[], freq="all")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(p.adjust)

df_GO_enrich_MF$Description <- gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen", "oxidoreductase activity (a)", df_GO_enrich_MF$Description)
df_GO_enrich_MF$Description <- gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen", "oxidoreductase activity (b)", df_GO_enrich_MF$Description)

GO_list_MF <- df_GO_enrich_MF$Description

GO_list_MF_plotpoint <- data.frame(Description=GO_list_MF, plotpoint=length(GO_list_MF):1)

df_GO_enrich_MF_sum <- df_GO_enrich_MF %>%
  dplyr::filter(Description %in% GO_list_MF) %>%
  dplyr::left_join(., GO_list_MF_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_class_total_MF <- df_GO_enrich_MF_sum %>%
  dplyr::distinct(Description, class_gene_total, plotpoint) %>%
  dplyr::arrange(-plotpoint)

df_GO_enrich_MF_sum$freq <- factor(df_GO_enrich_MF_sum$freq, levels = c("all","control_arms"), labels = c("Divergent regions","Non-divergent arms"))

plot_GO_MF <- df_GO_enrich_MF_sum %>%
  ggplot(.) +
  geom_vline(xintercept = -log10(0.05), color='blue', size=0.4) +
  geom_point(alpha=0.9) +
  aes(x=-log10(p.adjust), y=plotpoint, shape=freq, size=EnrichRatio, fill=Count, alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_MF):1, labels = GO_list_MF, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_MF):1, labels = GO_list_MF, name = "",sec.axis = dup_axis(labels = df_class_total_MF$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22)) +
  scale_size_continuous(range = c(2,4), breaks = c(2,3,4)) +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=11, color='black'), 
        legend.title = element_text(size=10, color='black'), legend.text = element_text(size=9, color='black'), 
        legend.position = c(0.77,0.32), panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=2, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="-log10(Corrected P-value)", shape = "Frequency", size = "Fold enrichment", fill = "Gene counts") 

plot_GO_MF

plot_GO_BP_MF <- cowplot::plot_grid(plot_GO_BP, plot_GO_MF, ncol=1, rel_heights = c(1,1), labels = c("a", "b"), align = "hv", label_size = 12)

plot_GO_BP_MF

ggsave("Plots/supplementary_figures/plot_GO_BP_MF.pdf", plot = plot_GO_BP_MF, width = 7.5, height =9, units = "in")

save(GO_list_BP, GO_list_BP_plotpoint, df_GO_enrich_BP_sum, GO_list_MF, GO_list_MF_plotpoint, df_GO_enrich_MF_sum, file="Processed_Data/GO_enrich.RData")
