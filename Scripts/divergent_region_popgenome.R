library(dplyr)
library(tidyr)
library(ggplot2)
library(scatterpie)

### pi, theta, tajima's D for divergent regions

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))


chr <- c("I","II","III","IV","V","X")

for(sm in chr) {
temp_neutrality_df <- read.csv(glue::glue("Processed_Data/popgenome/whole_pop/chrom{sm}_td.csv")) %>%
  dplyr::select(-X) %>%
  dplyr::rename(Pi=pi, theta_Watterson=theta, Tajima.D=td, CHROM=chrom, startWindow=start, endWindow=end) %>%
  dplyr::mutate(Pi=Pi, theta_Watterson=theta_Watterson) %>%
  tidyr::gather(-CHROM,-startWindow,-endWindow, key="statistic", value="value")

if(!exists("neutrality_df")){
  neutrality_df <- temp_neutrality_df
} else {
  neutrality_df <- dplyr::bind_rows(neutrality_df, temp_neutrality_df)
}
}

load("Processed_Data/divergent_classification.RData")
load("Processed_Data/Ce_Genome-wide_Neutrality_stats.Rda")

df_divergent_freq_simple <- df_divergent_final %>%
  dplyr::mutate(startWindow=START_BIN+1, endWindow=END_BIN) %>%
  dplyr::distinct(CHROM, startWindow, endWindow,freq_bin)

neutrality_df_div <- neutrality_df %>%
  dplyr::left_join(., df_divergent_freq_simple, by=c('CHROM','startWindow','endWindow')) %>%
  dplyr::mutate(freq_bin=ifelse(is.na(freq_bin), "Non-divergent", as.character(freq_bin)))

neutrality_df_div$freq_bin <- factor(neutrality_df_div$freq_bin, levels = c("Non-divergent","Rare","Intermediate", "Common"))
neutrality_df_div$statistic <- factor(neutrality_df_div$statistic, levels = c("theta_Watterson","Pi","Tajima.D"), labels = c("Watterson's estimator (θ)","Nucleotide diversity (π)","Tajima's D"))

plot_div_popgenome <- neutrality_df_div %>%
  na.omit() %>%
  ggplot(.) +
  geom_violin() +
  #geom_boxplot(outlier.alpha = 0) +
  #geom_jitter(size=0.1, alpha = 0.2) +
  aes(x=freq_bin, y=value) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.5, hjust=1), axis.text.y = element_text(size=10, color='black')) +
  facet_wrap(~statistic, scales = 'free', nrow=1)

plot_div_popgenome

ggsave("Plots/plot_div_popgenome.pdf", plot_div_popgenome, width =  6, height = 6)

non_divergent_value <- dplyr::filter(neutrality_df_div, statistic == "Nucleotide diversity (π)", freq_bin == "Non-divergent")$value
divergent_value <- dplyr::filter(neutrality_df_div, statistic == "Nucleotide diversity (π)", freq_bin != "Non-divergent")$value
common_divergent_value <- dplyr::filter(neutrality_df_div, statistic == "Nucleotide diversity (π)", freq_bin == "Common")$value

mean(divergent_value)/mean(non_divergent_value)
mean(common_divergent_value)/mean(non_divergent_value)



tajimad_div <- neutrality_df_div %>%
  dplyr::filter(statistic=="Tajima's D") %>%
  na.omit()

plot_div_TajD1 <- tajimad_div %>%
  na.omit() %>%
  ggplot(.) +
  #geom_violin() +
  geom_boxplot(outlier.alpha = 0) +
  #geom_jitter(size=0.1, alpha = 0.2) +
  aes(x=freq_bin, y=value) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=11, color='black'),
        axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.5, hjust=1), 
        axis.text.y = element_text(size=10, color='black')) +
  labs(y="Tajima's D") +
  coord_cartesian(ylim=c(-3,3.5))

plot_div_TajD1

aov_Taj <- aov(tajimad_div$value ~ tajimad_div$freq_bin) 
summary(aov_Taj)
TukeyHSD(aov_Taj)

View(pairwise.t.test(tajimad_div$value, tajimad_div$freq_bin, p.adj = "bonf")$p.value)

### distribution of pi theta

neutrality_df_div2 <- neutrality_df %>%
  dplyr::left_join(., df_divergent_freq_simple, by=c('CHROM','startWindow','endWindow')) %>%
  dplyr::mutate(freq_bin=ifelse(is.na(freq_bin), "Non-divergent", as.character(freq_bin)))

neut_spreaad_div <- neutrality_df_div2 %>%
  dplyr::filter(statistic %in% c("Pi", "theta_Watterson")) %>%
  tidyr::spread(key="statistic", value="value")

neut_spreaad_div %>%
  ggplot(.) +
  geom_point(size=0.5, alpha = 0.5) +
  aes(x=Pi, y=theta_Watterson) +
  theme_bw() +
  theme(panel.grid=element_blank())






plot_div_TajD2 <- tajimad_div %>%
  dplyr::mutate(class=ifelse(freq_bin == "Non-divergent", "Non-divergent","Divergent")) %>%
  na.omit() %>%
  ggplot(.) +
  geom_boxplot() +
  #geom_boxplot(outlier.alpha = 0) +
  #geom_jitter(size=0.1, alpha = 0.2) +
  aes(x=class, y=value) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.5, hjust=1), axis.text.y = element_text(size=10, color='black')) +
  facet_wrap(~statistic, scales = 'free', nrow=1)

plot_div_TajD2

### Enrichment of divergent regions in balanced loci

quantile(tajimad_div$value, probs = 0.99, na.rm = TRUE)

tajimad_div_balanced <- tajimad_div %>%
  dplyr::filter(value >= 2)
  #dplyr::filter(value >= quantile(tajimad_div$value, probs = 0.99, na.rm = TRUE))

tajimad_div_balanced %>%
  ggplot(.) +
  geom_bar() +
  aes(x=freq_bin) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.5, hjust=1), axis.text.y = element_text(size=10, color='black')) +
  labs(y="Size (kb)")

nrow(dplyr::filter(tajimad_div_balanced, freq_bin != "Non-divergent"))/nrow(tajimad_div_balanced)


df_divergent_window <- df_div_cluster_gapjoined_size %>%
  dplyr::filter(!is.na(div_class)) %>%
  dplyr::mutate(startWindow=START_BIN+1, endWindow=END_BIN) %>%
  dplyr::distinct(CHROM, startWindow, endWindow,div_class)

tajimad_div_balanced_window <- tajimad_div_balanced %>%
  dplyr::left_join(., df_divergent_window, by=c('CHROM', 'startWindow', 'endWindow'))

tajimad_div_balanced_window %>%
  ggplot(.) +
  geom_bar() +
  aes(x=div_class) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=10, color='black', vjust=0.5, hjust=1), axis.text.y = element_text(size=10, color='black')) +
  labs(x="bin class", y="Size (kb)")

nrow(dplyr::filter(tajimad_div_balanced_window, !is.na(div_class)))/nrow(tajimad_div_balanced_window)
