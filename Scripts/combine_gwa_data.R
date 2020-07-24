library(tidyverse)


setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("data/20150908_GWAS6residuals.Rda")
g6 <- pheno %>%
  dplyr::ungroup() %>%
  dplyr::select(strain, condition, trait, phenotype)%>%
  dplyr::mutate(condition_trait = paste(condition,trait,sep="_"))%>%
  dplyr::select(trait = condition_trait, strain, phenotype)
load("data/20150908_GWAS7residuals.Rda")
g7 <- pheno%>%
  dplyr::ungroup() %>%
  dplyr::select(strain, condition, trait, phenotype)%>%
  dplyr::mutate(condition_trait = paste(condition,trait,sep="_"))%>%
  dplyr::select(trait = condition_trait, strain, phenotype)
load("data/20150908_GWAS8residuals.Rda")
g8 <- pheno%>%
  dplyr::ungroup() %>%
  dplyr::select(strain, condition, trait, phenotype)%>%
  dplyr::mutate(condition_trait = paste(condition,trait,sep="_"))%>%
  dplyr::select(trait = condition_trait, strain, phenotype)



gwas <- dplyr::bind_rows(list(g6,g7,g8))%>%
  # head() %>%
  dplyr::mutate(trait = gsub("\\.","=", trait)) %>%
  dplyr::group_by(trait, strain)%>%
  dplyr::summarise(phenotype = mean(phenotype,na.rm=T))%>%
  tidyr::spread(trait, phenotype)%>%
  dplyr::ungroup()


save(gwas, file = "data/20200708_Combined_Processed_GWAS_Phenotypes.Rda")

bact_conds <- c("DA837","DA1877","HB101","HT115","JUb22","JUb23","JUb62","JUb67","JUb68","JUb71","JUb85","JUb86","JUb87","K12","OP50","OP50-B12","OP50-TRP")

bact_gwas <- dplyr::bind_rows(list(g6,g7,g8))%>%
  dplyr::mutate(trait = gsub("\\.","=", trait)) %>%
  dplyr::filter(grepl(paste(bact_conds,collapse = "|"), trait)) %>%
  dplyr::group_by(trait, strain)%>%
  dplyr::summarise(phenotype = mean(phenotype,na.rm=T))%>%
  tidyr::spread(trait, phenotype)%>%
  dplyr::ungroup()

save(bact_gwas, file = "data/20200708_Combined_Processed_GWAS_Bacterial_Phenotypes.Rda")
write.table(bact_gwas, file = "data/cegwas_bact_gwa_input.tsv", col.names = T, row.names = F, quote = F, sep = "\t")


# nextflow main.nf --traitfile=traits/cegwas_bact_gwa_input.tsv --vcf=bin/WI.20180527.impute.vcf.gz --p3d=FALSE --sthresh=EIGEN --out bacterial_gwas --freqUpper=0.5 --minburden=1

bact_conds <- c("JUb22","JUb23","JUb62","JUb67","JUb68","JUb71","JUb85","JUb86","JUb87")

bact_gwas <- dplyr::bind_rows(list(g6,g7,g8))%>%
  dplyr::mutate(trait = gsub("\\.","=", trait)) %>%
  dplyr::filter(grepl(paste(bact_conds,collapse = "|"), trait)) %>%
  dplyr::group_by(trait, strain)%>%
  dplyr::summarise(phenotype = mean(phenotype,na.rm=T))%>%
  tidyr::spread(trait, phenotype)%>%
  dplyr::ungroup()

# save(bact_gwas, file = "data/20200708_Combined_Processed_GWAS_Bacterial_Phenotypes.Rda")
write.table(bact_gwas, file = "data/cegwas_JUb_gwa_input.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
