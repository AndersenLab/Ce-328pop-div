#!/usr/bin/env Rscript
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
# setwd("..")

require(PopGenome)
require(WhopGenome)
require(data.table)
require(tidyverse)

# 1 - chromosome
# 2 - start of window
# 3 - end of window
# 4 - VCF file
# 5 - GFF file
# 6 - Save Label
# & - Strain file
# test args
# args <- c("II","7598325","8210489","Ce330_annotated.vcf.gz","WS245_exons.gff", "Arsenic", "249_samples.txt")

# args <- c("IV", "2", "18000000", "cbriggsae_cutter.vcf.gz", "../WS245_exons.gff", "Cbr", "cbr_strains.txt")
# Running script
# Rscript --vanilla Interval_Popgen.R II 10050022 12062611 Ce330_annotated.vcf.gz WS245_exons.gff Etoposide 249_samples.txt
args <- commandArgs(TRUE)

system(glue::glue("echo Initializing PopGenome Parameters"))

chr1 <- c(1,15072434)
chr2 <- c(1,15279421)
chr3 <- c(1,13783801)
chr4 <- c(1,17493829)
chr5 <- c(1,20924180)
chr6 <- c(1,17718942)

chr.lengths <- list(chr1,chr2,chr3,chr4,chr5,chr6)
chroms <- c("I","II","III","IV","V","X")

ANALYSIS_CHROM <- as.character(args[1])
CHROM_START <- chr.lengths[which(chroms == as.character(args[1]))][[1]][1]
CHROM_END <- chr.lengths[which(chroms == as.character(args[1]))][[1]][2]

SLIDE_DISTANCE <- 1000
WINDOW_SIZE <-  10000

REGION_START <- as.numeric(args[2])
REGION_END <- as.numeric(args[3])

OUTGROUP <- "JU1341"

system(glue::glue("echo Done Initializing PopGenome Parameters - WindowSize = {WINDOW_SIZE}, StepSize = {SLIDE_DISTANCE}, Whole Population, Chromosome = {ANALYSIS_CHROM}"))

system(glue::glue("bcftools view -S {args[7]} -r {ANALYSIS_CHROM}:{args[2]}-{args[3]} -v snps {args[4]} | bcftools plugin fixploidy - -- --force-ploidy 2 | sed 's/0\\/0/0|0/g' | sed 's/1\\/1/1|1/g' | sed 's/0\\/1/1|1/g' | sed 's/1\\/0/1|1/g' | sed 's/.\\/./.|./g' | bcftools filter -i 'F_MISSING < 0.1' | bcftools view -Oz -o tempRegion.vcf.gz"))
system(glue::glue("tabix tempRegion.vcf.gz"))
system(glue::glue("bcftools query -l tempRegion.vcf.gz > samples.txt"))

samples <- data.table::fread("samples.txt",header = F) %>%
  dplyr::pull(V1)

POPGENOME_VCF <- "tempRegion.vcf.gz"

system(glue::glue("echo PopGenome - Reading VCF file {ANALYSIS_CHROM}")) 

vcf_handle <- WhopGenome::vcf_open(POPGENOME_VCF)

GENOME_OBJECT <- Whop_readVCF(
  vcf_handle, 
  numcols = 100, 
  tid = ANALYSIS_CHROM, 
  from = REGION_START,
  to = REGION_END,
  include.unknown = TRUE)

system(glue::glue("echo PopGenome - Setting Outgroup and Defining Window Size"))

GENOME_OBJECT <- PopGenome::set.populations(GENOME_OBJECT, list(samples), diploid = FALSE)

GENOME_OBJECT <- PopGenome::set.outgroup(GENOME_OBJECT, OUTGROUP,  diploid = FALSE)

GENOME_OBJECT <- sliding.window.transform(GENOME_OBJECT, WINDOW_SIZE, SLIDE_DISTANCE, type=2)

GENOME_OBJECT <- PopGenome::F_ST.stats(GENOME_OBJECT, mode = "nucleotide", detail = T,FAST = T)

system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Detail Stats"))
GENOME_OBJECT <- PopGenome::detail.stats(GENOME_OBJECT)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Neutrality Stats"))
GENOME_OBJECT <- PopGenome::neutrality.stats(GENOME_OBJECT, detail = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Diversity Stats"))
GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT, pi = TRUE)


system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

save(GENOME_OBJECT, file = glue::glue("{args[6]}_{args[1]}_Statistics.Rda"))

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Neutrality and Diversity Stats"))

n_df <- data.frame(get.neutrality(GENOME_OBJECT)[[1]]) %>%
  tibble::rownames_to_column() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CHROM = args[1],
                startWindow = as.numeric(strsplit(rowname,split = " - ")[[1]][1]),
                endWindow = as.numeric(strsplit(strsplit(rowname,split = " - ")[[1]][2], split = " ")[[1]][1])) %>%
  dplyr::select(CHROM:endWindow, Tajima.D:Zeng.E)

d_df <- data.frame(get.diversity(GENOME_OBJECT)[[1]]) %>%
  tibble::rownames_to_column() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CHROM = args[1],
                startWindow = as.numeric(strsplit(rowname,split = " - ")[[1]][1]),
                endWindow = as.numeric(strsplit(strsplit(rowname,split = " - ")[[1]][2], split = " ")[[1]][1])) %>%
  dplyr::select(CHROM:endWindow, nuc.diversity.within:Pi)

d2_df <- dplyr::left_join(n_df,d_df, by = c("CHROM", "startWindow", "endWindow")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(theta_Watterson = c(GENOME_OBJECT@theta_Watterson)) %>%
  tidyr::gather(statistic, value, -(CHROM:endWindow))

save(d2_df, file = glue::glue("{args[6]}_{args[1]}_{args[2]}-{args[3]}_Diversity_Statistics.Rda"))

system("rm tempRegion.vcf.gz")
system("rm tempRegion.vcf.gz.tbi")
system("rm samples.txt")
