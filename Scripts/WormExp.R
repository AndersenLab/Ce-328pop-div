library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

### WormExp - microbes

microbes_div <- read.csv(file="Processed_Data/divergent_regions/genes/Wormexp_divergent_microbes.csv") %>%
  dplyr::mutate(is_divergent = "Divergent")
microbes_nondiv <- read.csv(file="Processed_Data/divergent_regions/genes/Wormexp_nondivergent_microbes.csv") %>%
  dplyr::mutate(is_divergent = "Non-divergent")

WormExp_microbe <- rbind(microbes_div,microbes_nondiv)

write.csv(WormExp_microbe, file="Processed_Data/divergent_regions/genes/Wormexp_microbes.csv") 

WormExp_microbe_strain <- read.csv(file="Processed_Data/divergent_regions/genes/Wormexp_microbes_strain.csv") %>%
  tidyr::separate(Term, into=c('direction', 'etc'))

WormExp_microbe_strain %>%
  dplyr::filter(direction=="UP") %>%
  dplyr::filter(Bonferroni < 0.05) %>%
  ggplot(.) +
  geom_bar() +
  aes(x=Strain, fill=is_divergent, position="identity") +
  theme_bw() +
  xlab("Microbes")

### WormExp - chemical/stress

chemical_stress_div <- read.csv(file="Processed_Data/divergent_regions/genes/Wormexp_divergent_chemical_stress.csv") %>%
  dplyr::mutate(is_divergent = "Divergent")
chemical_stress_nondiv <- read.csv(file="Processed_Data/divergent_regions/genes/Wormexp_nondivergent_chemical_stress.csv") %>%
  dplyr::mutate(is_divergent = "Non-divergent")

WormExp_chemical_stress <- rbind(chemical_stress_div,chemical_stress_nondiv) %>%
  dplyr::mutate(Term2=Term) %>%
  tidyr::separate(Term2, into=c('direction', 'etc')) %>%
  dplyr::mutate(class=ifelse(direction %in% c("Up", "UP"), "Up-regulated",
                             ifelse(direction %in% c("Down", "down"), "Down-regulated", NA)))

WormExp_chemical_stress %>%
  dplyr::filter(Bonferroni < 0.05, !is.na(class)) %>%
  ggplot(.) +
  geom_bar() +
  aes(x=is_divergent, position="identity") +
  theme_bw() +
  xlab("chemical_stress") +
  facet_grid(~class)

