library(dplyr)
library(tidyr)
library(ggplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load(file="Raw/WI_328_list.RData")
load("Processed_Data/indep_strain_info.RData")

### Geographic grouping

### Geographic distribution of 351 independent collection / at least one degree apart from each other (latitude: at least 110.57 km away, longitude: 111.32 km away at the equator, 85.39 km away at latitude 40 degree)

world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica

ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat, fill = as.factor(release)), shape =21, size = 1.2, alpha = 0.9) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  labs(fill = "Release")

########### Plotting for the identification of geographic structure

#1. North America

plot_north_america <- ggplot()+ geom_map(data=world, map=world,
                                         aes(x=long, y=lat, map_id=region),
                                         color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill='red') +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-140, -70), y=c(10,90))

plot_north_america

df_na_strain <- indep_strain_info %>%
  dplyr::filter(long > -140  & long < -70) %>%
  dplyr::filter(lat > 10 & lat < 90)

na_strain <- as.character(df_na_strain$strain)
na_isotype <- unique(as.character(df_na_strain$isotype))

#2. South America

plot_south_america <- ggplot()+ geom_map(data=world, map=world,
                                         aes(x=long, y=lat, map_id=region),
                                         color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill='red') +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-100, -30), y=c(-70,10))

plot_south_america

df_sa_strain <- indep_strain_info %>%
  dplyr::filter(long > -100  & long < -30) %>%
  dplyr::filter(lat > -70 & lat < 10)

sa_strain <- as.character(df_sa_strain$strain)
sa_isotype <- unique(as.character(df_sa_strain$isotype))

#3. Hawaii

plot_hw <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill = 'red') +
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

#4. New zealand

plot_NZ <- ggplot()+ geom_map(data=world, map=world,
                                         aes(x=long, y=lat, map_id=region),
                                         color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill='red') +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(165, 180), y=c(-50,-25))

plot_NZ

df_nz_strain <- indep_strain_info %>%
  dplyr::filter(long > 165  & long < 180) %>%
  dplyr::filter(lat > -50 & lat < -25)

nz_strain <- as.character(df_nz_strain$strain)
nz_isotype <- unique(as.character(df_nz_strain$isotype))

#5. Australia

plot_AU <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill = 'red') +
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

#6. Atlantic - Azores/Madeira, Sao Tome

plot_at <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill = 'red') +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-40, -10), y=c(20,50)) 

plot_at

df_at_strain <- indep_strain_info %>%
  dplyr::filter(state %in% c("Azores", "Madeira", "Sao Tome"))

at_strain <- as.character(df_at_strain$strain)
at_isotype <- unique(as.character(df_at_strain$isotype))

#7. Africa

## African strains

plot_af <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill = 'red') +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-20, 60), y=c(-40,30)) 

plot_af

df_af_strain <- indep_strain_info %>%
  dplyr::filter(!strain %in% at_strain) %>%
  dplyr::filter(long > -20  & long < 60) %>%
  dplyr::filter(lat > -40 & lat < 30) 

af_strain <- as.character(df_af_strain$strain)
af_isotype <- unique(as.character(df_af_strain$isotype))

#8. Europe

plot_EU <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill = 'red') +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-10, 25), y=c(35, 57)) 

plot_EU

df_eu_strain <- indep_strain_info %>%
  dplyr::filter(long > -10  & long < 25) %>%
  dplyr::filter(lat > 35 & lat < 57) 

eu_strain <- as.character(df_eu_strain$strain)
eu_isotype <- unique(as.character(df_eu_strain$isotype))

#9. Asia

df_as_strain <- indep_strain_info %>%
  dplyr::filter(country %in% c("China", "Japan", "Taiwan", "Israel"))

as_strain <- as.character(df_as_strain$strain)
as_isotype <- unique(as.character(df_as_strain$isotype))

#10. others

df_etc_strain <- indep_strain_info %>%
  dplyr::filter(!strain %in% c(af_strain, at_strain, au_strain, eu_strain, 
                               hw_strain, na_strain, nz_strain,
                               sa_strain, as_strain))

etc_strain <- as.character(df_etc_strain$strain)
etc_isotype <- unique(as.character(df_etc_strain$isotype))

### Merge                  

indep_strain_info_geo <- indep_strain_info %>%
  dplyr::mutate(geo = ifelse(strain %in% hw_strain, "Hawaii", 
                             ifelse(strain %in% na_strain, "N. America", 
                                    ifelse(strain %in% au_strain, "Australia", 
                                           ifelse(strain %in% af_strain, "Africa", 
                                                  ifelse(strain %in% at_strain, "Atlantic",
                                                                       ifelse(strain %in% eu_strain, "Europe", 
                                                                              ifelse(strain %in% sa_strain, "S. America", 
                                                                                     ifelse(strain %in% nz_strain, "New Zealand", 
                                                                                                   ifelse(strain %in% as_strain, "Asia", "Unknown")))))))))) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(n_geo=n_distinct(geo)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(geo_isotype = ifelse(n_geo==1, geo, "Multiple"))

geo.colours <- c("Hawaii"="lightskyblue2", "Australia"="#0047ab", "N. America"="#C24CF6",
                 "Europe"="#FF0009","S. America"="maroon", "Africa"="lightpink2", "Atlantic"="turquoise", 
                 "New Zealand"="darkorange1", "Asia" = "black", "Unknown" = "burlywood3", "Multiple"="springgreen4")


eu.colours <- c("Greece"="lightskyblue2", "France"="#0047ab", "Czechia"="#C24CF6", "Germany"="yellow3",
                 "Spain"="#FF0009","Portugal"="maroon", "Switzerland"="lightpink2", "United Kingdom"="#9b351b", 
                 "Netherlands"="darkorange1", "Belgium" = "black", "Italy"="springgreen4")

## Big Island strain

plot_bi <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="gray51", fill="white", size=0.2, alpha = 0.8)+
  geom_point(data = df_SNV_geo, aes(x=long, y=lat, size=homo_alt/1e5, fill=geo), shape =21, alpha = 0.7) +
  geom_text_repel(data=dplyr::filter(df_SNV_geo, long < -154 & long > -156),
                  label=dplyr::filter(df_SNV_geo, long < -154 & long > -156)$strain, 
                  aes(x=long,y=lat), size = 3) +
  theme_map()+
  theme(legend.position='none')+
  scale_fill_manual(values = geo.colours) +
  scale_size_continuous(range = c(1,8), trans='exp', guide=FALSE) +
  ggsn::scalebar(data = dplyr::filter(indep_strain_info, !is.na(lat)), dist = 2000, dist_unit = "km", location = "bottomleft",
                 transform = TRUE, model = "WGS84", height = 0.02, st.size = 2.3, st.dist	= 0.03, border.size = 0.3)  +
  guides(fill= guide_legend(nrow=3)) +
  coord_cartesian(x=c(-157, -154), y=c(19,20.5)) 

plot_bi

df_bi_strain <- indep_strain_info %>%
  dplyr::filter(long > -157  & long < -154) %>%
  dplyr::filter(lat > 19 & lat < 20.5) 

bi_strain <- as.character(df_bi_strain$strain)
bi_isotype <- unique(as.character(df_bi_strain$isotype))

save(indep_strain_info_geo,geo.colours,eu.colours,bi_strain, file="Processed_Data/strain_geo.RData")


## EU strains - divide FR, DE, UK, NE


#8. Iberia

plot_ib <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="white", size=0.2)+
  geom_point(data = indep_strain_info, aes(x=long, y=lat), shape =21, size =3, alpha = 0.9, fill = 'red') +
  geom_label_repel(aes(long, lat, label = strain), data = indep_strain_info, fontface = 'bold', color = 'black', size = 2, box.padding = 0.25, point.padding = 0.3, segment.color = 'black', segment.alpha = 0.6, fill = 'red', nudge_y = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_map()+
  theme(text = element_text(size=12)) +
  lims(x=c(-10, 4), y=c(34,43)) 

plot_ib

df_ib_strain <- indep_strain_info %>%
  dplyr::filter(long > -10  & long < 4) %>%
  dplyr::filter(lat > 34 & lat < 43) 

ib_strain <- as.character(df_ib_strain$strain)
ib_isotype <- unique(as.character(df_ib_strain$isotype))


## FR strain

df_fr_strain <- indep_strain_info %>%
  dplyr::filter(country=="France")

fr_strain <- as.character(df_fr_strain$strain)
fr_isotype <- unique(as.character(df_fr_strain$isotype))

## DE strain

df_de_strain <- indep_strain_info %>%
  dplyr::filter(country=="Germany")

de_strain <- as.character(df_de_strain$strain)
de_isotype <- unique(as.character(df_de_strain$isotype))


