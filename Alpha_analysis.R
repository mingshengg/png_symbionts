#=======================#
# Investigate if there are any sig diff in diversity metrics across sites and species
# Last updated: May 2023
#=======================#

library(dplyr)
library(vegan)
library(car)
library(ggplot2)
library(reshape2)

se <- function (x) { sd(x)/sqrt(length(x))}

otu_table <- readxl::read_xlsx('./type_profiles.xlsx') %>% as.data.frame()
rownames(otu_table) <- otu_table$SampleID
metadata <- readxl::read_xlsx('./comb_metadata_permanova.xlsx')

# occurrence across the OTUs
otu <- otu_table[4:ncol(otu_table)]
otu_presence <- ifelse(otu > 0, 1, 0)
occurrence <- data.frame(occurrence = colSums(otu_presence))

otu_table %>% 
  select(`C66/C3/C91-C1-1667_C-1748_C`) %>%
  filter(`C66/C3/C91-C1-1667_C-1748_C`>0)

ggplot(occurrence, aes(x = occurrence)) +
  geom_histogram() +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  xlab('Occurrence') +
  ylab('Frequency')

# count the number of types that only appeared  in a single sample
count_1 <- ifelse(occurrence == 1, 1, 0)
sum(count_1)

count_2 <- ifelse(occurrence == 2, 1, 0)
sum(count_2)

count_less_than_5 <- ifelse(occurrence <= 5, 1, 0)
sum(count_less_than_5)

## alpha diversity of indiv species
occ_table <- cbind(otu_table[1:3], otu_presence)
#pacuta
pacuta <- occ_table[occ_table$Species == 'Pocillopora acuta',]
nrow(pacuta)
pacuta[4:ncol(pacuta)] %>% rowSums() %>% mean()
pacuta[4:ncol(pacuta)] %>% rowSums() %>% se()

pacuta[4:ncol(pacuta)] %>% colSums() %>% sort(decreasing = TRUE)

#diplo
diplo <- occ_table[occ_table$Species == 'Diploastrea heliopora',]
nrow(diplo)
diplo[4:ncol(diplo)] %>% rowSums() %>% mean()
diplo[4:ncol(diplo)] %>% rowSums() %>% se()

diplo[4:ncol(diplo)] %>% colSums() %>% sort(decreasing = TRUE)

#pachy
pachy <- occ_table[occ_table$Species == 'Pachyseris speciosa',]
nrow(pachy)
pachy[4:ncol(pachy)] %>% rowSums() %>% mean()
pachy[4:ncol(pachy)] %>% rowSums() %>% se()

pachy[4:ncol(pachy)] %>% colSums() %>% sort(decreasing = TRUE)

#plutea
plutea <- occ_table[occ_table$Species == 'Porites lutea',]
nrow(plutea)
plutea[4:ncol(plutea)] %>% rowSums() %>% mean()
plutea[4:ncol(plutea)] %>% rowSums() %>% se()

plutea[4:ncol(plutea)] %>% colSums() %>% sort(decreasing = TRUE)

# how many sampled individuals had only one type profile found?
temp <- occ_table[4:ncol(occ_table)] %>% rowSums()
no_profile <- cbind(occ_table[1:3], temp)
no_profile %>% filter(Species == 'Diploastrea heliopora') %>% 
  filter(temp == 1) %>% nrow()

no_profile %>% filter(Species == 'Pachyseris speciosa') %>% 
  filter(temp == 1) %>% nrow()

no_profile %>% filter(Species == 'Pocillopora acuta') %>% 
  filter(temp == 1) %>% nrow()

no_profile %>% filter(Species == 'Porites lutea') %>% 
  filter(temp == 1) %>% nrow()

## number of type profiles per species and site
occ_table %>% melt() %>% filter(value > 0) %>%
  group_by(Species, Site) %>%
  summarise(avg = length(unique(variable))) %>%
  print(n = 24)

library(microbiome)
# to investigate if type profiles were found across samples within each coral host and site
species = c(unique(occ_table$Species))
site = unique(occ_table$Site)
for (i in species){
  for (j in site){
    print(c(i,j))
    temp <- occ_table %>% melt() %>% filter(value >0) %>%
      filter(Species == i & Site == j) %>% count(variable) 
    rows = occ_table %>% melt() %>% filter(value >0) %>%
      filter(Species == i & Site == j) %>% nrow()
    print(cbind(temp,rows))
  }
}


p_speciosa_profiles <- occ_table %>% melt() %>% filter(Species == 'Pachyseris speciosa') %>% filter(value > 0) %>% .$variable %>% unique()
p_lutea_profiles <- occ_table %>% melt() %>% filter(Species == 'Porites lutea') %>% filter(value > 0) %>% .$variable %>% unique()
p_acuta_profiles <- occ_table %>% melt() %>% filter(Species == 'Pocillopora acuta') %>% filter(value > 0) %>% .$variable %>% unique()
d_heliopora_profiles <- occ_table %>% melt() %>% filter(Species == 'Diploastrea heliopora') %>% filter(value > 0) %>% .$variable %>% unique()

p_speciosa_unique <- p_speciosa_profiles[!(p_speciosa_profiles %in% p_lutea_profiles) & 
                                           !(p_speciosa_profiles %in% p_acuta_profiles) &
                                           !(p_speciosa_profiles %in% d_heliopora_profiles)]

p_lutea_unique <- p_lutea_profiles[!(p_lutea_profiles %in% p_speciosa_profiles) & 
                                           !(p_lutea_profiles %in% p_acuta_profiles) &
                                           !(p_lutea_profiles %in% d_heliopora_profiles)]

p_acuta_unique <- p_acuta_profiles[!(p_acuta_profiles %in% p_lutea_profiles) & 
                                           !(p_acuta_profiles %in% p_speciosa_profiles) &
                                           !(p_acuta_profiles %in% d_heliopora_profiles)]

d_helipora_unique <- d_heliopora_profiles[!(d_heliopora_profiles %in% p_lutea_profiles) & 
                                           !(d_heliopora_profiles %in% p_acuta_profiles) &
                                           !(d_heliopora_profiles %in% p_speciosa_profiles)]


species_site <- occ_table %>% melt() %>% filter(Species == 'Diploastrea heliopora') %>%
  filter(value > 0) %>% group_by(variable) %>%
  count(Site)

# no. of type profiles further exclusive to a single site
species_site %>% filter(variable %in% d_helipora_unique) %>% count(variable) %>% filter(n == 1) %>% nrow()
