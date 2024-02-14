##########################
# Highly diverse Symbiodiniaceae communities hosted by corals in a
# global hotspot of marine biodiversity
#
# Prepared by Ming Sheng Ng - 20/12/2023
##########################

library(ggplot2); packageVersion('ggplot2') # 3.4.4
library(dplyr); packageVersion('dplyr') # 1.1.3
library(stringr); packageVersion('stringr') # 1.5.0
library(reshape2); packageVersion('reshape2') # 1.4.4
library(vegan); packageVersion('vegan') # 2.6.4
library(pairwiseAdonis); packageVersion('pairwiseAdonis') # 0.4.1
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Code for reproducing figures and results

##### fig 2 Symbiodiniaceae genera composition ######
clade_table <- read.table('./datafiles/clade_profiles.txt', header = T, sep = '\t')
comp <- clade_table[4:7] %>%
  apply(1,function(x){x/sum(x) * 100}) %>% 
  t() %>% as.data.frame()

clade_composition <- cbind(clade_table[1:3],comp) %>% melt()

clade_composition <- clade_composition[order(clade_composition$Species),]

pal <- RColorBrewer::brewer.pal(4, 'Set1')

fig2 <- ggplot(aes(x = SampleID, y = value, fill = variable), data = clade_composition) +
  geom_bar(position = 'stack', stat='identity') +
  facet_wrap(Species ~ Site, scales = 'free_x', nrow = 4, ncol = 6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y = 'Relative abudance (%)') + 
  scale_fill_manual(values = pal)

###################
##### fig 3a Symbiodiniaceae type profile composition #####
otu_table <- readxl::read_xlsx('./datafiles/type_profiles.xlsx') %>% as.data.frame()
otu_comp <- otu_table[4:ncol(otu_table)] %>%
  apply(1,function(x){x/sum(x) * 100}) %>% 
  t() %>% as.data.frame()

count_prev <- function(x) { # custom function to order by prevalence
  ifelse(x>0,1,0) %>% sum()
}

otu_melt <- cbind(otu_table[1:3], otu_comp) %>% melt()
# order by prevalence
otu_melt$variable <- with(otu_melt, reorder(variable, value, FUN = count_prev, decreasing = TRUE))

source('./iwanthue_palettes.R') # palette of colours generated from https://medialab.github.io/iwanthue/

fig3a <- ggplot(data = otu_melt, aes(x = SampleID, y = value, fill = variable)) +
  geom_bar(position = 'stack', stat = 'identity', colour = 'black') +
  facet_wrap(Species ~ Site, scales ="free_x", nrow = 4, ncol = 6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(legend.position = 'None') + #remove legends
  labs(y = 'Relative abundance (%)') +
  scale_fill_manual(values = version2)

##################

##### fig 3b Symbiodiniaceae type profile composition in majority sequence #####
majority <- read.table('./datafiles/majority sequences.txt', sep = '\t', header = TRUE, check.names = FALSE)

majority_comp <- majority[7:ncol(majority)] %>%
  apply(1,function(x){x/sum(x) * 100}) %>% 
  t() %>% as.data.frame()

maj <- cbind(majority[2:5], majority_comp[2:ncol(majority_comp)]) %>% melt()
maj$variable <- as.character(maj$variable)

#changes duplicate names back - e.g. C15.1 to C15
for (i in 1:nrow(maj)){
  if (str_detect(maj[i,"variable"], '\\.')){
    var = str_split(maj[i,"variable"], '\\.')[[1]][1]
    maj[i,"variable"] = var
    print(var) #to observe progress - remove to prevent printing of progress
  }
}

# order by prevalence
maj$variable <- with(maj, reorder(variable, value, FUN = count_prev, decreasing = TRUE))

fig3b <- ggplot(data = maj, aes(x = SampleID, y = value, fill = variable)) +
  geom_bar(position = 'stack', stat = 'identity', colour = 'black') +
  facet_wrap(Species ~ Site, scales ="free_x", nrow = 4, ncol = 6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = version3) +
  labs(y= 'Relative abundance (%)') +
  theme(legend.position = 'none') # remove legends

##### 1. sequence beta diversity #####
##### between-types + figure 4a #####
##### cladocopium #####
c_unifrac_pcoa <- read.csv('./datafiles/cladocopium_unifrac_sample_PCoA.csv')
type_profiles <- readxl::read_xlsx('./datafiles/type_profiles_raw.xlsx')

# remove last row
c_proportion_explained <- c_unifrac_pcoa[nrow(c_unifrac_pcoa), 3:ncol(c_unifrac_pcoa)]
c_unifrac_pcoa <- c_unifrac_pcoa[-nrow(c_unifrac_pcoa),]

c_unifrac_pcoa$Species <- rep('coral', nrow(c_unifrac_pcoa))
c_unifrac_pcoa$Site <- rep('site', nrow(c_unifrac_pcoa))

# add species and site metadata 
for (row in 1:nrow(c_unifrac_pcoa)){
  uid = c_unifrac_pcoa[row,]$sample_uid
  c_unifrac_pcoa[row, 'Species'] <- type_profiles[type_profiles$`ITS2 type profile` == uid, 'Species']
  c_unifrac_pcoa[row, 'Site'] <- type_profiles[type_profiles$`ITS2 type profile` == uid, 'Site']
}

fig4a_clad <- ggplot(data = c_unifrac_pcoa, aes(x = PC1, y = PC2, colour = Species)) +
  geom_point(size = 2) + 
  theme_bw() # fig 4a


##### durusdinium #####
d_unifrac_pcoa <- read.csv('./datafiles/durusdinium_unifrac_sample_PCoA.csv')

# remove last row
d_proportion_explained <- d_unifrac_pcoa[nrow(d_unifrac_pcoa), 3:ncol(d_unifrac_pcoa)]
d_unifrac_pcoa <- d_unifrac_pcoa[-nrow(d_unifrac_pcoa),]

d_unifrac_pcoa$Species <- rep('coral', nrow(d_unifrac_pcoa))
d_unifrac_pcoa$Site <- rep('site', nrow(d_unifrac_pcoa))

for (row in 1:nrow(d_unifrac_pcoa)){
  uid = d_unifrac_pcoa[row,]$sample_uid
  d_unifrac_pcoa[row, 'Species'] <- type_profiles[type_profiles$`ITS2 type profile` == uid, 'Species']
  d_unifrac_pcoa[row, 'Site'] <- type_profiles[type_profiles$`ITS2 type profile` == uid, 'Site']
}

fig4a_durus <- ggplot(data = d_unifrac_pcoa, aes(x = PC1, y = PC2, colour = Species)) +
  geom_point(size = 2) + 
  theme_bw() 


##### between-sample + figure 4b ######
cladocopium_types <- read.csv('./datafiles/unifrac_profiles_PCoA_coords_C.csv')
ctype_proportion_explained <- cladocopium_types[nrow(cladocopium_types),]

cladocopium_types <- cladocopium_types[-nrow(cladocopium_types),]

fig4b_clad <- ggplot(data = cladocopium_types, aes(x = PC1, y = PC2, label = dominant)) + 
  geom_point(aes(colour = dominant)) + 
  geom_text() +
  labs(colour = 'Majority sequence') + 
  theme_bw()

durusdinium_types <- read.csv('./datafiles/unifrac_profiles_PCoA_coords_D.csv')
dtype_proportion_explained <- durusdinium_types[nrow(durusdinium_types),]

durusdinium_types <- durusdinium_types[-nrow(durusdinium_types),]

fig4b_durus <- ggplot(data = durusdinium_types, aes(x = PC1, y = PC2, label = dominant)) + 
  geom_point(aes(colour = dominant)) + 
  geom_text() +
  labs(colour = 'Majority sequence') + 
  theme_bw()


##### 2. ecological beta diversity #####
types_table <- read.csv('./datafiles/comb_abs_abund.csv', check.names = FALSE, header = T)
row.names(types_table) = otu_table$sampleID

## form robust.aitchison distance matrix
types <- types_table[,-c(1,2,3,4)]
comb_ait <- vegdist(types, 'robust.aitchison')

## is host a significant driver of symb community?
# PERMANOVA (combined)
table2 <- adonis2(comb_ait ~ types_table$Species*otu_table$Site)

# overall PCA
comb_pca <- prcomp(comb_ait)
ordiplot(comb_pca,)
summary(comb_pca)

comb_x <- comb_pca$x %>% as.data.frame() %>% select('PC1','PC2')
comb_x$species <- otu_table$Species

# PCA with Robust Aitchison distances (fig 5a)
fig5a <- ggplot(comb_x) +
  geom_point(aes(x = PC1, y = PC2, colour = species, pch = species), size = 2) +
  scale_shape_manual(values = c(0,1,2,4)) + 
  stat_ellipse(aes(x = PC1, y = PC2, colour = species), lty = 'dashed', linewidth = 0.7) +
  theme_bw() + theme(legend.position = c(.87,.87)) + 
  xlab('PC1 (69.22%)') + ylab('PC2 (20.31%)')



##### PERMANOVA (separated by host species) ######
# Since species was significant in symbd structuring - look at individual host species

##### diplo PCA, permanova, betadisper #####

metadata <- readxl::read_xlsx('./datafiles/comb_metadata_permanova.xlsx')
diplo_rows <- metadata$species == 'Diploastrea heliopora'
diplo_meta <- metadata %>% filter(diplo_rows)
diplo_meta$sampleID <- diplo_meta$sampleID %>% str_replace('Diplo_','')

diplo_otu <- readxl::read_xlsx('./datafiles/diplo_seq_absabund.xlsx') %>% as.data.frame()
row.names(diplo_otu) <- diplo_otu[,1]
diplo_otu <- diplo_otu[-1]

#ensuring that sampleIDs are the same
diplo_meta$sampleID == row.names(diplo_otu)

diplo_ait <- vegdist(diplo_otu[-1], 'robust.aitchison')

diplo_pca <- prcomp(diplo_ait, center = T)
biplot(diplo_pca, )
summary(diplo_pca)

diplo_envfit <- envfit(diplo_pca ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data=diplo_meta,perm=999)

diplo.scores <- as.data.frame(scores(diplo_envfit, display = 'vectors'))
diplo.scores <- cbind(diplo.scores, Env = rownames(diplo.scores))

en_coord_cont = as.data.frame(scores(diplo_envfit, 'vectors')) * ordiArrowMul(diplo_envfit)
en_coord_cat = as.data.frame(scores(diplo_envfit, 'factors')) * ordiArrowMul(diplo_envfit)

diplo_x <- diplo_pca$x %>% as.data.frame() %>% select('PC1','PC2')
diplo_x$site <- diplo_meta$collection_island

fig5b_diplo <- ggplot(diplo_x) +
  geom_point(aes(x = PC1, y = PC2, colour = site)) +
  stat_ellipse(aes(x = PC1, y = PC2, colour = site), lty = 2) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               data = en_coord_cont, colour = "#838383", arrow = arrow(length = unit(.3, 'cm'))) +
  geom_text(data = en_coord_cont, aes(x = PC1 + 0.05, y = PC2 + 0.10), colour = "#838383", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme_bw() + 
  theme(legend.position = c(0.89,0.85)) + 
  xlab('PC1 (70.1%)') + ylab('PC2 (9.98%)') + ggtitle('Diploastrea heliopora') 


diplo_permanova_env <- adonis2(as.dist(diplo_ait) ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data = diplo_meta)
diplo_permanova_env

diplo_betadisper <- betadisper(diplo_ait, group = diplo_meta$collection_island)
anova(diplo_betadisper)
plot(diplo_betadisper)

diplo_permanova_site <- adonis2(as.dist(diplo_ait) ~ diplo_meta$collection_island)
diplo_permanova_site
pairwise.adonis(as.dist(diplo_ait), diplo_meta$collection_island)


##### plutea PCA, permanova, betadisper #####

plutea_rows <- metadata$species == 'Pocillopora lutea'
plutea_meta <- metadata %>% filter(plutea_rows)
plutea_meta$sampleID <- plutea_meta$sampleID %>% str_replace('Plutea_','')

plutea_otu <- readxl::read_xlsx('./datafiles/plutea_seq_absabund.xlsx') %>% as.data.frame()
row.names(plutea_otu) <- plutea_otu$sampleID
plutea_otu <- plutea_otu[,-c(1,2,3)] #remove sample_uid and sampleID

plutea_otu <- plutea_otu[-50,]
plutea_meta <- plutea_meta[-50,]

plutea_meta$sampleID == row.names(plutea_otu)

plutea_ait <- vegdist(plutea_otu, 'robust.aitchison')

plutea_pca <- prcomp(plutea_ait, center = T)
biplot(plutea_pca)
summary(plutea_pca)

plutea_envfit <- envfit(plutea_pca ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data=plutea_meta,perm=999)

plutea.scores <- as.data.frame(scores(plutea_envfit, display = 'vectors'))
plutea.scores <- cbind(plutea.scores, Env = rownames(plutea.scores))

en_coord_cont = as.data.frame(scores(plutea_envfit, 'vectors')) * ordiArrowMul(plutea_envfit)
en_coord_cat = as.data.frame(scores(plutea_envfit, 'factors')) * ordiArrowMul(plutea_envfit)

plutea_x <- plutea_pca$x %>% as.data.frame() %>% select('PC1','PC2')
plutea_x$site <- plutea_meta$collection_island

fig5b_plutea <- ggplot(plutea_x) +
  geom_point(aes(x = PC1, y = PC2, colour = site)) +
  stat_ellipse(aes(x = PC1, y = PC2, colour = site)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               data = en_coord_cont, colour = "#838383") +
  geom_text(data = en_coord_cont, aes(x = PC1 + 0.05, y = PC2 + 0.10), colour = "#838383", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme_bw() + 
  theme(legend.position = 'none') +
  xlab('PC1 (87.8%)') + ylab('PC2 (1.85%)') + ggtitle('Porites lutea')


plutea_permanova_env <- adonis2(as.dist(plutea_ait) ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data = plutea_meta)
plutea_permanova_env

plutea_betadisper <- betadisper(plutea_ait, group = plutea_meta$collection_island)
anova(plutea_betadisper)
plot(plutea_betadisper)
permutest(plutea_betadisper)
TukeyHSD(plutea_betadisper)

plutea_permanova_site <- adonis2(as.dist(plutea_ait) ~ plutea_meta$collection_island)
plutea_permanova_site
pairwise.adonis(as.dist(plutea_ait), plutea_meta$collection_island)


##### pacuta PCA, permanova, betadisper #####

pacuta_rows <- metadata$species == 'Pocillopora acuta'
pacuta_meta <- metadata %>% filter(pacuta_rows)
pacuta_meta$sampleID <- pacuta_meta$sampleID %>% str_replace('Pacuta_','')

pacuta_otu <- readxl::read_xlsx('./datafiles/pacuta_seq_absabund.xlsx') %>% as.data.frame()
row.names(pacuta_otu) <- pacuta_otu$sampleID
pacuta_otu <- pacuta_otu[,-c(1,2,3)] #remove sample_uid, sampleID, and Type

pacuta_meta$sampleID == row.names(pacuta_otu)

pacuta_ait <- vegdist(pacuta_otu, 'robust.aitchison')

pacuta_pca <- prcomp(pacuta_ait, center = T) #transformed to aitchison so no need to scale
summary(pacuta_pca)

pacuta_envfit <- envfit(pacuta_pca ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data=pacuta_meta,perm=999)

ordiplot(pacuta_pca)
plot(pacuta_envfit)

pacuta.scores <- as.data.frame(scores(pacuta_envfit, display = 'vectors'))
pacuta.scores <- cbind(pacuta.scores, Env = rownames(pacuta.scores))

en_coord_cont = as.data.frame(scores(pacuta_envfit, 'vectors')) * ordiArrowMul(pacuta_envfit)
en_coord_cat = as.data.frame(scores(pacuta_envfit, 'factors')) * ordiArrowMul(pacuta_envfit)

pacuta_x <- pacuta_pca$x %>% as.data.frame() %>% select('PC1','PC2')
pacuta_x$site <- pacuta_meta$collection_island

fig5b_pacuta <- ggplot(pacuta_x) +
  geom_point(aes(x = PC1, y = PC2, colour = site)) +
  stat_ellipse(aes(x = PC1, y = PC2, colour = site)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               data = en_coord_cont, colour = "#838383") +
  geom_text(data = en_coord_cont, aes(x = PC1 + 0.05, y = PC2 + 0.10), colour = "#838383", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme_bw() + 
  theme(legend.position = 'none') +
  xlab('PC1 (61.8%)') + ylab('PC2 (14.0%)') + ggtitle('Pocillopora acuta')


pacuta_permanova_env <- adonis2(as.dist(pacuta_ait) ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data = pacuta_meta)
pacuta_permanova_env

pacuta_betadisper <- betadisper(pacuta_ait, group = pacuta_meta$collection_island)
anova(pacuta_betadisper)
plot(pacuta_betadisper)
TukeyHSD(pacuta_betadisper)

pacuta_permanova_site <- adonis2(as.dist(pacuta_ait) ~ pacuta_meta$collection_island)
pacuta_permanova_site
pairwise.adonis(as.dist(pacuta_ait), pacuta_meta$collection_island)


##### pachy PCA, permanova, betadisper #####

pachy_rows <- metadata$species == 'Pachyseris speciosa'
pachy_meta <- metadata %>% filter(pachy_rows)
pachy_meta$sampleID <- pachy_meta$sampleID %>% str_replace('Pachyseris_','')

pachy_otu <- readxl::read_xlsx('./datafiles/pachy_seq_absabund.xlsx') %>% as.data.frame()
row.names(pachy_otu) <- pachy_otu$sampleID
pachy_otu <- pachy_otu[,-c(1,2,3)] #remove sample_uid, sampleID, and Type

pachy_meta$sampleID == row.names(pachy_otu)

pachy_ait <- vegdist(pachy_otu, 'robust.aitchison')

pachy_pca <- prcomp(pachy_ait, center = T) #transformed to aitchison so no need to scale
summary(pachy_pca)

pachy_envfit <- envfit(pachy_pca ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data=pachy_meta,perm=999)

ordiplot(pachy_pca)
plot(pachy_envfit)

pachy.scores <- as.data.frame(scores(pachy_envfit, display = 'vectors'))
pachy.scores <- cbind(pachy.scores, Env = rownames(pachy.scores))

en_coord_cont = as.data.frame(scores(pachy_envfit, 'vectors')) * ordiArrowMul(pachy_envfit)
en_coord_cat = as.data.frame(scores(pachy_envfit, 'factors')) * ordiArrowMul(pachy_envfit)

pachy_x <- pachy_pca$x %>% as.data.frame() %>% select('PC1','PC2')
pachy_x$site <- pachy_meta$collection_island

fig5b_pachy <- ggplot(pachy_x) +
  geom_point(aes(x = PC1, y = PC2, colour = site)) +
  stat_ellipse(aes(x = PC1, y = PC2, colour = site)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               data = en_coord_cont, colour = "#838383") +
  geom_text(data = en_coord_cont, aes(x = PC1 + 0.05, y = PC2 + 0.10), colour = "#838383", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme_bw() + 
  theme(legend.position = 'none') +
  xlab('PC1 (59.7%)') + ylab('PC2 (24.3%)') + ggtitle('Pachyseris speciosa')

pachy_pcaplot

pachy_permanova_env <- adonis2(as.dist(pachy_ait) ~ BO_parmean + BO_ph + BO21_tempmean_ss + BO21_salinitymean_ss, data = pachy_meta)
pachy_permanova_env

pachy_betadisper <- betadisper(pachy_ait, group = pachy_meta$collection_island)
anova(pachy_betadisper)
plot(pachy_betadisper)
TukeyHSD(pachy_betadisper)

pachy_permanova_site <- adonis2(as.dist(pachy_ait) ~ pachy_meta$collection_island)
pachy_permanova_site
pairwise.adonis(as.dist(pachy_ait), pachy_meta$collection_island)
