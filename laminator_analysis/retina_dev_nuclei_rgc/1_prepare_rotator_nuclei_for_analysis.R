# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 25.01.22


setwd('/home/charmel/4i_organoid_pipeline/retina/Rotator/Analysis')
source('../utils.R')

# load 4i nuclei datasets
list_4i_nuclei <- readRDS('/local2/USERS/charmel/integrated_seurat_4i_scRNA_scATAC_v3.rds')
names(list_4i_nuclei) <- c('adult','week12','week18','week24','week39','week6')

df_subset <- map(list_4i_nuclei, function (x){x$seu_4i@meta.data %>% mutate(organoid=as.character(organoid)) %>% as_tibble()}) %>%
  bind_rows() %>% filter(age=='week12', cell_type=='RGC') %>% mutate(x=round(X), y=round(Y))
subset <- df_subset %>% distinct(organoid) %>% pull()

df_intensity <- read_radial_profiles(subset=subset, df_subset_positions = df_subset) %>% z_scoring_age_sample()

dir <- '/local2/USERS/charmel/rotator_dataset_nuclei/'

features <- df_intensity %>% distinct(gene) %>% pull()
for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:',features[feature]))
  df_intensity %>% filter(gene==features[feature]) %>% write.csv(paste0(dir,features[feature],'.csv'))
}