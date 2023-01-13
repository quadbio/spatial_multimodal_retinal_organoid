# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 25.01.22


setwd('laminator_analysis/retina_dev_nuclei_rgc')
source('../utils.R')

# load 4i nuclei datasets
list_4i_nuclei <- readRDS('/local2/USERS/charmel/integrated_seurat_4i_scRNA_scATAC_v3.rds') # TODO change path
names(list_4i_nuclei) <- c('adult','week12','week18','week24','week39','week6')


# protein stains
df_subset <- map(list_4i_nuclei, function (x){x$seu_4i@meta.data %>% mutate(organoid=as.character(organoid)) %>% as_tibble()}) %>%
  bind_rows() %>% filter(age=='week12', cell_type=='RGC') %>% mutate(x=round(X), y=round(Y))
subset <- df_subset %>% distinct(organoid) %>% pull()

df_intensity <- read_radial_profiles(subset=subset, df_subset_positions = df_subset) %>% z_scoring_age_sample()

dir <- 'data/processed/4i/laminator/dataset_nuclei/'

features <- df_intensity %>% distinct(gene) %>% pull()
for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:',features[feature]))
  df_intensity %>% filter(gene==features[feature]) %>% write.csv(paste0(dir,features[feature],'.csv'))
}

# for mtus
dir <- '/data/processed/4i/laminator/dataset_nuclei_mtu/'

df_intensity <- read_radial_cluster_profiles(subset=subset, df_subset_positions=df_subset) %>%
  scale_cluster_profiles()

df_intensity <- df_intensity %>%
  ungroup() %>%
  mutate(stain=str_c("cluster_",stain)) %>%
  group_by(stain) %>%
  mutate(intensity=rescale(intensity)) %>% ungroup()


features <- df_intensity %>% distinct(stain) %>% pull()

for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:', features[feature]))
  df_intensity %>%
    filter(stain == features[feature]) %>%
    write.csv(paste0(dir, features[feature], '.csv'))
}





