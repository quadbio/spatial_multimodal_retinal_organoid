# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 11.10.21

setwd('/home/charmel/4i_organoid_pipeline/retina/Rotator')
source('../utils.R')

df_intensity <- read_cluster_profiles(timepoints = c("week6","week12","week18","week24","week39","adult")) %>%
  unite('id',c('organoid','window'), remove = FALSE) %>%
  scale_cluster_profiles()

message('Loading meta files...')
df_meta <- readRDS('/local2/USERS/charmel/rotator_dataset_2/df_meta.rds')

selected_wedges <- df_meta %>% filter(!selected) %>% pull(id)

message('Filtering intensity profiles for selected windows...')
df_intensity <- df_intensity %>%
  filter(id %in% selected_wedges) %>%
  ungroup() %>%
  mutate(stain=str_c("cluster_",stain)) %>%
  select(id, organoid, window, radial_position, stain, intensity_scaled)

dir <- '/local2/USERS/charmel/rotator_dataset_2/'

features <- df_intensity %>% distinct(stain) %>% pull()

for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:', features[feature]))
  df_intensity %>%
    filter(stain == features[feature]) %>%
    write.csv(paste0(dir, features[feature], '.csv'))
}






