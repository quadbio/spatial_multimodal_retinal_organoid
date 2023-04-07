# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 13.09.21

setwd('/home/charmel/4i_organoid_pipeline/retina/Rotator/Analysis')
source('../utils.R')


# protein stains
df_intensity <- read_intensity_profiles(timepoints = c("week6","week12","week18","week24","week39", "adult")) %>%
  unite('id',c('organoid','window'),remove=FALSE) %>% z_scoring_age_sample()

message('Loading meta files and filtering windows that originate from cut areas in the images...')
df_meta <- load_meta_files() %>%
  unite('id',c('organoid','window'),remove=FALSE) %>%
  group_by(organoid, contour) %>%
  filter(n()>12,
         !(organoid=="72"&contour==1),
         !(organoid =="72"&(window <10))) %>%
  detect_cut_areas() %>%
  mutate(selected=ifelse(organoid %in% c("43","6","39"), FALSE, selected))

df_meta <- bind_rows(df_meta %>% filter(organoid != "48"),
                     df_meta %>% filter(organoid=="48") %>% mutate(selected = !(window %in% 17:167)))

selected_wedges <- df_meta %>% filter(!selected) %>% pull(id)

message('Filtering intensity profiles for selected windows...')
df_intensity <- df_intensity %>%
  filter(id %in% selected_wedges) %>%
  ungroup() %>% rename(stain=gene) %>%
  select(id, organoid, window, radial_position, stain, intensity_scaled)

dir <- '/local2/USERS/charmel/rotator_dataset_2/'

saveRDS(df_meta,paste0(dir,'df_meta.rds'))

features <- df_intensity %>% distinct(stain) %>% pull()
for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:',features[feature]))
  df_intensity %>% filter(stain==features[feature]) %>% write.csv(paste0(dir,features[feature],'.csv'))
}

# mtu clusters

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




