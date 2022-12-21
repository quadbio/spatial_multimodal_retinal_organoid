# Title     : run distatis on distance matrices
# Objective : unifies distance matrices with distatis resulting into an compromised distance matrix of features and observations
# Created by: harmelc
# Created on: 16.09.21

library(abind)
library(DistatisR)
library(tidyverse)


# load datasets
df_pseudotime <- readRDS('/local2/USERS/charmel/pseudotime_trajectory_final/pseudotime_trajectory_v2.rds')  %>%
  mutate(dpt_rank=max(dpt_rank)-dpt_rank)
wedge_ids <- df_pseudotime %>% distinct(id) %>% pull()


message('Loading distance matrices and reshaping array...')
feature_cube <- readRDS('/local2/USERS/charmel/rotator_analysis_results/distance_matrices_fourier_2.rds') %>%
  abind(along=3)
feature_cube <- feature_cube[wedge_ids, wedge_ids,]

message('Running distatis...')
distatis_res <- distatis(feature_cube)

message('Saving distatis results...')
saveRDS(distatis_res, '/local2/USERS/charmel/rotator_analysis_results/distatis_from_fourier_2_v2.rds')
