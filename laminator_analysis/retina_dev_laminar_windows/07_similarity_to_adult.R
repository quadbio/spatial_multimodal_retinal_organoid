
library(abind)
library(DistatisR)
library(tidyverse)

# load datasets
df_pseudotime <- readRDS('data/processed/4i/laminator/analysis_results/laminar_window_dpt_trajectory.rds')  %>%
  mutate(dpt_rank=max(dpt_rank)-dpt_rank)
wedge_ids <- df_pseudotime %>% distinct(id) %>% pull()


message('Loading distance matrices and reshaping array...')
feature_cube <- readRDS('data/processed/4i/laminator/analysis_results/distance_matrices_fourier_protein.rds') %>%
  abind(along=3)
feature_cube <- feature_cube[wedge_ids, wedge_ids,]

message('Running distatis...')
distatis_res <- distatis(feature_cube)

message('Saving distatis results...')
saveRDS(distatis_res, 'data/processed/4i/laminator/analysis_results/distatis_from_fourier.rds')


# similiarity to adults
wedge_ids <- df_pseudotime %>% distinct(id) %>% pull()
wedge_ids_adult <- df_pseudotime %>% filter(type=='Adult') %>% distinct(id)  %>% pull(id)
wedge_ids_end <- df_pseudotime %>% filter(dpt_rank >= quantile(dpt_rank,.95)) %>% distinct(id) %>% pull(id)
wedge_ids_start <- df_pseudotime %>% filter(dpt_rank <= quantile(dpt_rank,.05)) %>% distinct(id) %>% pull(id)

calculate_maturation_score <- function(wedge_ids_start, wedge_ids_end, distatis_res, df_pseudotime){
  df_start <-  apply(distatis_res$res4Splus$Splus[wedge_ids,wedge_ids_start], 1, mean) %>%
    enframe() %>%
    rename(id=name, distance=value) %>%
    mutate(maturation_point='start') %>%
    left_join(df_pseudotime %>%
                select(id, dpt_rank, type), by='id')

  df_end <-  apply(distatis_res$res4Splus$Splus[wedge_ids,wedge_ids_end], 1, mean) %>%
    enframe() %>%
    rename(id=name, distance=value) %>% mutate(maturation_point='end') %>%
    left_join(df_pseudotime %>% select(id, dpt_rank, type), by='id')

  df_maturation <- bind_rows(df_start, df_end) %>%
    group_by(maturation_point) %>%
    mutate(distance = (distance-mean(distance))/sd(distance)) %>%
     ungroup() %>%
     pivot_wider(id_cols = c('id','type','dpt_rank'), values_from = 'distance', names_from = 'maturation_point') %>%
     mutate(score=rescale(end-start, to=c(-1,1)))
  return(df_maturation)
}

df_maturation <- calculate_maturation_score(wedge_ids_start, wedge_ids_end, distatis_res, df_pseudotime)
saveRDS(distatis_res, 'data/processed/4i/laminator/analysis_results/df_maturation_score.rds')



