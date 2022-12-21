# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 15.10.21

suppressPackageStartupMessages({
  library(TSdist)
  library(tidyverse)
  library(furrr)
  library(FBN)
})

setwd('/local2/USERS/charmel/rotator_dataset_2/')

# list files of feature datasets
files <- list.files(pattern='.csv')
files <- files[str_detect(files, 'cluster')]
names(files) <- str_remove(files, '.csv')

message("Setting up the future...")
future::plan('multisession', workers = length(files))


smooth_and_down_sample_profile <- function(x, window_size=20, sampling_factor=2){
  x <- meanFilter(x, windowSize = window_size)
  t <- 1:length(x)
  f <- approxfun(t,x)
  x_sampled <- f(seq(from=1,to=length(x),by=sampling_factor))
  return(x_sampled)
}


# define distance function
get_distances_by_feature <- function(file, measure='fourier') {
  df <- read_csv(file)
  df <- df  %>% select(id, radial_position, intensity_scaled)  %>%
    pivot_wider(id_cols = id, values_from = intensity_scaled, names_from = "radial_position")
  M <- as.matrix(df %>% select(-id))
  rownames(M) <- df$id
  M <- t(apply(M, 1, smooth_and_down_sample_profile))
  M_dist <-  TSDatabaseDistances(M, distance=measure, diag=TRUE, upper=TRUE) %>% as.matrix()
  return(M_dist)
}

message("Calculating fourier distances with future...")
distance_matrices <- future_map(files,get_distances_by_feature, measure='fourier')
message("Saving fourier distance matrices...")
saveRDS(distance_matrices, '/local2/USERS/charmel/rotator_analysis_results/distance_matrices_cluster_fourier_2.rds')

#message("Calculating cross correlation distances with future...")
#distance_matrices <- future_map(files,get_distances_by_feature, measure='ccor')
#message("Saving cross correlation distance matrices...")
#saveRDS(distance_matrices, '/local2/USERS/charmel/rotator_analysis_results/distance_matrices_ccor_2.rds')