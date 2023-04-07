# load libraries
library(TSdist)
library(tidyverse)
library(furrr)
library(FBN)

# set working directory
setwd('data/processed/4i/laminator/dataset_fish')

# list files of feature datasets
files <- list.files(pattern='.csv')
names(files) <- str_remove(files, '.csv')

message("Setting up the future...")
future::plan('multisession', workers = 20)

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
  df <- df  %>% select(id, radial_position, intensity)  %>%
    pivot_wider(id_cols = id, values_from = intensity, names_from = "radial_position")
  M <- as.matrix(df %>% select(-id))
  rownames(M) <- df$id
  M <- t(apply(M, 1, smooth_and_down_sample_profile))
  M_dist <-  TSDatabaseDistances(M, distance=measure, diag=TRUE, upper=TRUE) %>% as.matrix()
  return(M_dist)
}

message("Calculating fourier distances with future...")
distance_matrices <- future_map(files,get_distances_by_feature, measure='fourier')
message("Saving fourier distance matrices...")
dir_results <- 'data/processed/fish/laminator/analysis_results/'
dir.create(dir_results)
saveRDS(distance_matrices, paste0(dir_results,'distance_matrices_fourier_fish.rds'))

