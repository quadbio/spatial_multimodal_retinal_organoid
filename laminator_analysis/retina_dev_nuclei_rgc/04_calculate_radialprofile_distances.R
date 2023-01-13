suppressPackageStartupMessages({
  library(TSdist)
  library(tidyverse)
  library(furrr)
  library(FBN)
})

setwd('data/processed/4i/laminator/dataset_nuclei/')

# list files of feature datasets
files <- list.files(pattern='.csv')
names(files) <- str_remove(files, '.csv')

# set result direcetory
dir_results <-  'data/processed/4i/laminator/analysis_results/'
dir.create(dir_results)

message("Setting up the future...")
future::plan('multisession', workers = length(files))


smooth_and_down_sample_profile <- function(x, window_size=10, sampling_factor=2){
  x <- meanFilter(x, windowSize = window_size)
  t <- 1:length(x)
  f <- approxfun(t,x)
  x_sampled <- f(seq(from=1,to=length(x),by=sampling_factor))
  return(x_sampled)
}


# define distance function
get_distances_by_feature <- function(file, measure='fourier', max_radius=500) {
  df <- read_csv(file, show_col_types = FALSE)
  df <- df  %>% select(id, position, intensity_scaled)  %>% filter(position <= max_radius) %>%
    pivot_wider(id_cols = id, values_from = intensity_scaled, names_from = "position")
  M <- as.matrix(df %>% select(-id))
  rownames(M) <- df$id
  M <- t(apply(M, 1, smooth_and_down_sample_profile))
  M_dist <-  TSDatabaseDistances(M, distance=measure, diag=TRUE, upper=TRUE) %>% as.matrix()
  return(M_dist)
}


calculate_distances <- function(max_radius, files){
  message(paste0("Calculating fourier distances with future for maximal radius of ",max_radius," ..."))
  distance_matrices <- future_map(files,get_distances_by_feature, measure='fourier', max_radius=max_radius)
  message("Saving fourier distance matrices...")
  path <- paste0(dir_results,'distance_matrices_fourier_nuclei_',max_radius,'.rds')
  saveRDS(distance_matrices, path)
  message(paste('Saved to:',path))
}

radial_profile_sizes <- c(20,40,60,80,100,250,500)
map(radial_profile_sizes, calculate_distances, files=files)
