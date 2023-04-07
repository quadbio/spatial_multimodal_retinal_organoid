setwd('/home/charmel/4i_organoid_pipeline/retina/Rotator/Analysis')
source('../utils.R')

# Load libraries
library(tidyverse)
library(destiny)
library(scales)
library(abind)
library(umap)
library(patchwork)
set.seed(42)
library(RColorBrewer)
library(furrr)

run_pseudotime <- function(feature_cube){
  # Retrieve wedge ids and join with age information
  df <- tibble(position_mat = 1:dim(feature_cube)[1],
                      id = dimnames(feature_cube)[[1]])

  # Calculate mean of distances across features
  mean_distances <- apply(feature_cube, c(1,2), mean)

  # Run diffusion map from distance matrix
  dimnames(mean_distances) <- NULL
  dm <- DiffusionMap(as.data.frame(df), distance = as.dist(mean_distances))
  # run pseudotime analysis
  dpt <- DPT(dm)

  # Add louvain clusters
  df <- df %>% mutate(louvain_transitions=get_louvain_clusters(dm@transitions))

  # Add pseudotime and pseudotime rank and DCs to dfs
  dpt_c <-dpt$dpt
  df <- df %>% mutate(dpt=dpt_c, dpt_rank=rank(dpt)) %>%
                left_join(as_tibble(dm@eigenvectors) %>%
                mutate(position_mat=df$position_mat), by='position_mat')

  return(df)
}

# Generate umap embedding from eigenvectors of the diffusion map and add to df_wedges
run_umap <- function(df, dims){
  df_tmp <- df %>% select(DC1:DC20)
  umap.embedding <- umap(df_tmp[,1:dims])
  df <- df %>% mutate(UMAP_1=umap.embedding$layout[,1],
                      UMAP_2=umap.embedding$layout[,2])
  return(df)
}

run_kmeans <- function(df, n=10){
  df$kmean <- as.factor(kmeans(df[,c("UMAP_1","UMAP_2")],n, nstart=100,iter.max=1000)[['cluster']])
  return(df)
}

run_louvain <- function(df, threshold=1.5){
  dist = df %>% select(UMAP_1, UMAP_2) %>% dist()
  dist[dist > threshold] = 0
  df <- df %>%
    mutate(louvain=dist %>% as.matrix() %>% get_louvain_clusters())
  df$louvain <- as.factor(df$louvain)
  return(df)
}

df_meta <- list.files('data/processed/4i/laminator/dataset_nuclei/', full.names = 1)[1] %>% read_csv()

run_spatial_heterogeneity_analysis <- function(max_radius, df_meta){
  # Load distance matrices
  message(paste('Started spatial heterogeneity analysis for maximal radius of',max_radius,'...'))
  message('Loading distance matrices...')
  file <- paste0('data/processed/4i/laminator/analysis_results/distance_matrices_fourier_nuclei_',max_radius,'.rds')
  M.list <- readRDS(file)
  # Scale distance matrices and bind to NxNxfeature array
  M.list <- map(M.list, function(x){rescale(log(x+1), to=c(0,1))})
  feature_cube <- abind(M.list, along=3)
  message('Running diffusion maps analysis...')
  df <- run_pseudotime(feature_cube) %>% run_umap(dims=5) %>% run_louvain()
  df <- df %>% left_join(df_meta %>% mutate_at(vars(id), as.character) %>% select(x,y,id,age,cell_type,organoid),by='id') %>%
    mutate(organoid = as.character(organoid))
  return(df)
}

radial_profile_sizes <- c(40,60,80,100,250,500)
names(radial_profile_sizes) <- radial_profile_sizes

message("Setting up the future...")
future::plan('multisession', workers = length(radial_profile_sizes))
results <- future_map(radial_profile_sizes, run_spatial_heterogeneity_analysis, df_meta=df_meta)

saveRDS(results, 'data/processed/4i/laminator/analysis_results/dpt_heterogeneity_w12_rgc.rds')
