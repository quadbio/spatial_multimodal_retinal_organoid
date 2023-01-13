source('laminator_analysis/utils.R')

# Load libraries
library(tidyverse)
library(destiny)
library(scales)
library(abind)
library(umap)

run_pseudotime <- function(timepoints=NULL, wedge_ids=NULL, df_meta, feature_cube){
  # Retrieve wedge ids and join with age information
  df_wedges <- tibble(position_mat = 1:dim(feature_cube)[1],
                      id = dimnames(feature_cube)[[1]]) %>%
    separate(id, into = c('organoid', 'window'), remove=FALSE) %>%
    left_join(df_meta %>% distinct(organoid, type), by='organoid')

  if(!is.null(timepoints)){
   df_wedges <- df_wedges %>% filter(type %in% timepoints)
  }
  if(!is.null(wedge_ids)){
    df_wedges <- df_wedges %>% filter(id %in% wedge_ids)
  }
  df_wedges %>% pull(id) -> filtered_profiles

  # Calculate mean of distances across features
  mean_distances <- apply(feature_cube[filtered_profiles,filtered_profiles,], c(1,2), mean)

  # Run diffusion map from distance matrix
  dimnames(mean_distances) <- NULL
  dm <- DiffusionMap(as.data.frame(df_wedges), distance = as.dist(mean_distances))
  # Run pseudotime analysis
  dpt <- DPT(dm)

  # Add louvain clusters
  df_wedges <- df_wedges %>% mutate(louvain_transitions=get_louvain_clusters(dm@transitions))

  # Add pseudotime and pseudotime rank and DCs to df_wedges
  dpt_c <-dpt$dpt
  df_wedges <- df_wedges %>% mutate(dpt=dpt_c, dpt_rank=rank(dpt)) %>%
                left_join(as_tibble(dm@eigenvectors) %>%
                mutate(position_mat=df_wedges$position_mat), by='position_mat')

  df_wedges$type <- factor(df_wedges$type, levels = c('week6', 'week12', 'week18', 'week24', 'week39','adult'))

  return(df_wedges)
}

run_umap <- function(df, dims){
  set.seed(42)
  df_tmp <- df %>% select(DC1:DC20)
  umap.embedding <- umap(df_tmp[,1:dims])
  df <- df %>% mutate(UMAP_1=umap.embedding$layout[,1],
                      UMAP_2=umap.embedding$layout[,2])
  return(df)
}

run_kmeans <- function(df, n=10){
  set.seed(42)
  df$kmean <- as.factor(kmeans(df[,c("UMAP_1","UMAP_2")],n, nstart=100,iter.max=1000)[['cluster']])
  return(df)
}

run_louvain <- function(df, threshold=1.5){
  set.seed(42)
  dist <- df %>% select(UMAP_1, UMAP_2) %>% dist()
  dist[dist > threshold] <- 0
  df <- df %>%
    mutate(louvain=dist %>% as.matrix() %>% get_louvain_clusters())
  df$louvain <- as.factor(df$louvain)
  return(df)
}

filter_wedges_by_bg_count <- function(df, max_frac_bg){
  df_bg <- read_csv('data/processed/4i/laminator/analysis_results/df_bg_counts.csv') %>%
  mutate(id=str_remove(wedge, '.png'),
         bg_frac = bg_count / (1000*100))

  df <- df %>% left_join(df_bg, by='id')

  df <- df %>% filter(bg_frac <= max_frac_bg)

  return(df)
}

# Load meta files
df_meta <- readRDS('data/processed/4i/laminator/dataset_protein/df_meta.rds')

# Load distance matrices
M_list <- readRDS('data/processed/4i/laminator/analysis_results/distance_matrices_fourier_protein.rds')

# Scale distance matrices and bind to NxNxfeature array
M_list <- map(M_list, function(x){rescale(log(x+1), to=c(0,1))})
feature_cube <- abind(M_list, along=3)
timepoints <- c('week6', 'week12', 'week18', 'week24', 'week39','adult')
names(timepoints) <- timepoints

pseudotime_res <- map(timepoints,run_pseudotime, df_meta=df_meta, feature_cube=feature_cube)
ids <- df_meta %>% filter_wedges_by_bg_count(max_frac_bg = 0.01) %>% distinct(id) %>% pull()
pseudotime_res_filtered <- map(timepoints,run_pseudotime, df_meta=df_meta, feature_cube=feature_cube, wedge_ids=ids)

# add UMAP embedding and cluster based on embedding
pseudotime_res <- map(pseudotime_res, run_umap, dims=10)
pseudotime_res <- map(pseudotime_res, run_kmeans, n=8)
pseudotime_res <- map(pseudotime_res, run_louvain)

# save results
saveRDS(pseudotime_res, 'data/processed/4i/laminator/analysis_results/diffusion_cluster.rds')

# rerun pseudotime with selected clusters
filter_selected_profile_clusters <- function(df_cluster) {

  exclude_list <- list(week39 = c(4, 5, 7, 11:14, 16:18, 20:22, 25:29),
                       week24 = c(1, 4, 6:9, 11:13, 15),
                       week18 = c(1, 2, 5, 7:9, 13:21, 23),
                       week12 = c(1, 2, 4, 8, 11, 12, 14, 19, 20),
                       week6 = c(2:4, 11:13, 16))
  excluded_clusters <- exclude_list %>%
    enframe() %>%
    unnest() %>%
    mutate(cluster_type = str_c(name, value, sep = '_')) %>%
    pull(cluster_type)

  df_cluster <- df_cluster %>%
    mutate(cluster_type = str_c(type, louvain, sep = '_')) %>%
    filter(!(cluster_type %in% excluded_clusters))

  return(df_cluster)
}

df_diffusion_res <-  pseudotime_res %>%
  bind_rows() %>% mutate(louvain=ifelse(type=='adult',1,louvain)) %>%
  mutate(cluster_type = str_c(type,louvain, sep='_')) %>%
  select(id, type, louvain, cluster_type, organoid) %>%
  filter_selected_profile_clusters()

selected_ids <- df_diffusion_res  %>% distinct(id) %>% pull()

df_pseudotime_refined <- run_pseudotime(timepoints=NULL,
                                        wedge_ids=selected_ids,
                                        df_meta=df_meta,
                                        feature_cube=feature_cube)

df_pseudotime_refined <- df_pseudotime_refined %>%
  run_umap(dims=10)

df_pseudotime_refined <- df_pseudotime_refined %>%
  left_join(df_diffusion_res %>%
              select(id, cluster_type), by='id')

# save results
saveRDS(df_pseudotime_refined,'data/processed/4i/laminator/analysis_results/laminar_window_dpt_trajectory.rds')
