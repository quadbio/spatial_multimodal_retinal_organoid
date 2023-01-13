# Load libraries
library(tidyverse)
library(destiny)
library(scales)
library(abind)
library(umap)

source('laminator_analysis/utils.R')
source('laminator_analysis/retina_fish_laminar_windows/utils_compare.R')

run_pseudotime <- function(wedge_ids=NULL,samples=NULL, df=df_meta, cube=feature_cube){
  # Retrieve wedge ids and join with age information
  df_wedges <- tibble(position_mat = 1:dim(cube)[1],
                      id = dimnames(cube)[[1]]) %>%
    separate(id, into = c('organoid', 'window'), remove=FALSE, sep='_') %>%
    left_join(df %>% distinct(organoid, type), by='organoid')

  if(!is.null(samples)){
   df_wedges <- df_wedges %>% filter(organoid %in% samples)
  }
  if(!is.null(wedge_ids)){
    df_wedges <- df_wedges %>% filter(id %in% wedge_ids)
  }
  df_wedges %>% pull(id) -> filtered_profiles

  # Calculate mean of distances across features
  mean_distances <- apply(cube[filtered_profiles,filtered_profiles,], c(1,2), mean)

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

  df_wedges$type <- factor(df_wedges$type, levels = c('week13', 'week32'))

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

run_louvain <- function(df, threshold=1.5){
  set.seed(42)
  dist <- df %>% select(UMAP_1, UMAP_2) %>% dist()
  dist[dist > threshold] <- 0
  df <- df %>%
    mutate(louvain=dist %>% as.matrix() %>% get_louvain_clusters())
  df$louvain <- as.factor(df$louvain)
  return(df)
}

# Load meta files
df_meta <- read_rds('data/processed/4i/laminator/dataset_protein/df_meta.rds') %>%
  mutate(type=ifelse(str_detect(organoid,'A|B'),'week32', 'week13'))

df_meta <- df_meta  %>% group_by(organoid) %>% filter(window %in% selected_windows[[unique(organoid)]]) %>% ungroup()
# Load distance matrices
M.list <- readRDS('data/processed/fish/laminator/analysis_results/distance_matrices_fourier_fish.rds')

# Scale distance matrices and bind to NxNxfeature array
M.list <- map(M.list %>% remove_incomplete_stains(), function(x){rescale(log(x+1), to=c(0,1))})
feature_cube <- abind(M.list, along=3)

ids <- df_meta %>% pull(id)

samples <- list('week13'=df_meta %>% filter(type=='week13') %>% pull(id),
                'week32'=df_meta %>% filter(type=='week32') %>% pull(id))


pseudotime_res <- map(samples,run_pseudotime)

# add UMAP embedding and cluster based on embedding
pseudotime_res <- map(pseudotime_res,run_umap,dims=10)
pseudotime_res <- map(pseudotime_res, run_louvain)
pseudotime_res <- bind_rows(pseudotime_res)

# save results
saveRDS(pseudotime_res, 'data/processed/fish/laminator/analysis_results/diffusion_cluster_fish.rds')
