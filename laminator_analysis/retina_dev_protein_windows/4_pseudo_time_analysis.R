# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 17.09.21

setwd('/home/charmel/4i_organoid_pipeline/retina/Rotator/Analysis')
source('../utils.R')

# Load libraries
library(tidyverse)
library(destiny)
library(scales)
library(abind)
library(umap)
library(ggridges)
library(egg)

set.seed(42)

# Load meta files
df_meta <- readRDS('/local2/USERS/charmel/rotator_dataset_2/df_meta.rds')

# Load distance matrices
M.list <- readRDS('/local2/USERS/charmel/rotator_analysis_results/distance_matrices_fourier_2.rds')

# Scale distance matrices and bind to NxNxfeature array
M.list <- map(M.list, function(x){rescale(log(x+1), to=c(0,1))})
feature_cube <- abind(M.list, along=3)

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

timepoints <- c('week6', 'week12', 'week18', 'week24', 'week39','adult')
names(timepoints) <- timepoints

pseudotime_res <- map(timepoints,run_pseudotime, df_meta=df_meta, feature_cube=feature_cube)
ids <- df_meta %>% filter_wedges_by_bg_count(max_frac_bg = 0.01) %>% distinct(id) %>% pull()
pseudotime_res_filtered <- map(timepoints,run_pseudotime, df_meta=df_meta, feature_cube=feature_cube, wedge_ids=ids)

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
  dist <- df %>% select(UMAP_1, UMAP_2) %>% dist()
  dist[dist > threshold] <- 0
  df <- df %>%
    mutate(louvain=dist %>% as.matrix() %>% get_louvain_clusters())
  df$louvain <- as.factor(df$louvain)
  return(df)
}

# add UMAP embedding and cluster based on embedding
pseudotime_res <- map(pseudotime_res, run_umap, dims=10)
pseudotime_res <- map(pseudotime_res, run_kmeans, n=8)
pseudotime_res <- map(pseudotime_res, run_louvain)

# plotting
bind_rows(pseudotime_res) %>% ggplot(aes(x=UMAP_1, y=UMAP_2, col=louvain)) +
  geom_point(size=1)+
  coord_equal() +
  theme_void() +
  geom_text(data = bind_rows(pseudotime_res) %>%
    group_by(louvain, type) %>%
                summarise(x=mean(UMAP_1), y=mean(UMAP_2)),
              aes(x=x, y=y, label=louvain), col='black') +
  facet_wrap(~type)

ggplot(data=NULL) +
  geom_point(data=pseudotime_res[['week39']],aes(x=UMAP_1, y=UMAP_2, col=louvain))+
  coord_equal() +
  theme_void()+
  geom_text(data = pseudotime_res[['week39']] %>%
    group_by(louvain) %>%
                summarise(x=mean(UMAP_1), y=mean(UMAP_2)),
              aes(x=x, y=y, label=louvain))

# save results
saveRDS(pseudotime_res, '/local2/USERS/charmel/rotator_analysis_results/diffusion_cluster.rds')

# rerun pseudotime with selected clusters
diffusion_res_list <- readRDS('/local2/USERS/charmel/rotator_analysis_results/diffusion_cluster_filtered.rds')

df_diffusion_res <-  diffusion_res_list %>%
  bind_rows() %>% mutate(louvain=ifelse(type=='adult',1,louvain)) %>%
  mutate(cluster_type = str_c(type,louvain, sep='_')) %>% select(id, type, louvain, cluster_type, organoid) %>%
  filter_selected_profile_clusters(type='manual_pw_include_w24')

selected_ids <- df_diffusion_res  %>% distinct(id) %>% pull()

df_pseudotime_refined <- run_pseudotime(timepoints=NULL, wedge_ids=selected_ids, df_meta=df_meta, feature_cube=feature_cube)

df_pseudotime_refined <- df_pseudotime_refined %>% run_umap(dims=10)

df_pseudotime_refined <- df_pseudotime_refined %>% left_join(df_diffusion_res %>% select(id, cluster_type), by='id')

# save results
saveRDS(df_pseudotime_refined,'/local2/USERS/charmel/pseudotime_trajectory_final/pseudotime_trajectory.rds')

# plotting
p1 <- ggplot(df_pseudotime_refined, aes(x=UMAP_1, y=UMAP_2, col=type)) +
  geom_point(size=.5)+
  coord_equal() +
  theme_article() +
  guides(col=guide_legend(title="Age"))+
  scale_color_viridis_d()
p2 <- ggplot(df_pseudotime_refined, aes(x=UMAP_1, y=UMAP_2, col=dpt_rank)) +
  geom_point(size=.5) +
  coord_equal() +
  theme_article() +
  labs(color="Pseudotime rank")+
  scale_color_viridis_c(direction = -1, option="magma")
p3 <- ggplot(df_pseudotime_refined, aes(x=DC1, y=DC2, col=type)) +
  geom_point() +
  theme_article() +
  guides(col=guide_legend(title="Age"))+
  scale_color_viridis_d()

p4 <- ggplot(df_pseudotime_refined, aes(x=max(dpt_rank)-dpt_rank, y=type, fill=type)) +
  geom_density_ridges(alpha=0.75) +
  theme_article() +
  guides(fill=guide_legend(title="Age"))+
  scale_fill_viridis_d()

ggarrange(plots = list(p1,p3,p2,p4),
          heights = c(1,1),
          widths = c(1,1),
          nrow = 2, ncol = 2,
          top = textGrob(label = 'Pseudotime selected clusters all timepoints',
                         gp = gpar(fontsize = 12)))

df_pseudotime_refined$cluster_type <- factor(df_pseudotime_refined$cluster_type,
                                             levels = df_pseudotime_refined %>%
                                               group_by(cluster_type) %>%
                                               summarise(med_rank = median(dpt_rank)) %>%
                                               arrange(-med_rank) %>%
                                               pull(cluster_type))

df_pseudotime_refined %>%
  ggplot(aes(x=max(dpt_rank)-dpt_rank, y=cluster_type ,fill=cluster_type)) +
  geom_density_ridges(alpha=0.75) + theme_article()

map_contour_to_image(df_pseudotime_refined %>% left_join(df_meta %>% select(id, x, y), by='id'),
                     point = "43", annotation='dpt_rank', return_gg = TRUE)

# rerun with updated clustering for week 39
distatis_res <- readRDS('/local2/USERS/charmel/rotator_analysis_results/distatis_from_fourier_2.rds')

diffusion_res_list_v2 <- readRDS('/local2/USERS/charmel/rotator_analysis_results/diffusion_cluster_filtered.rds')
diffusion_res_list_v2$week39 <- run_louvain(diffusion_res_list$week39, threshold=0.75)

map_contour_to_image(diffusion_res_list$week39 %>% left_join(df_meta %>% select(id, x, y), by='id'),
                     point = "43", annotation='louvain', return_gg = TRUE)

df_diffusion_res <-  diffusion_res_list_v2 %>%
  bind_rows() %>% mutate(louvain=ifelse(type=='adult',1,louvain)) %>%
  mutate(cluster_type = str_c(type,louvain, sep='_')) %>% select(id, type, louvain, cluster_type, organoid) %>%
  filter_selected_profile_clusters(type='manual_pw_include_w24')






