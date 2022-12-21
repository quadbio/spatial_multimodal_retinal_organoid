# Title     : 5_graph_building.R
# Objective : Building a graph network of timewise wedge clusters gained from diffusion maps
# Created by: harmelc
# Created on: 30.09.21


# load libraries
library(tidyverse)
library(tidygraph)
library(ggraph)
library(abind)
library(DistatisR)
library(furrr)
library(scales)

# read distance matrices and convert to array
M.dist <- readRDS('/local2/USERS/charmel/rotator_analysis_results/distance_matrices_fourier_2.rds')
feature_cube <- abind(M.dist, along=3)
rm(M.dist)

# read in diffusion map results
diffusion_res <- readRDS( '/local2/USERS/charmel/rotator_analysis_results/diffusion_cluster_filtered.rds') %>%
  bind_rows() %>% mutate(louvain=ifelse(type=='adult',1,louvain)) %>%
  mutate(cluster_type = str_c(type,louvain, sep='_'))


filter_selected_profile_clusters <- function(df_cluster){
  exclude.list <- list(week39=c(5,6,9,11,12,13),
                       week24=NULL,
                       week18=c(8,9),
                       week12=c(4,5,8,9),
                       week6=1)
  exclude.list %>% enframe() %>%
    unnest() %>%
    mutate(cluster_type=str_c(name, value, sep='_')) %>%
    pull(cluster_type) -> excluded_clusters

  df_cluster %>% mutate(cluster_type=str_c(type,louvain,sep='_')) %>%
    filter(!(cluster_type %in% excluded_clusters))
}

#diffusion_res <- diffusion_res %>% filter_selected_profile_clusters()


# define nodes
df_nodes <- diffusion_res %>%
  group_by(type) %>% mutate(n_type=n()) %>% ungroup() %>%
  group_by(cluster_type) %>% mutate(n_cluster_type=n()) %>%
  distinct(cluster_type, n_type, n_cluster_type) %>%
  summarise(size_norm=n_cluster_type/n_type)

# define edges
df_edges <- diffusion_res %>% distinct(louvain, type) %>% unite(cluster_type) %>%
  expand(cluster_type, cluster_type_to=cluster_type) %>%
  rename(from = cluster_type,
         to= cluster_type_to)

# add edge ids
df_edges$edge_id <- df_edges %>% group_by(from, to) %>% group_indices()

# store edge ids
edge_ids <- df_edges %>% group_by(edge_id) %>% nest() %>% deframe()


# define weighting function
get_edge_weight <- function(edge, df, distances, distatis_input=FALSE){
  print(edge)
  from_ids <- df %>% filter(cluster_type == edge$from) %>% pull(id)
  to_ids <- df %>% filter(cluster_type == edge$to) %>% pull(id)
  all_ids <- c(from_ids, to_ids)
  if(!distatis_input){
    distances <- distances[all_ids,all_ids,]
    distatis_res <- distatis(distances)
    dist_comp <- distatis_res$res4Splus$Splus
  } else {
    dist_comp <- distances[all_ids, all_ids]
  }
  diag(dist_comp) <- NaN
  dist_comp <- dist_comp[from_ids,to_ids]
  mean_dist <- mean(dist_comp, na.rm = 1)
  sd_dist <- sd(dist_comp, na.rm = 1)
  return(list(mean=mean_dist,sd=sd_dist))
}

future::plan('multisession', workers = 20)
options(future.globals.maxSize= 891289600*100)
weights <- future_map(edge_ids, get_edge_weight, df=diffusion_res, distances=feature_cube)
weights <- future_map(edge_ids, get_edge_weight, df=diffusion_res, distances=distatis_res$res4Splus$Splus, distatis_input=TRUE)

df_weights <- weights %>% bind_rows(.id='edge_id') %>% mutate(edge_id = as.numeric(edge_id))

df_edges_combined <- full_join(df_edges, df_weights) %>%
 separate(from, remove = FALSE, into = c('type_from','cluster_from')) %>%
  separate(to, remove = FALSE, into = c('type_to','cluster_to'))


filter_edges <- function(df_edges, k=7, timepoints = c('week6','week12','week18','week24','week39','adult')){
  df_edges <- df_edges %>% mutate(weights = rescale(mean))
  df_filtered <- tibble()
  for ( i in seq_along(timepoints[-6])){
    df_tmp <- df_edges %>%
      filter(type_from == timepoints[i],
             !(from==to))
    if(i>1){
      df_tmp <- df_tmp %>% filter(type_to %in% timepoints[c(i-1,i,i+1)]) %>%
        group_by(from) %>% arrange(from, -mean) %>% slice(1:k)
    } else {
      df_tmp <- df_tmp %>% filter(type_to %in% timepoints[c(i,i+1)]) %>%
        group_by(from) %>% arrange(from, -mean) %>% slice(1:k)
    }
    df_filtered <- bind_rows(df_filtered, df_tmp)
  }
  df_tmp <- df_edges %>% filter(type_from=="adult", !(from==to), type_to=='week39') %>%
    group_by(from) %>% arrange(from, -mean) %>% slice(1:k)

  return(bind_rows(df_filtered, df_tmp))
}


random_walk_frequencies <- function(df_edges, start_clusters='week6', n=10) {
  df_results <- tibble()
  df_edges <- df_edges %>% mutate(weights=rescale(mean))
  weights_timepoint <- c('week6' =1,
                        'week12'=2,
                        'week18'=3,
                        'week24'=4,
                        'week39'=5,
                        'adult'= 6) %>% rescale(c(0.0001,1))
  for(i in 1:n){
    print(i)
    df_counts <- df_edges %>% distinct(from) %>% mutate(count=0)

    position <- df_counts %>% filter(str_detect(from,start_clusters)) %>%
      sample_n(1) %>% pull(from)

    while (position!='adult_1') {
      df_counts <- df_counts %>% mutate(count = ifelse(from == position, count+1, count))
      df_tmp1 <- df_edges %>% filter(from==position)
      if(nrow(df_tmp1)>0){
        position_1 <- df_tmp1 %>%
          mutate(weights_time = weights_timepoint[type_to]) %>%
          sample_n(1, weight = weights + weights_time) %>% pull(to)
      } else {
        position_1 <- NA
      }
      df_tmp2 <- df_edges %>% filter(to==position)
      if(nrow(df_tmp2)>0){
        position_2 <- df_tmp2 %>%
          mutate(weights_time = weights_timepoint[type_from]) %>%
          sample_n(1,  weight = weights + weights_time) %>% pull(from)
      } else {
        position_2 <- NA
      }
      position <- sample(na.omit(c(position_1, position_2)),1)
    }
    df_counts <- df_counts %>% mutate(count = ifelse(from == position, count+1, count))
    df_counts <- df_counts %>% mutate(count=ifelse(count>0,1,0))
    df_results <- bind_rows(df_results, df_counts %>% mutate(freq_iter=count, iter=i))
  }

  df_results <- df_results %>% group_by(from) %>%
    summarise(freq=mean(freq_iter)) %>% rename(cluster_type=from)

  return(df_results)
}


df_edges_2 <- df_edges_combined  %>%
  filter_edges(k=5)
df_nodes_2 <- random_walk_frequencies(df_edges_2, n=10)
df_nodes_2 <- left_join(df_nodes, df_nodes_2, by='cluster_type')
graph <- tbl_graph(df_nodes_2 %>% separate(cluster_type,
                                           remove = FALSE,
                                           into = c('type','cluster')),
                                           df_edges_2)

p1 <- ggraph(graph , layout="stress", weights=weights) +
  geom_edge_link(aes(alpha=weights)) + geom_node_point(aes(col=freq, size=size_norm)) +
  geom_node_text(aes(label=cluster), col="white") +
  theme_void() + scale_size(range=c(5,15)) +
  scale_alpha_continuous(range = c(0,0.5))+
  scale_color_viridis_c(option='magma')+
  guides(size='none',
         alpha= 'none')

p2 <- ggraph(graph %>% activate(edges) , layout = "stress", weights = weights) +
  geom_edge_link(aes(alpha = rescale(weights))) + geom_node_point(aes(col = type, size = size_norm)) +
  geom_node_text(aes(label = cluster), col = "white") +
  theme_void() + scale_size(range = c(5,15)) +
  scale_alpha_continuous(range = c(0,0.5))+
  guides(size = 'none',
         alpha = 'none')

fig <- image_graph()
p2
dev.off()
image_browse(fig)
graph <- tbl_graph(df_nodes %>%
                     separate(cluster_type, remove = FALSE, into = c('type','cluster')),
                   df_edges_combined %>% group_by(from) %>% arrange(from, -mean) %>% slice(1:5), directed = TRUE)

graph <- tbl_graph(df_nodes %>%
                     separate(cluster_type, remove = FALSE, into = c('type','cluster')) %>%
                     filter(type!='week24'),
                   df_edges_combined %>%
                     filter_edges(k=5, timepoints = c('week6','week12','week18','week39','adult')),
                   directed = TRUE)



graph <- tbl_graph(df_nodes %>%
                     separate(cluster_type, remove = FALSE, into = c('type','cluster')),
                   df_edges_combined %>%
                     group_by(from) %>% arrange(from, -mean) %>%
                     slice(2:4), directed = TRUE) %>% activate(edges) %>%
          mutate(weights = rescale(mean))


p2 <- ggraph(graph, layout="stress", weights= weights) +
  geom_edge_link(alpha=0.25) + geom_node_point(aes(col=type, size=size_norm)) +
  geom_node_text(aes(label=cluster), col="white") +
  theme_void() + scale_size(range=c(5,15)) +
  scale_alpha_continuous(range = c(0,0.5))+
  guides(size='none',
         alpha= 'none')


get_edge_weight_check <- function(edge, df, distances){
  from_ids <- df %>% filter(cluster_type == edge$from) %>% pull(id)
  to_ids <- df %>% filter(cluster_type == edge$to) %>% pull(id)
  all_ids <- c(from_ids, to_ids)
  distances <- distances[all_ids,all_ids,]
  dist_comp_2 <- apply(distances, c(1,2), mean)
  diag(dist_comp_2) <- NaN
  distatis_res <- distatis(distances)
  dist_comp <- distatis_res$res4Splus$Splus
  diag(dist_comp) <- NaN
  list('cross_av'= dist_comp_2[from_ids,to_ids] %>% as.vector(),
       'from_av'=dist_comp_2[from_ids,from_ids] %>% as.vector(),
       'to_av'=dist_comp_2[to_ids,to_ids] %>% as.vector(),
       'cross_c'= dist_comp[from_ids,to_ids] %>% as.vector(),
       'from_c'=dist_comp[from_ids,from_ids] %>% as.vector(),
       'to_c'=dist_comp[to_ids,to_ids] %>% as.vector())
}


get_edge_weight_check(edge_ids[['1613']], diffusion_res, feature_cube) %>%
  enframe() %>% unnest() %>%
  separate(name, into = c('type','method'), sep ='_') %>%
  ggplot(aes(x=value, fill=type)) +
  geom_density(alpha=0.25) +
  facet_wrap(~method, scales="free")