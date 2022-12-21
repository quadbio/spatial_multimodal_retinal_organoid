# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 26.01.22

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


df_meta <- list.files('/local2/USERS/charmel/rotator_dataset_nuclei', full.names = 1)[1] %>% read_csv()


run_spatial_heterogeneity_analysis <- function(max_radius, df_meta){
  # Load distance matrices
  message(paste('Started spatial heterogeneity analysis for maximal radius of',max_radius,'...'))
  message('Loading distance matrices...')
  file <- paste0('/local2/USERS/charmel/rotator_analysis_results/distance_matrices_fourier_nuclei_',max_radius,'.rds')
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

results <- readRDS('/local2/USERS/charmel/dpt_heterogeneity_w12_rgc.rds') %>% bind_rows(.id='max_radius') %>% distinct()



organoids <- df %>% distinct(organoid) %>% pull()
organoid_cols <- brewer.pal(length(organoids),'Set3')
names(organoid_cols) <- organoids
clusters <- df %>% distinct(louvain) %>% pull()
cluster_cols <- brewer.pal(length(clusters),'Set3')
names(cluster_cols) <- clusters


# plotting
p1 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, col=louvain)) +
  geom_point(size=1)+
  coord_equal() +
  theme_void() +
  geom_text(data = df %>%
    group_by(louvain, max_radius) %>%
                summarise(x=mean(UMAP_1), y=mean(UMAP_2)),
              aes(x=x, y=y, label=louvain), col='black') +
  scale_color_manual(values=cluster_cols) +
  facet_wrap(~max_radius)
p2 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, col=organoid)) +
  geom_point(size=1)+
  coord_equal() +
  theme_void() +
  scale_color_manual(values=organoid_cols)+
  facet_wrap(~max_radius)



ggsave(plot=p1 | p2,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots_v2/tissue_heterogeneity_rgc_umaps.png'),
       dpi = 250, unit='mm', width = 400.05, height = 215.22)


plot_feature_on_organoids <- function(max_radius, df, organoids, cluster_cols){
  df <- df %>% filter(max_radius == as.character(max_radius))
  for (i in seq_along(organoids)){
    p <- map_contour_to_image(df, point = organoids[i], annotation='louvain', return_gg = TRUE) +
      scale_color_manual(values=cluster_cols)
    ggsave(plot=p,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots_v2/rgc_heterogeneity/',max_radius,'_',organoids[i],'.png'),
       dpi = 250, unit='mm', width = 400.05, height = 215.22)
  }
}

map(radial_profile_sizes, plot_feature_on_organoids, df=df, organoids=organoids, cluster_cols=cluster_cols)

library(magick)



crop_circle <- function(im,size){
 # create a new image with white background and black circle
    radius = size * 100
    fig <- image_draw(image_blank(radius, radius))
    symbols(radius/2, radius/2, circles=(radius/2)-3, bg='black', inches=FALSE, add=TRUE)
    dev.off()
    fig <- image_scale(fig, str_c(size,size,sep = 'x'))
    # create an image composite using both images
    im2 <- image_composite(im, fig, operator='copyopacity')

  # set background as white
  return(image_background(im2, 'black'))
}



crop_string_df <- function(x_min, x_max, y_min, y_max){
  window_x <- x_max - x_min
  window_y <- y_max - y_min
  center_x <- x_min + window_x/2
  center_y <- y_min + window_y/2
  offset_x <- center_x - (window_x / 2)
  offset_x <- paste("+", offset_x, sep = "")
  offset_y <- center_y - (window_y / 2)
  offset_y <- paste("+", offset_y, sep = "")
  crop_string <- paste0(paste(window_x), "x", paste(window_y), paste(offset_x), paste(offset_y))
  return(crop_string)
}

create_neighborhood_collage <- function(radius=40, df, n=12, img_type='nuclei'){
  dir_list <- list(nuclei='/links/groups/treutlein/DATA/imaging/charmel/shiny_input/',
                   mtu='/nas/groups/treutlein/DATA/imaging/charmel/wedges_from_clustering/composite_imgs/')
  dir <- dir_list[[img_type]]
  df_positions <- df %>% filter(max_radius == as.character(radius)) %>% select(max_radius, x, y, organoid, louvain, id)  %>%
    mutate(x_min = x - radius,
           x_max = x + radius,
           y_max = y + radius,
           y_min = y - radius,
           crop_string = crop_string_df(x_min, x_max, y_min, y_max))
  organoids <- df_positions %>% distinct(organoid) %>% pull()
  names(organoids) <- organoids
  file= list(nuclei='/channel_hoechst.tif',
               mtu='.png')
  df_image_info <- map(organoids, function(x){str_c(dir,x,file[[img_type]]) %>% image_read() %>% image_info()}) %>% bind_rows(.id='organoid')
  df_positions <- df_positions %>% left_join(df_image_info %>% select(organoid, width, height), by='organoid') %>% group_by(id) %>%
    mutate(bound_x_min = x_min<0, bound_x_max=x_max > width, bound_y_min = y_min < 0, bound_y_max = y_max > height,
           out_of_bounds = any(bound_x_min, bound_x_max, bound_y_min, bound_y_max)) %>% filter(out_of_bounds == 0) %>%
    group_by(louvain) %>% mutate_at(vars(louvain),as.character) %>%
    sample_n(n) %>% ungroup()
  clusters <- df_positions %>% distinct(louvain) %>% pull()
  clusters_list <- vector('list', length=length(clusters))
  for (i in seq_along(clusters)){
    organoids_list <- vector('list', length = length(organoids))
    for (j in seq_along(organoids)){
      df_tmp <- df_positions %>% filter(organoid ==organoids[j], louvain == clusters[i])
      if(nrow(df_tmp)==0){
        next
      }
      image <- str_c(dir,organoids[j],file[[img_type]]) %>% image_read() %>% image_normalize()
      crop_strings <- df_tmp %>% pull(crop_string)
      images <- map(crop_strings, function(x){image_crop(image, geometry = x)})
      organoids_list[[j]] <- images
    }
    clusters_list[[i]] <- organoids_list
  }
  clusters_list <- map(clusters_list, image_join)
  names(clusters_list) <- clusters
  clusters_list <- map(clusters_list, crop_circle, size=2*radius)
  clusters_list <- map(clusters_list,image_montage, tile='2x6', bg='black') %>% map(image_border, geometry='2x2')
  map(clusters, function(x){clusters_list[[x]] %>% image_annotate(as.character(x), size = 25, color = "white", gravity = 'center')}) %>%
    image_join() %>% image_append() %>%
    image_write(paste0('/links/groups/treutlein/DATA/imaging/charmel/plots_v2/collage_',img_type,'_',radius,'.png'), format = 'png')
}

create_neighborhood_collage(radius=40, df, n=12)
create_neighborhood_collage(radius=60, df, n=12)
create_neighborhood_collage(radius=80, df, n=12)
create_neighborhood_collage(radius=100, df, n=12)


create_neighborhood_collage(radius=40, df, n=12, img_type='mtu')
create_neighborhood_collage(radius=60, df, n=12, img_type='mtu')
create_neighborhood_collage(radius=80, df, n=12, img_type='mtu')
create_neighborhood_collage(radius=100, df, n=12, img_type='mtu')