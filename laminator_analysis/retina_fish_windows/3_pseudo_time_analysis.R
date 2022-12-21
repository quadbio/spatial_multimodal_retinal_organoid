# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 17.09.21

setwd('/home/charmel/4i_organoid_pipeline/retina/rotator_analysis/retina_fish_windows')
source('../utils.R')
source('0_utils_compare.R')
# Load libraries
library(tidyverse)
library(destiny)
library(scales)
library(abind)
library(umap)
library(ggridges)
library(egg)
library(magick)
library(RColorBrewer)

remove_incomplete_stains <- function(dist.list){
  dim <- map(dist.list, nrow)
  dim <- dim[!map(dim, is.null) %>% unlist()]
  dims_unique <- unlist(dim) %>% unique() %>% sort(decreasing = TRUE)
  dim_keep <- dims_unique[1]
  dim <- dim[map(dim, function(x){x==dim_keep}) %>% unlist()]
  return(dist.list[names(dim)])
}
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
  df_tmp <- df %>% select(DC1:DC20)
  umap.embedding <- umap(df_tmp[,1:dims])
  df <- df %>% mutate(UMAP_1=umap.embedding$layout[,1],
                      UMAP_2=umap.embedding$layout[,2])
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

set.seed(42)

# Load meta files
df_meta <- read_rds('/local2/USERS/charmel/rotator_dataset_fish_2/df_meta.rds') %>%
  mutate(type=ifelse(str_detect(organoid,'A|B'),'week32', 'week13'))

selected_windows <- list('A1-1' = 0:136,
                         'A1-2' = c(0:2, 28:132, 300:334),
                         'A1-3' = c(0:30, 134:244),
                         'A2-1' = 0:149,
                         'A2-2' = c(0:154, 188:239),
                         'A2-3' = c(0:20, 39:82, 104:161),
                         'B1-1' = 0:166,
                         'B1-2' = c(0:4,27:131,375:458),
                         'B1-3' = c(0:8,30:116,246:252),
                         'B2-1' = c(0:8,32:54,68:168),
                         'B2-2' = c(0:11,29:35,99:214),
                         'B2-3' = 328:420,
                         'D1-2' = c(0:9,38:67,102:176,194:290),
                         'D1-3' = c(0:26,49:125,166:184,206:321),
                         'C1-1' = c(0:1,22:105,131:276),
                         'C1-2' = c(1:45,144:233),
                         'C1-3' = c(0:8,27:33,72:112, 625:702,717:760,788:800,827:900),
                         'C2-1' = c(65:116,132:266,290:330,540:638),
                         'C2-2' = c(0:4,35:235),
                         'C2-3' = c(0:304,323:361)
                         )

# df_meta  %>% group_by(organoid) %>% filter(window %in% selected_windows[[unique(organoid)]]) %>% ungroup() %>% ggplot(aes(x=x, y=y)) + geom_point() + facet_wrap(~organoid)
#
# plot_contour_positions_on_organoid('B2-3')
# df_meta_select <- df_meta %>% filter(organoid == 'C2-1')
# p <- ggplot(df_meta_select %>% mutate(Y =  - y, X = x - 1), aes(X,Y, col=window)) + geom_point() + scale_color_viridis_c()
# ggplotly(p)



df_meta <- df_meta  %>% group_by(organoid) %>% filter(window %in% selected_windows[[unique(organoid)]]) %>% ungroup()
# Load distance matrices
M.list <- readRDS('/local2/USERS/charmel/rotator_analysis_results_fish/distance_matrices_2.rds')

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

# image mapping function
map_to_image <- function(df, image) {

  image <- image_read(image)  %>% image_normalize()
  width <- image_info(image) %>% pull(width)
  height <- image_info(image) %>% pull(height)

  df_plot <- df  %>% mutate(Y = height - y, X = x - 1)
  if (unique(df$type) == 'week13'){
    image <- image_extent(image, '17152x15008', gravity = 'southwest', color='black')
  } else {
    image <- image_extent(image, '8576x8576', gravity = 'southwest', color='black')
  }

  width <- image_info(image) %>% pull(width)
  height <- image_info(image) %>% pull(height)
  # set image as background in the plot
  bg <- rasterGrob(image, width = unit(1, 'npc'), height = unit(1, 'npc'), interpolate = TRUE)

  p1 <- ggplot(df_plot, aes(x = X, y = Y, col=louvain)) +
    coord_fixed() +
    annotation_custom(grob = bg,
                      xmin = 0,
                      xmax = width,
                      ymin = 0,
                      ymax = height) +
    scale_x_continuous(limits = c(0, width),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, height),
                       expand = c(0, 0)) +
    geom_point(size=0.5)+
    scale_color_manual(values = cluster_cols) +
    xlab('') + ylab('')  #transparent legend panel
  return(p1)
}

dir_dapi <- '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission/DAPI-lower exposure_submission'

image_files <- list.files(dir_dapi)
names(image_files) <- image_files %>% str_remove('32955-slide928-2_') %>% str_remove('_DAPI.tiff')
cluster_cols <- brewer.pal(pseudotime_res %>% distinct(louvain) %>% nrow(),'Set3')
names(cluster_cols) <- pseudotime_res %>% distinct(louvain) %>% pull()

plot_contour_cluster_on_organoid <- function (sample){
  p <- pseudotime_res %>%
    left_join(df_meta %>% select(x, y, id), by = 'id') %>%
    filter(organoid == sample) %>%
    map_to_image(image = paste0(dir_dapi, '/', image_files[sample]))
  p <- p + theme_void()
  return(p)
}
samples <- pseudotime_res %>% distinct(organoid) %>% pull() %>% set_names()
plots_samples <- map(samples, plot_contour_cluster_on_organoid)

# plotting
p1 <- pseudotime_res %>% ggplot(aes(x=UMAP_1, y=UMAP_2, col=louvain)) +
  geom_point(size=1)+
  coord_equal() +
  theme_void() +
  geom_text(data = pseudotime_res %>%
    group_by(louvain, type) %>%
                summarise(x=mean(UMAP_1), y=mean(UMAP_2)),
              aes(x=x, y=y, label=louvain), col='black') +
  scale_color_manual(values = cluster_cols) +
  facet_wrap(~type)

p2 <- pseudotime_res %>% ggplot(aes(x=UMAP_1, y=UMAP_2, col=organoid)) +
  geom_point(size=1)+
  coord_equal() +
  theme_void() +
  facet_wrap(~type)

library(patchwork)
p10 <- (plots_samples[['A1-1']] | plots_samples[['A1-2']] | plots_samples[['A1-3']])/(plots_samples[['A2-1']] | plots_samples[['A2-2']] | plots_samples[['A2-3']])/(plots_samples[['B1-1']] | plots_samples[['B1-2']] | plots_samples[['B1-3']])/(plots_samples[['B2-1']] | plots_samples[['B2-2']] | plots_samples[['B2-3']]) + plot_layout(guides='collect')
p11 <- (plots_samples[['C1-1']] | plots_samples[['C1-2']] | plots_samples[['C1-3']] | plots_samples[['C2-1']]) / (plots_samples[['C2-2']] | plots_samples[['C2-3']] | plots_samples[['D1-2']] | plots_samples[['D1-3']]) + plot_layout(guides='collect')

ggsave(plot=p10,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/tissue_heterogeneity_fish_2_week32.png'),
       dpi = 250, unit='mm', width =400.05 , height = 215.22)
ggsave(plot=p11,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/tissue_heterogeneity_fish_2_week13.png'),
       dpi = 250, unit='mm', width =400.05 , height = 215.22)
p12 <-  p1/p2
ggsave(plot=p12,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/tissue_heterogeneity_fish_embeddings_2.png'),
       dpi = 250, unit='mm', width =400.05 , height = 215.22)
# save results
saveRDS(pseudotime_res, '/local2/USERS/charmel/rotator_analysis_results_fish/diffusion_cluster_final_2.rds')

dir_dapi <- '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission/aligned'
image_files <- list.files(dir_dapi, full.names=1)
names(image_files) <- list.files(dir_dapi) %>% str_remove('.tif')

#image_files <- list.files(dir_dapi)
#names(image_files) <- image_files %>% str_remove('32955-slide928-2_') %>% str_remove('_DAPI.tiff')
cluster_cols <- brewer.pal(pseudotime_res %>% distinct(louvain) %>% nrow(),'Set3')
names(cluster_cols) <- pseudotime_res %>% distinct(louvain) %>% pull()


pseudotime_res %>% filter(organoid =='C1-2') %>% left_join(df_meta %>% select(x, y, id), by = 'id') %>% na.omit() %>% map_to_image(image = image_files['C1-2'])


# image mapping function
map_to_image <- function(df, image) {

  image <- image_read(image)  %>% image_normalize()
  width <- image_info(image) %>% pull(width)
  height <- image_info(image) %>% pull(height)

  df_plot <- df  %>% mutate(Y = height - y, X = x - 1)

  width <- image_info(image) %>% pull(width)
  height <- image_info(image) %>% pull(height)
  # set image as background in the plot
  bg <- rasterGrob(image, width = unit(1, 'npc'), height = unit(1, 'npc'), interpolate = TRUE)

  p1 <- ggplot(df_plot, aes(x = X, y = Y, col=louvain)) +
    coord_fixed() +
    annotation_custom(grob = bg,
                      xmin = 0,
                      xmax = width,
                      ymin = 0,
                      ymax = height) +
    scale_x_continuous(limits = c(0, width),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, height),
                       expand = c(0, 0)) +
    geom_point(size=0.5)+
    scale_color_manual(values = cluster_cols) +
    xlab('') + ylab('')  #transparent legend panel
  return(p1)
}

plot_contour_cluster_on_organoid <- function (sample){
  p <- pseudotime_res %>%
    left_join(df_meta %>% select(x, y, id), by = 'id') %>%
    filter(organoid == sample) %>%
    map_to_image(image = image_files[sample])
  p <- p + theme_void() + geom_rect(data = NULL, xmin = 100, xmax = 825, ymin = 100, ymax = 200, col='white', fill='white')
  return(p)
}
samples <- pseudotime_res %>% distinct(organoid) %>% pull() %>% set_names()
plots_samples <- map(c('B1-3','C1-2') %>% set_names(), plot_contour_cluster_on_organoid)


p <- plots_samples[['B1-3']]  | plots_samples[['C1-2']]
ggsave(plot=p,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/window_cluster_overlay_examples.png'),
       dpi = 250, unit='mm', width =400.05/2 , height = 215.22/2)
