library(TSdist)
library(tidyverse)
library(Seurat)
library(furrr)
library(FBN)
library(pracma)
library(secr)
library(abind)
library(scales)
library(furrr)
library(magick)
library(grid)
library(gridExtra)

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
                         'C1-2' = c(0:119,144:233),
                         'C1-3' = c(0:8,27:33,72:112, 625:702,717:760,788:800,827:900),
                         'C2-1' = c(65:116,132:266,290:330,540:638),
                         'C2-2' = c(0:4,35:235),
                         'C2-3' = c(0:304,323:361)
                         )

plot_contour_positions_on_organoid <- function (sample, range=NULL){
  if(!is.null(range)){
    p <- df_meta  %>% filter(organoid==sample) %>% filter(window %in% range) %>%
         map_to_image(image = paste0(dir_dapi,'/',image_files[sample]))
  } else {
    p <- df_meta  %>% filter(organoid==sample) %>%
         map_to_image(image = paste0(dir_dapi,'/',image_files[sample]))
  }

  p <- p + theme_void()
  return(p)
}
# image mapping function
map_to_image <- function(df, image) {

  image <- image_read(image)  %>% image_normalize()
  width <- image_info(image) %>% pull(width)
  height <- image_info(image) %>% pull(height)

  df_plot <- df  %>% mutate(Y = height - y, X = x - 1)
  # set image as background in the plot
  bg <- rasterGrob(image, width = unit(1, 'npc'), height = unit(1, 'npc'), interpolate = TRUE)

  p1 <- ggplot(df_plot, aes(x = X, y = Y, col=window)) +
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
    geom_point()+
    scale_color_viridis_c() +
    xlab('') + ylab('')  #transparent legend panel
  return(p1)
}

extract_celltypes_by_dpt_rank <- function(result.list){
  dpt_cell_list <- result.list$seu_4i@assays$DptRankWindow@data %>% apply(1, function(x){names(x[x==1])})
  df_metadata_dpt <- map(dpt_cell_list, function(x){result.list$seu_4i@meta.data[x,]}) %>%
    bind_rows(.id='dpt_rank') %>% mutate(dpt_rank = as.numeric(dpt_rank), organoid = as.numeric(organoid))
  return(df_metadata_dpt)
}

map_predicted_features_to_spatial_trajectory <- function(df_dpt, features, df_features){
  df_dpt %>% left_join(df_features %>% select(c('age','meta_cluster',features)), by=c('age','meta_cluster')) %>%
  pivot_longer(cols=features, values_to = 'value', names_to = 'feature') %>%
  group_by(feature,age) %>%
  mutate(value=rescale(value)) %>% ungroup()
}

map_features_to_spatial_trajectory <- function(df_dpt, features, df_features, count_amplification=10){
  df_dpt %>% left_join(df_features %>% select(c('age','meta_cluster',features)), by=c('age','meta_cluster')) %>%
  pivot_longer(cols=features, values_to = 'value', names_to = 'feature') %>%
  group_by(feature,age) %>%
  mutate(value=rescale(value)) %>%
  group_by(feature) %>%
  mutate(density = round(value * count_amplification)) %>%
  select(distance_to_contour, dpt_rank, density, feature) %>%
  uncount(density) %>% ungroup()
}

get_norm_density <- function(x, df){
    df <- df %>% filter(dpt_rank == x)
    if (nrow(df)==0){
      return(tibble(position=0:1000,intensity=0))
    } else{
      density_est <- df %>% pull(distance_to_contour)%>% stats::density(from = 0, to=1000, n=1001, bw=10)
      density_est$y <- density_est$y/sum(density_est$y)
      return(tibble(position=density_est$x,intensity=density_est$y))
    }
}

# calculate density estimate and weight by overall expression
generate_density_matrix <- function(transcript,df){
  dpt_ranks <- df %>% distinct(dpt_rank) %>% pull()
  df_summary <- df %>% group_by(dpt_rank) %>% summarise(n=n())
  df <- df %>% filter(feature==transcript)
  names(dpt_ranks) <- dpt_ranks
  df_densities <- map(dpt_ranks, get_norm_density, df=df) %>% bind_rows(.id='dpt_rank')
  df_densities <- df_densities %>% left_join(df_summary,by="dpt_rank") %>% mutate(intensity=intensity*n) %>%
    mutate(intensity_scaled=rescale(intensity)) %>% select(-n, -intensity) %>%
    pivot_wider(id_cols = dpt_rank, names_from = position, values_from = intensity_scaled)
  M <- as.matrix(df_densities %>% select(-dpt_rank))
  rownames(M) <- df_densities %>% pull(dpt_rank)
  return(M)
}

get_wedge_coordinates <- function(x, y, angle, width, height, offset=0){
  width <- width/2
  coords <- cbind(c(x+width, x+width, x-width, x-width),
                  c(y+offset, y-height,y-height, y+offset))
  coords_rotated <- rotate(coords, degrees=(angle-180)*-1, centrexy=c(x,y))
  colnames(coords_rotated) <- c('x','y')
  return(coords_rotated)
}

map_spots_to_wedge <- function(wedge_id, df_meta, df_spots, offset=0){
  df_meta <- df_meta %>% filter(id == wedge_id)
  organoid_id <- df_meta %>% distinct(organoid) %>% pull()
  df_spots <- df_spots %>% filter(organoid == organoid_id)
  coords_wedge <- get_wedge_coordinates(df_meta$x, df_meta$y, df_meta$angle, df_meta$slice_width, df_meta$window_size, offset=offset)
  df_spots$in_wedge <- inpolygon(df_spots %>% pull(x),
                                  df_spots %>% pull(y),
                                  coords_wedge[,1],
                                  coords_wedge[,2])
  df_spots <- df_spots %>% filter(in_wedge)
  return(df_spots)
}

euclidean <- function(a, b) {
  sqrt(sum((a - b)^2))
}

get_distance_to_contour <- function(x,y, df_contours){
  M <- df_contours  %>% select(x,y) %>% as.matrix()
  positions <-  cbind(x,y)
  apply(positions, 1,
        function(x){apply(M, 1, euclidean, b=x) %>% min() %>% unique()}
  )
}

calculate_distances_to_contour <- function (point, df_nuclei, df_contours){
  df_contours <- df_contours %>% filter(organoid==point)
  df_nuclei <- df_nuclei %>%
    filter(organoid==point) %>%
    mutate(distance_to_contour = get_distance_to_contour(ids = ID, x = X, y = Y, df_contours = df_contours))
  return(df_nuclei)
}

map_cells_and_distance_to_window <- function(id,df_meta, df_spots,df_contours){
   df_spots_sub <- map_spots_to_wedge(id, df_meta, df_spots)
   df_spots_contour <- map_spots_to_wedge(id, df_meta %>% mutate(window_size=round(window_size/2)), df_contours, offset=50)
   df_spots_sub <- df_spots_sub %>% mutate(distance_to_contour = get_distance_to_contour(x = x, y = y, df_contours = df_spots_contour))
   return(df_spots_sub)
}

generate_intensity_matrix <- function(transcript,df_dpt){
  M_feature <- matrix(0, nrow=df_dpt %>% distinct(dpt_rank) %>% nrow(), ncol=1001)
  rownames(M_feature) <- df_dpt %>% distinct(dpt_rank) %>% pull()
  df_dpt <- df_dpt %>% filter(feature==transcript)
  dpt_ranks <- df_dpt %>% distinct(dpt_rank) %>% pull()
  for (i in seq_along(dpt_ranks)){
    df_tmp <- df_dpt %>% filter(dpt_rank==dpt_ranks[i])
    positions <- df_tmp %>% pull(distance_int)
    values <- df_tmp %>% pull(value)
    M_feature[dpt_ranks[i],positions+1] <- values
  }
  return(M_feature)
}

smooth_and_down_sample_profile <- function(x, window_size=50, sampling_factor=2){
  x <- meanFilter(x, windowSize = window_size)
  t <- 1:length(x)
  f <- approxfun(t,x)
  x_sampled <- f(seq(from=1,to=length(x),by=sampling_factor))
  return(x_sampled)
}

map_fish_windows_to_trajectory_bulk_zones <- function(fish_window, df_fish_windows, df_4i_windows, features, max_distance=100){
  # set up variables and dataframes
  df_fish_selected <- df_fish_windows %>% filter(window==fish_window) %>%
    mutate(zone=as.character(floor(rescale(distance_to_contour, to=c(0,5), from = c(0,1000)))))
  windows_4i <- df_4i_windows %>% distinct(dpt_rank) %>% pull()
  dist_tmp <- 100000
  tmp_results <- NULL
  all_dists <- vector('list',length = length(windows_4i)) %>% set_names(windows_4i)
  # find best matching 4i window
  for(i in seq_along(windows_4i)){
    # select 4i window
    df_4i_selected <-  df_4i_windows %>% filter(dpt_rank==windows_4i[i]) %>%
      mutate(zone=as.character(floor(rescale(distance_to_contour, to=c(0,5), from = c(0,1000)))))

    zones_4i <- unique(df_4i_selected$zone) %>% sort()
    zones_fish <- unique(df_fish_selected$zone) %>% sort()

    # check zones in both windows and skip if they do not match
    if (any(!zones_4i %in% zones_fish) | any(!zones_fish %in% zones_4i)){
      next
    } else {
      # generate correlation matrix
      df_fish_bulks <- df_fish_selected %>% mutate(zone=str_c(zone,'fish', sep='_')) %>%
        select(zone,all_of(features)) %>% group_by(zone) %>% summarise_at(features,mean)
      M_fish_expression <- df_fish_bulks %>% select(-zone) %>% as.matrix()
      rownames(M_fish_expression) <- df_fish_bulks$zone
      df_4i_bulks <- df_4i_selected %>% mutate(zone=str_c(zone,'4i', sep='_')) %>%
        select(zone,all_of(features)) %>% group_by(zone) %>% summarise_at(features,mean)
      M_4i_expression <- df_4i_bulks %>% select(-zone)%>% as.matrix()
      rownames(M_4i_expression) <- df_4i_bulks$zone

      names(zones_4i) <- str_c(zones_4i,'4i',sep='_')
      names(zones_fish) <- str_c(zones_fish,'fish',sep='_')

      M_dist_expression <- rbind(M_fish_expression[names(zones_fish),],M_4i_expression[names(zones_4i),]) %>%
        dist() %>% as.matrix()

      df <- M_dist_expression %>% as_tibble(rownames = 'zone_fish')
      df <- df %>% pivot_longer(cols = colnames(df)[2:(ncol(df))], names_to = 'zone_4i') %>%
        mutate(zone_fish=str_remove(zone_fish,'_fish'), zone_4i=str_remove(zone_4i,'_4i')) %>%
        filter(zone_fish == zone_4i)

      mean_dist <- mean(df$value)
      all_dists[i] <- mean_dist
      if(is.nan(mean_dist)|is.na(mean_dist))
        next
      if(mean_dist < dist_tmp){
        tmp_results <- list(M_dist_expression = M_dist_expression,
                            dists = df, mean_dist=mean_dist,
                            df_4i_match = df_4i_selected,
                            dpt_rank=windows_4i[i],
                            organoid=unique(df_4i_selected$organoid))
        dist_tmp <- mean_dist
    }
    }
  }
  if(is.null(tmp_results)){
    return(NULL)
  } else {
    tmp_results$all_dists <- all_dists
    # select 4i window
    df_4i_selected <- tmp_results$df_4i_match

    # generate correlation matrix
    M_fish_expression <- df_fish_selected %>% select(all_of(features)) %>% as.matrix()
    rownames(M_fish_expression) <- df_fish_selected$id
    M_4i_expression <- df_4i_selected %>% select(all_of(features)) %>% as.matrix()
    rownames(M_4i_expression) <- df_4i_selected$ID
    M_cor_expression <- cor(t(M_fish_expression),t(M_4i_expression), method = 'pearson')

    # generate spatial contraint matrix
    distances_to_contour_fish <- df_fish_selected %>% pull(distance_to_contour)
    names(distances_to_contour_fish) <- df_fish_selected$id
    distances_to_contour_4i <- df_4i_selected %>% pull(distance_to_contour)
    names(distances_to_contour_4i) <- df_4i_selected$ID
    M_dist_spatial <- dist(c(distances_to_contour_fish, distances_to_contour_4i)) %>% as.matrix()
    M_dist_spatial <- M_dist_spatial[rownames(M_cor_expression),colnames(M_cor_expression)]

    # apply constrain to correlation matrix and match nuclei
    M_cor_expression[M_dist_spatial<=max_distance] <- 0
    nuclei_fish <- rownames(M_cor_expression)
    max_nuclei <- vector(mode='character', length = length(nuclei_fish))
    max_cor <- vector(mode='numeric', length = length(nuclei_fish))
    for (i in seq_along(nuclei_fish)){
       index_max <- which.max(M_cor_expression[nuclei_fish[i],])
       nucleus_max <- colnames(M_cor_expression)[index_max][1]
       max_nuclei[i] <- nucleus_max
       max_cor[i] <-M_cor_expression[nuclei_fish[i],index_max]
    }
    df <- tibble(id_fish=nuclei_fish, id_4i=max_nuclei, organoid_4i=tmp_results$organoid, max_cor=max_cor, fish_window=fish_window, dpt_rank=tmp_results$dpt_rank)
    tmp_results$df_matched_nuclei <- df
    return(tmp_results)
  }

}
map_fish_windows_to_trajectory_bulk <- function(fish_window, df_fish_windows, df_4i_windows, features){
  # set up variables and dataframes
  df_fish_selected <- df_fish_windows %>% filter(window==fish_window)
  windows_4i <- df_4i_windows %>% distinct(dpt_rank) %>% pull()
  dist_tmp <- 100000
  tmp_results <- NULL
  all_dists <- vector('list',length = length(windows_4i)) %>% set_names(windows_4i)
  # find best matching 4i window
  for(i in seq_along(windows_4i)){
    # select 4i window
    df_4i_selected <-  df_4i_windows %>% filter(dpt_rank==windows_4i[i])
    # generate correlation matrix
    df_fish_bulks <- df_fish_selected %>% select(all_of(features)) %>% summarise_at(features,mean)
    M_fish_expression <- df_fish_bulks  %>% as.matrix()
    df_4i_bulks <- df_4i_selected %>% select(all_of(features))  %>% summarise_at(features,mean)
    M_4i_expression <- df_4i_bulks %>% as.matrix()

    M_dist_expression <- rbind(M_fish_expression,M_4i_expression) %>% dist() %>% as.matrix()
    M_dist_expression <- M_dist_expression[zones_fish, zones_4i]
    dists <- diag(M_dist_expression)
    mean_dist <- mean(dists)
    all_dists[i] <- mean_dist
    if(mean_dist < dist_tmp){
      tmp_results <- list(M_dist_expression = M_dist_expression,
                          dists = dists, mean_dist=mean_dist,
                          df_4i_match = df_4i_selected,
                          dpt_rank=windows_4i[i],
                          organoid=unique(df_4i_selected$organoid))
      dist_tmp <- mean_dist
    }
  }
  if(is.null(tmp_results)){
    return(NULL)
  } else {
    tmp_results$all_dists <- all_dists
    return(tmp_results)
  }

}




map_fish_windows_to_trajectory <- function(fish_window, df_fish_windows, df_4i_windows, features, max_distance=100){
  # set up variables and dataframes
  df_fish_selected <- df_fish_windows %>% filter(window==fish_window)
  windows_4i <- df_4i_windows %>% distinct(dpt_rank) %>% pull()
  cor_tmp <- -1
  tmp_results <- NULL
  max_cor_results <- vector('list',length = length(windows_4i)) %>% set_names(windows_4i)
  # find best matching 4i window
  for(i in seq_along(windows_4i)){
    # select 4i window
    df_4i_selected <-  df_4i_windows %>% filter(dpt_rank==windows_4i[i])

    # generate correlation matrix
    M_fish_expression <- df_fish_selected %>% select(all_of(features)) %>% as.matrix()
    rownames(M_fish_expression) <-df_fish_selected$id
    M_4i_expression <- df_4i_selected %>% select(all_of(features)) %>% as.matrix()
    rownames(M_4i_expression) <- df_4i_selected$ID
    M_cor_expression <- cor(t(M_fish_expression),t(M_4i_expression), method = 'pearson')

    # generate spatial contraint matrix
    distances_to_contour_fish <- df_fish_selected %>% pull(distance_to_contour)
    names(distances_to_contour_fish) <- df_fish_selected$id
    distances_to_contour_4i <- df_4i_selected %>% pull(distance_to_contour)
    names(distances_to_contour_4i) <- df_4i_selected$ID
    M_dist_spatial <- dist(c(distances_to_contour_fish, distances_to_contour_4i)) %>% as.matrix()
    M_dist_spatial <- M_dist_spatial[rownames(M_cor_expression),colnames(M_cor_expression)]

    # apply constrain to correlation matrix and summarise results
    M_cor_expression[M_dist_spatial<=max_distance] <- 0
    max_cors <- apply(M_cor_expression, 1 ,max, na.rm=1)
    mean_cor <- mean(max_cors)
    max_cor_results[[as.character(windows_4i[i])]] <- max_cors
    if(mean_cor > cor_tmp){
      tmp_results <- list(M_cor_expression = M_cor_expression,
                          M_dist_spatial = M_dist_spatial,
                          max_cors = max_cors, mean_cor=mean_cor,
                          dpt_rank=windows_4i[i],
                          organoid=unique(df_4i_selected$organoid))
      cor_tmp <- mean_cor
    }
  }
  if(is.null(tmp_results)){
    return(NULL)
  } else {
    tmp_results$all_max_cors <- max_cor_results
    # match nuclei
    M_cor_expression <- tmp_results$M_cor_expression
    nuclei_fish <- rownames(M_cor_expression)
    max_nuclei <- vector(mode='character', length = length(nuclei_fish))
    max_cor <- vector(mode='numeric', length = length(nuclei_fish))
    for (i in seq_along(nuclei_fish)){
       index_max <- which.max(M_cor_expression[nuclei_fish[i],])
       nucleus_max <- colnames(M_cor_expression)[index_max][1]
       max_nuclei[i] <- nucleus_max
       max_cor[i] <-M_cor_expression[nuclei_fish[i],index_max]
    }
    df <- tibble(id_fish=nuclei_fish, id_4i=max_nuclei, organoid_4i=tmp_results$organoid, max_cor=max_cor, fish_window=fish_window, dpt_rank=tmp_results$dpt_rank, mean_cor=tmp_results$mean_cor)
    tmp_results$df <- df
    return(tmp_results)
  }

}

match_laminar_nuclei <- function (fish_window, df_fish_windows, df_4i_windows, df_mapped_windows, features){
  df_fish_selected <- df_fish_windows %>% filter(window==fish_window)
  selected_dpt_rank <- df_mapped_windows %>% filter(window==fish_window) %>% distinct(dpt_rank)%>% pull(dpt_rank)
  df_4i_selected <-  df_4i_windows %>% filter(dpt_rank==as.numeric(selected_dpt_rank))

  M_fish <- df_fish_selected %>% select(features) %>% as.matrix()
  rownames(M_fish) <-df_fish_selected$id

  M_4i <- df_4i_selected %>% select(features) %>% as.matrix()
  rownames(M_4i) <- df_4i_selected$ID

  M_cor <- cor(t(M_fish),t(M_4i), method = 'spearman')

  # TODO: add spatial constrains
  nuclei_fish <- rownames(M_cor)
  max_nuclei <- vector(mode='character', length = length(nuclei_fish))
  max_cor <- vector(mode='numeric', length = length(nuclei_fish))
  for (i in seq_along(nuclei_fish)){
     index_max <- which.max(M_cor[nuclei_fish[i],])
     nucleus_max <- colnames(M_cor)[index_max][1]
     max_nuclei[i] <- nucleus_max
     max_cor[i] <-M_cor[nuclei_fish[i],index_max]
  }
  df <- tibble(id_fish=nuclei_fish, id_4i=as.numeric(max_nuclei), max_cor=max_cor, fish_window=fish_window, dpt_rank=selected_dpt_rank)
  return(df)
}

extract_dpt_rank_and_4i_intensities <- function(result.list){
  dpt_cell_list <- result.list$seu_4i@assays$DptRankWindow@data %>% apply(1, function(x){names(x[x==1])})
  df_metadata_dpt <- map(dpt_cell_list, function(x){result.list$seu_4i@meta.data[x,]}) %>%
    bind_rows(.id='dpt_rank') %>% mutate(dpt_rank = as.numeric(dpt_rank), organoid = as.numeric(organoid))
  if(result.list$seu_4i@meta.data %>% distinct(age) %>% pull() == 'adult'){
     df_expression <- result.list$seu_4i@assays$RNA@scale.data %>% as.matrix() %>% t() %>% as_tibble(rownames='id')%>%left_join(df_metadata_dpt%>%as_tibble(rownames='id'),by='id')
  } else {
    df_expression <- result.list$seu_4i@assays$integrated@scale.data %>% as.matrix() %>% t() %>% as_tibble(rownames='id')%>%left_join(df_metadata_dpt%>%as_tibble(rownames='id'),by='id')
  }

  return(df_expression)
}

extract_expression_profile <- function(age_x,ID_x,organoid_x){
  cell <- list_4i_nuclei[[age_x]][['seu_4i']]@meta.data %>% filter(age == age_x, organoid == organoid_x, ID == ID_x) %>% rownames()
  if(age_x == 'adult'){
    return(list(list_4i_nuclei[[age_x]][['seu_4i']]@assays$RNA@scale.data[,cell] %>% as_tibble(rownames='gene')))
  } else {
    return(list(list_4i_nuclei[[age_x]][['seu_4i']]@assays$integrated@scale.data[,cell] %>% as_tibble(rownames='gene')))
  }
}