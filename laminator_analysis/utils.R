# Load libraries
library(tidyverse)
library(scales)
library(FBN)

# define helper functions
rename_cols <- function(x, df_stains){
  df_stains <-  df_stains %>% mutate(exclude=str_detect(`gene ID`, 'excluded')) %>% filter(exclude==0)
  if (as.integer(x) %in% df_stains$stain){
    return(df_stains %>% filter(stain==as.integer(x)) %>% pull(`gene ID`))
  } else {
    return(x)
  }
}

na.omit.list <- function(y) {return(y[!vapply(y, function(x) all(is.na(x)), logical(1))])}

z_scoring_condition_sample <- function(df){
    message("Z-scoring across stain, condition and sample...")
    df_age <- df %>% group_by(condition, stain) %>% summarize(mean_condition = mean(intensity, na.rm = 1),
                                             sd_condition = sd(intensity, na.rm = 1))
    df_sample <- df %>% group_by(stain, organoid) %>% summarize(mean_sample = mean(intensity, na.rm = 1),
                                             sd_sample = sd(intensity, na.rm = 1))
    message("Merging datasets and calculating intensities...")
    df <- df  %>% left_join(df_age, by=c('stain','condition')) %>%
      left_join(df_sample, by=c('stain','organoid')) %>%
      mutate(intensity_scaled = ((intensity - mean_sample)/sd_sample) * sd_condition + mean_condition)

    message('Scaling data...')
    df <- df %>% group_by(stain) %>% mutate(intensity_scaled = rescale(intensity_scaled, to = c(0,1))) %>%
      ungroup()

    return(df)
}

z_scoring_age_sample <- function(df){
    message("Z-scoring across stain, age and sample...")
    df_age <- df %>% group_by(age, gene) %>% summarize(mean_age = mean(intensity),
                                             sd_age = sd(intensity))
    df_sample <- df %>% group_by(gene, organoid) %>% summarize(mean_sample = mean(intensity),
                                             sd_sample = sd(intensity))
    message("Merging datasets and calculating intensities...")
    df <- df  %>% left_join(df_age, by=c('gene','age')) %>%
      left_join(df_sample, by=c('gene','organoid')) %>%
      mutate(intensity_scaled = ((intensity - mean_sample)/sd_sample) * sd_age + mean_age)

    message('Scaling data...')
    df <- df %>% group_by(gene) %>% mutate(intensity_scaled = rescale(intensity_scaled, to = c(0,1))) %>%
      ungroup()

    return(df)
}

scale_cluster_profiles <- function(df) {
    message('Scaling data...')
    df %>% group_by(stain) %>%
      mutate(intensity_scaled = rescale(intensity, to = c(0,1))) %>%
      ungroup()
}

detect_cut_areas <- function(df, direction_thresh=0.125){
 df %>% group_by(organoid, contour) %>%
  arrange(organoid, contour, position) %>%
  mutate(direction = get_directions(x,y),
         score = detect_straight_lines(direction),
         keep_direction = ifelse(direction <= direction_thresh | direction >= (pi/2)-direction_thresh,
                      FALSE, TRUE),
         contour_part = detect_consecutive_areas(keep_direction)) %>%
  group_by(organoid, contour, contour_part) %>% mutate(length = n(), keep_score = sum(keep_direction)) %>% ungroup() %>%
  mutate(straight_area= (length > 6 & keep_score == 0)) %>%
  group_by(organoid, contour) %>%
  mutate(straight_part = detect_consecutive_areas(straight_area)) %>%
  group_by(organoid, contour,straight_part) %>%
  mutate(length_straight_part = n()) %>% ungroup() %>%
  mutate(add = (length_straight_part < 15 & straight_area == FALSE),
         selected = straight_area | add)
}

detect_straight_lines <- function (angles){
  x <- angles %>% abs() %>% medianFilter(windowSize = 5)
  res <- vector('numeric', length = length(x))
  for (i in seq_along(x)){
    res[i] <- lead(x,i-1)[1:3] %>% sd()
  }
  res_rev <- vector('numeric', length = length(x))
  for (i in seq_along(x)){
    res_rev[i] <- lead(rev(x),i-1)[1:5] %>% sd()
  }
  scores <- tibble(forward = res, reverse = rev(res_rev), position = seq_along(res)) %>%
    group_by(position) %>% summarise(score = mean(c(forward, reverse), na.rm = TRUE)) %>% pull(score)
  return(scores)
}

get_directions <- function (x, y){
  res <- vector('numeric', length = length(x))
  for (i in seq_along(x)){
    x_pair <- lead(x,i-1)[1:2]
    y_pair <- lead(y,i-1)[1:2]
    res[i] <- abs(atan((x_pair[1]-x_pair[2])/(y_pair[1]-y_pair[2])))
  }
  res[length(res)] <- res[length(res)-1]
  return(res)
}

detect_consecutive_areas <- function(x){
  contour <- vector(mode = 'numeric', length = length(x))
  contour_part <- 1
  state_contour_part <- x[1]
  for (i in seq_along(x)) {
      state_position <- x[i]
      if (state_position==state_contour_part) {
       contour[i] <- contour_part
      } else {
        contour_part <- contour_part + 1
        contour[i] <- contour_part
        state_contour_part <- x[i]
      }

  }
  return(contour)
}

get_louvain_clusters <- function(transitions) {
 	graph <- igraph::graph_from_adjacency_matrix(transitions, 'undirected', weighted = TRUE)
 	as.integer(unclass(igraph::membership(igraph::cluster_louvain(graph))))
}
