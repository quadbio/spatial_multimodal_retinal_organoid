library(pracma)
library(secr)

# set working directory
setwd('integration_4i_sc_rna_atac')

# set directory
dir <- 'data/processed/4i/seq_integration'

# load functions
source('utils.R')

# load meta files
source('load_meta_files.R')

files <- list.files(dir, full.names = 1)
names(files) <- list.files(dir)


list_4i_nuclei <- map(files, readRDS)
names(list_4i_nuclei) <- c('adult','week_12','week_18','week_24','week_39','week_6')


euclidean <- function(a, b) {
  sqrt(sum((a - b)^2))
}

get_distance_to_contour <- function(ids,x,y, df_contours){
  M <- df_contours  %>% select(x,y) %>% as.matrix()
  positions <-  cbind(x,y)
  rownames(positions) <- ids
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

df_contours <- read_csv('data/raw/4i/df_contours.csv') %>%
  filter(!(organoid==72 & contour == 1))

add_distance_to_contour <- function (result.list){
  points <- result.list$seu_4i@meta.data %>% distinct(organoid) %>% pull(organoid)
  names(points) <- points
  message(result.list$seu_4i@meta.data %>% distinct(age) %>% pull())
  # setup of the future
  options(future.globals.maxSize= 891289600*2)
  future::plan('multisession', workers = length(points), gc=TRUE)

  df_distances <- future_map(points,calculate_distances_to_contour, df_nuclei = result.list$seu_4i@meta.data, df_contours = df_contours) %>%
    bind_rows() %>% select(distance_to_contour)

  result.list$seu_4i@meta.data %>% mutate(rowname=rownames(.))%>% left_join(df_distances %>% mutate(rowname=rownames(.)), by='rowname') -> df
  rownames(df) <- df$rowname
  result.list$seu_4i@meta.data <- df %>% select(-rowname)

  return(result.list)
}

# load datasets
df_pseudotime <- readRDS('/local2/USERS/charmel/pseudotime_trajectory_final/pseudotime_trajectory.rds') %>%
  mutate(dpt_rank=max(dpt_rank)-dpt_rank)
wedge_ids <- df_pseudotime %>% distinct(id) %>% pull()
names(wedge_ids) <- wedge_ids

df_meta <- readRDS('/local2/USERS/charmel/rotator_dataset_2/df_meta.rds') %>% filter(id %in% wedge_ids)


get_wedge_coordinates <- function(x, y, angle, width, height){
  width <- width/2
  coords <- cbind(c(x+width, x+width, x-width, x-width),
                  c(y, y-height,y-height, y))
  coords_rotated <- rotate(coords, degrees=(angle-180)*-1, centrexy=c(x,y))
  colnames(coords_rotated) <- c('x','y')
  return(coords_rotated)
}

map_nuclei_to_wedge <- function(wedge_id, df_meta, df_nuclei, df_pseudotime){
  df_meta <- df_meta %>% filter(id == wedge_id)
  organoid_id <- df_meta %>% pull(organoid) %>% as.numeric()
  df_nuclei <- df_nuclei %>% filter(organoid == organoid_id)
  coords_wedge <- get_wedge_coordinates(df_meta$x, df_meta$y, df_meta$angle, df_meta$slice_width, df_meta$window_size)
  df_nuclei$in_wedge <- inpolygon(df_nuclei %>% pull(X),
                                  df_nuclei %>% pull(Y),
                                  coords_wedge[,1],
                                  coords_wedge[,2])
  df_nuclei <- df_nuclei %>% filter(in_wedge) %>%
    mutate(dpt_rank = df_pseudotime %>% filter(id==wedge_id) %>% pull(dpt_rank))
  return(df_nuclei)
}

# map nuclei to their respective wedges in pseudotime
add_wedge_dpt_assay <- function(result.list, wedge_ids, df_meta,  df_pseudotime){
  mapped_nuclei <- map(wedge_ids,map_nuclei_to_wedge,df_meta, result.list$seu_4i@meta.data, df_pseudotime)

  df_pseudotime_assay <- map(mapped_nuclei, function(x){rownames_to_column(x,var = 'cell_id')}) %>% bind_rows(.id='wedge_id') %>% as_tibble() %>%
    group_by(dpt_rank, cell_id) %>% summarise(n=n())
  cells_in_window <-  df_pseudotime_assay %>% distinct(cell_id) %>% pull()
  cells_not_in_window <- rownames(result.list$seu_4i@meta.data)[!(rownames(result.list$seu_4i@meta.data) %in% cells_in_window)]

  df_pseudotime_assay <- df_pseudotime_assay %>% pivot_wider(id_cols = dpt_rank, names_from = cell_id, values_from = n) %>%
    ungroup() %>%
    mutate_all(~replace(., is.na(.), 0))
  df_pseudotime_assay[,cells_not_in_window] <- 0

  M <- df_pseudotime_assay %>% select(-dpt_rank) %>% as.matrix()
  rownames(M) <- df_pseudotime_assay %>% select(dpt_rank) %>% pull() %>% as.character()

  result.list$seu_4i[['DptRankWindow']] <- CreateAssayObject(data=M)

  return(result.list)
}

structure_meta_data <- function(result.list){
  pred_scores <- colnames(result.list$seu_4i@meta.data)[str_starts(result.list$seu_4i@meta.data %>% colnames(), 'prediction.score.')]
  result.list$seu_4i@meta.data <- result.list$seu_4i@meta.data %>% select(-pred_scores)
  return(result.list)
}

list_4i_nuclei_updated <- map(list_4i_nuclei, function(x){x  %>%
  structure_meta_data() %>% add_wedge_dpt_assay(wedge_ids = wedge_ids, df_meta = df_meta, df_pseudotime = df_pseudotime) %>%
  add_distance_to_contour()})

saveRDS(list_4i_nuclei_updated, '../data/processed/4i/seq_integration/integrated_seurat_4i_scRNA_scATAC.rds')