# Load libraries
library(tidyverse)
library(scales)
library(destiny)
library(Seurat)
library(gridExtra)
library(grid)
library(abind)
library(dtw)
library(DistatisR)
library(TSdist)
library(FBN)
library(magick)
library(plotly)

# Define functions
read_intensity_profiles <- function(dir='/links/groups/treutlein/DATA/imaging/charmel/rotator_analysis_wo_scaling_2',
                                    subset=NULL,
                                    timepoints= c("week6","week12","week18","week24","week39","adult")){
    message(paste0("Input directory: ",dir))
    message("Loading meta files and marker informations")
    # set directories and get file paths
    directories <- list.files(dir, full.names = 1)
    names(directories) <- list.files(dir)
    if(!is.null(subset)){
      directories <- directories[subset]
    }
    file_paths <- map(directories, list.files, pattern = ".csv", full.names=1)
    files <- map(file_paths, 1)
    meta_paths <- map(file_paths,2)
    image_paths <- map(file_paths, 3)

    # link position in numpy stack via file paths to antibody stain
    df_image_files <- map_df(image_paths, read_csv, col_types = cols(), .id ='organoid') %>%
      rename(stain=index) %>%
      mutate(file=str_remove(file, '/links/groups/treutlein/DATA/imaging/charmel/shiny_input/'),
             stain=as.integer(stain)) %>%
      separate(file, into = c('sample','file'), sep = '/') %>% select(stain, file) %>%
      distinct(stain, file)

    df_channels  <- tibble(`color channel`=c('red','green','far_red'),
                       channel = c('0','1','2'))

    df_stains <- read_tsv("/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/ABs_table.txt", col_types = cols())

    df_stains <- left_join(df_stains, df_channels, by='color channel') %>%
      mutate(file = str_c('channel', channel, 'cycle', sep = '_') %>% str_c(cycle,'.tif')) %>%
      select(file, channel, AB, `gene ID`) %>% bind_rows(tibble(file='channel_hoechst.tif', AB="Hoechst", `gene ID`="Hoechst")) %>%
      left_join(df_image_files, by='file')

    # select stains to be excluded
    exclude <- df_stains %>% mutate(exclude=str_detect(`gene ID`, 'excluded')) %>% filter(exclude==1) %>% pull(stain)

    # load meta informations
    df_meta <- map_df(meta_paths, read_csv, col_types = cols(), .id ='organoid')

    # read intensity profiles and join with marker names and meta information
    message("Loading intensity profiles...")
    df <- map_df(files,read_csv,col_types = cols(),.id="organoid")

    message('Excluding stains that did not work...')
    df <- df %>% filter(!(stain %in% exclude))

    message("Joining intensity profiles with meta files and matching marker information...")
    df <- df  %>%
            left_join(df_stains, by='stain') %>%
            left_join(df_meta, by=c('window','organoid'))

    # adding age information
    df_conditions <- read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt', col_types = cols())
    permutations <- df_conditions %>% pull() %>% str_remove('_.tif') %>% str_extract('AB\\d') %>% str_extract('\\d') %>% as.numeric()
    points <- df_conditions %>% pull() %>% str_extract('s\\d*') %>% str_remove('s') %>% as.numeric()

    message("Adding age information...")
    df_conditions <- df_conditions %>% pull() %>% str_remove('_AB.*') %>%
                        str_remove('_c0_') %>% str_remove('s\\d*') %>%
                        str_remove('\\(.*') %>%
                        as_tibble() %>% mutate(point = points - 1, order = permutations) %>%
                        rename(type = value) %>% select(point,type,order) %>% filter(type != 'NA') %>%
                        rename(organoid=point) %>%
                        mutate(organoid=as.character(as.integer(organoid)))

    df <- df %>% left_join(df_conditions %>% select(-order), by='organoid') %>%
            rename(age=type, gene=`gene ID`) %>%
            select("organoid","radial_position","window","intensity","gene","age", "x", "y") %>%
            filter(age %in% timepoints)

    return(df)
}

read_cluster_profiles <- function (dir='/links/groups/treutlein/DATA/imaging/charmel/rotator_analysis_cluster_2',
                                   subset=NULL,
                                   timepoints= c("week6","week12","week18","week24","week39","adult")){
   message("Loading meta files and marker informations")
    # set directories and get file paths
    directories <- list.files(dir, full.names = 1)
    names(directories) <- list.files(dir)
    if(!is.null(subset)){
      directories <- directories[subset]
    }
    file_paths <- map(directories, list.files, pattern = ".csv", full.names=1)
    files <- map(file_paths, 1)
    meta_paths <- map(file_paths,2)
    image_paths <- map(file_paths, 3)

    # load meta informations
    df_meta <- map_df(meta_paths, read_csv, col_types = cols(), .id ='organoid')

    df_image_files <- map_df(image_paths, read.csv, .id ='organoid') %>% as_tibble() %>%
        rename(stain=X) %>%
        mutate(file=str_remove(file, '/links/groups/treutlein/DATA/imaging/charmel/clustering_results/fusion_znormalized/'),
              stain=as.integer(stain)) %>%
        separate(file, into = c('sample','folder','file'), sep = '/') %>%
        distinct(organoid,stain, file) %>% mutate(file=str_remove(file,'cluster') %>% str_remove('.png'))

    # read intensity profiles and join with marker names and meta information
    message("Loading cluster profiles...")
    df <- map_df(files,read_csv,col_types = cols(),.id="organoid")

    message("Joining intensity profiles with meta files...")
    df <- df %>% left_join(df_meta, by=c('window','organoid'))
    # add cluster informations
    df <- df %>% left_join(df_image_files, by=c('organoid','stain')) %>% mutate(stain=as.numeric(file))

    # adding age information
    message("Adding age information...")
    df_conditions <- read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt', col_types = cols())
    permutations <- df_conditions %>% pull() %>% str_remove('_.tif') %>% str_extract('AB\\d') %>% str_extract('\\d') %>% as.numeric()
    points <- df_conditions %>% pull() %>% str_extract('s\\d*') %>% str_remove('s') %>% as.numeric()

    df_conditions <- df_conditions %>% pull() %>% str_remove('_AB.*') %>%
                        str_remove('_c0_') %>% str_remove('s\\d*') %>%
                        str_remove('\\(.*') %>%
                        as_tibble() %>% mutate(point = points - 1, order = permutations) %>%
                        rename(type = value) %>% select(point,type,order) %>% filter(type != 'NA') %>%
                        rename(organoid=point) %>%
                        mutate(organoid=as.character(as.integer(organoid)))

    df <- df %>% left_join(df_conditions %>% select(-order), by='organoid') %>%
            rename(age=type) %>%
            select("organoid","radial_position","stain","window","intensity","age","x","y") %>%
            filter(age %in% timepoints)

  return(df)
}


rename_cols <- function(x, df_stains){
  df_stains <-  df_stains %>% mutate(exclude=str_detect(`gene ID`, 'excluded')) %>% filter(exclude==0)
  if (as.integer(x) %in% df_stains$stain){
    return(df_stains %>% filter(stain==as.integer(x)) %>% pull(`gene ID`))
  } else {
    return(x)
  }
}
na.omit.list <- function(y) {return(y[!vapply(y, function(x) all(is.na(x)), logical(1))])}

read_radial_profiles <- function(dir='/links/groups/treutlein/DATA/imaging/charmel/rotator_analysis_nuclei',
                                 subset=NULL, df_subset_positions=NULL){
    message(paste0("Input directory: ",dir))
    message("Loading meta files and marker informations")
    # set directories and get file paths
    directories <- list.files(dir, full.names = 1)
    names(directories) <- list.files(dir)
    if(!is.null(subset)){
      directories <- directories[subset] %>% na.omit.list()
    }
    file_paths <- map(directories, list.files, pattern = ".csv", full.names=1)
    files <- map(file_paths, 1)
    image_paths <- map(file_paths, 2)

    # link position in numpy stack via file paths to antibody stain
    df_image_files <- map_df(image_paths, read_csv, col_types = cols(), .id ='organoid') %>%
      rename(stain=index) %>%
      mutate(file=str_remove(file, '/links/groups/treutlein/DATA/imaging/charmel/shiny_input/'),
             stain=as.integer(stain)) %>%
      separate(file, into = c('sample','file'), sep = '/') %>% select(stain, file) %>%
      distinct(stain, file)

    df_channels  <- tibble(`color channel`=c('red','green','far_red'),
                       channel = c('0','1','2'))

    df_stains <- read_tsv("/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/ABs_table.txt", col_types = cols())

    df_stains <- left_join(df_stains, df_channels, by='color channel') %>%
      mutate(file = str_c('channel', channel, 'cycle', sep = '_') %>% str_c(cycle,'.tif')) %>%
      select(file, channel, AB, `gene ID`) %>% bind_rows(tibble(file='channel_hoechst.tif', AB="Hoechst", `gene ID`="Hoechst")) %>%
      left_join(df_image_files, by='file')

    # select stains to be excluded
    exclude <- df_stains %>% mutate(exclude=str_detect(`gene ID`, 'excluded')) %>% filter(exclude==1) %>% pull(stain)

    # read intensity profiles and join with marker names and meta information
    message("Loading intensity profiles...")
    df <- map_df(files,read_csv,col_types = cols(),.id="organoid")

    message('Excluding stains that did not work...')
    df <- df %>% select(-as.character(exclude))
    colnames(df) <- map_chr(colnames(df),rename_cols, df_stains=df_stains)

    if(!is.null(df_subset_positions)){
      message('Subsetting for selected nuclei...')
      df <- df %>% left_join(df_subset_positions %>% select(cell_type, x,y), by=c('x','y')) %>% na.omit()
    }

    # adding age information
    df_conditions <- read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt', col_types = cols())
    permutations <- df_conditions %>% pull() %>% str_remove('_.tif') %>% str_extract('AB\\d') %>% str_extract('\\d') %>% as.numeric()
    points <- df_conditions %>% pull() %>% str_extract('s\\d*') %>% str_remove('s') %>% as.numeric()

    message("Adding age information...")
    df_conditions <- df_conditions %>% pull() %>% str_remove('_AB.*') %>%
                        str_remove('_c0_') %>% str_remove('s\\d*') %>%
                        str_remove('\\(.*') %>%
                        as_tibble() %>% mutate(point = points - 1, order = permutations) %>%
                        rename(type = value) %>% select(point,type,order) %>% filter(type != 'NA') %>%
                        rename(organoid=point) %>%
                        mutate(organoid=as.character(as.integer(organoid)))

    df <- df %>% left_join(df_conditions %>% select(-order), by='organoid') %>%
            rename(age=type)
    message('Creating cell ids...')
    df <- df %>% mutate(id=df %>% group_by(x,y,organoid,cell_type) %>% group_indices())


    features <- df_stains %>% mutate(exclude=str_detect(`gene ID`, 'excluded')) %>% filter(exclude==0) %>% pull("gene ID")
    message('Reshaping dataframe...')
    df <- df %>% pivot_longer(cols=features, names_to = 'gene', values_to = 'intensity')
    return(df)
}

read_radial_cluster_profiles <- function(dir='/links/groups/treutlein/DATA/imaging/charmel/rotator_analysis_nuclei_mtu',
                                          subset=NULL, df_subset_positions=NULL){
    message(paste0("Input directory: ",dir))
    message("Loading meta files and marker informations")
    # set directories and get file paths
    directories <- list.files(dir, full.names = 1)
    names(directories) <- list.files(dir)
    if(!is.null(subset)){
      directories <- directories[subset] %>% na.omit.list()
    }
    file_paths <- map(directories, list.files, pattern = ".csv", full.names=1)
    files <- map(file_paths, 1)
    image_paths <- map(file_paths, 2)

    df_image_files <- map_df(image_paths, read.csv, .id ='organoid') %>% as_tibble() %>%
        rename(stain=index) %>%
        mutate(file=str_remove(file, '/links/groups/treutlein/DATA/imaging/charmel/clustering_results/fusion_znormalized/'),
              stain=as.integer(stain)) %>%
        separate(file, into = c('sample','folder','file'), sep = '/') %>%
        distinct(organoid,stain, file) %>% mutate(file=str_remove(file,'cluster') %>% str_remove('.png'))

    # read intensity profiles and join with marker names and meta information
    message("Loading cluster profiles...")
    df <- map_df(files,read_csv,col_types = cols(),.id="organoid")

    if(!is.null(df_subset_positions)){
      message('Subsetting for selected nuclei...')
      df <- df %>% left_join(df_subset_positions %>% select(cell_type, x,y), by=c('x','y')) %>% na.omit()
    }

    # add cluster informations
    df <- df %>% pivot_longer(cols = as.character(0:31), values_to = 'intensity', names_to = 'stain') %>%
      mutate(stain= as.numeric(stain)) %>%
      left_join(df_image_files, by=c('organoid','stain')) %>% mutate(stain=as.numeric(file))

    # adding age information
    message("Adding age information...")
    df_conditions <- read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt', col_types = cols())
    permutations <- df_conditions %>% pull() %>% str_remove('_.tif') %>% str_extract('AB\\d') %>% str_extract('\\d') %>% as.numeric()
    points <- df_conditions %>% pull() %>% str_extract('s\\d*') %>% str_remove('s') %>% as.numeric()

    df_conditions <- df_conditions %>% pull() %>% str_remove('_AB.*') %>%
                        str_remove('_c0_') %>% str_remove('s\\d*') %>%
                        str_remove('\\(.*') %>%
                        as_tibble() %>% mutate(point = points - 1, order = permutations) %>%
                        rename(type = value) %>% select(point,type,order) %>% filter(type != 'NA') %>%
                        rename(organoid=point) %>%
                        mutate(organoid=as.character(as.integer(organoid)))

    df <- df %>% left_join(df_conditions %>% select(-order), by='organoid') %>%
            rename(age=type)
    message('Creating cell ids...')
    df <- df %>% mutate(id=df %>% group_by(x,y,organoid,cell_type) %>% group_indices()) %>% select(-file)
    return(df)
}


load_meta_files <- function(dir='/links/groups/treutlein/DATA/imaging/charmel/rotator_analysis_wo_scaling_2', type='development'){
    # set directories and get file paths
    directories <- list.files(dir, full.names = 1)
    names(directories) <- list.files(dir)

    file_paths <- map(directories, list.files, pattern = ".csv", full.names=1)
    files <- map(file_paths, 1)
    meta_paths <- map(file_paths,2)

    df_meta <- map_df(meta_paths, read_csv, col_types = cols(), .id ='organoid')
    if(type=='development'){
     df_conditions <- read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt', col_types = cols())
     permutations <- df_conditions %>% pull() %>% str_remove('_.tif') %>% str_extract('AB\\d') %>% str_extract('\\d') %>% as.numeric()
     points <- df_conditions %>% pull() %>% str_extract('s\\d*') %>% str_remove('s') %>% as.numeric()


      df_conditions <- df_conditions %>% pull() %>% str_remove('_AB.*') %>%
                        str_remove('_c0_') %>% str_remove('s\\d*') %>%
                        str_remove('\\(.*') %>%
                        as_tibble() %>% mutate(point = points - 1, order = permutations) %>%
                        rename(type = value) %>% select(point,type,order) %>% filter(type != 'NA') %>%
                        rename(organoid=point) %>%
                        mutate(organoid=as.character(as.integer(organoid)))
      df_meta %>% left_join(df_conditions %>% select(-order), by='organoid')
    }
    if(type=='conditions'){
      df_conditions <- read_tsv('/links/groups/treutlein/DATA/imaging/4i_projects/Retina_4i_PW_new/accessory_data/conditions_by_sample.txt')
      df_meta %>% separate(organoid, into = c('well','region'),remove=FALSE) %>% left_join(df_conditions, by='well') %>%
        rename(condition=Condition, additional_id=ID)
    } else {
      df_meta
    }

}


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

filter_largest_cluster <- function(df){
   # filtering positions for largest cluster by organoid
    df_label <- df %>% group_by(organoid, label) %>% summarise(n=n()) %>% mutate(keep=n==max(n))
    df <- df %>% left_join(df_label %>% select(-n), by=c('organoid','label')) %>% filter(keep == TRUE)
    return(df)
}

detect_cut_areas_old <- function(df){
  message('Detecting cut areas on the contours...')
  df <- df %>% select(organoid, position, contour) %>% distinct()
  df <- df %>% group_by(organoid, contour) %>%
    mutate(lag_x=lag(x, order_by = window),
           lead_x=lead(x, order_by = window),
           lag_y=lag(y, order_by = window),
           lead_y=lead(y, order_by = window)) %>%
    mutate(lag_x = ifelse(is.na(lag_x), x[window==max(window)],lag_x),
           lead_x = ifelse(is.na(lead_x), x[window==min(window)],lead_x),
           lag_y = ifelse(is.na(lag_y), y[window==max(window)],lag_y),
           lead_y = ifelse(is.na(lead_y), y[window==max(window)],lead_y)) %>%
    mutate(direction = abs(atan((lead_x-lag_x)/(lead_y-lag_y)))) %>%
    mutate(keep= ifelse(direction <= 0.15 | direction >= (pi/2)-0.15, FALSE, TRUE))

  organoids <- df %>% distinct(organoid) %>% pull()
  df_result <- tibble()
  for (sample in organoids) {
    df_tmp <- df %>% filter(organoid == sample) %>% arrange(window)
    contour <- vector(mode="numeric", length = nrow(df_tmp))
    df_tmp <- df_tmp %>% mutate(window_pos = seq_along(contour))
    contour_part <- 0
    for (position in seq_along(contour)) {
      cut_area <- df_tmp %>% filter(window_pos==position) %>% pull(keep)
      if (!cut_area) {
       contour[position]=contour_part
      } else {
        contour_part <- contour_part + 1
        contour[position]=contour_part
      }
    }
    df_tmp$contour_part <- contour
    df_result <- bind_rows(df_result, df_tmp)
  }
  df_result <- df_result %>%
    group_by(contour_part, organoid) %>%
    mutate(n=n()) %>% ungroup()

  return(df_result)
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


df_conditions <- read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt', col_types = cols())
permutations <- df_conditions %>% pull() %>% str_remove('_.tif') %>% str_extract('AB\\d') %>% str_extract('\\d') %>% as.numeric()
points <- df_conditions %>% pull() %>% str_extract('s\\d*') %>% str_remove('s') %>% as.numeric()


df_conditions <- df_conditions %>% pull() %>% str_remove('_AB.*') %>%
                        str_remove('_c0_') %>% str_remove('s\\d*') %>%
                        str_remove('\\(.*') %>%
                        as_tibble() %>% mutate(point = points - 1, order = permutations) %>%
                        rename(type = value) %>% select(point,type,order) %>% filter(type != 'NA') %>%
                        rename(organoid=point) %>%
                        mutate(organoid=as.character(as.integer(organoid)))


run_position_inate_pseudotime <- function(df_feature_matrix){
  count_matrix <- df_feature_matrix %>% select(-organoid, -window, -cut_area) %>% as.matrix() %>% t()
  colnames(count_matrix) <- as.character(seq(1:nrow(df_feature_matrix)))
  seu_profiles <- CreateSeuratObject(count_matrix)
  seu_profiles@meta.data <- seu_profiles@meta.data %>%
      mutate(organoid = df_feature_matrix$organoid,
             window=df_feature_matrix$window) %>%
      left_join(df_conditions %>% select(-order), by="organoid")
  seu_profiles <- seu_profiles %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>%
      FindNeighbors(dims = 1:15) %>% FindClusters(resolution = 0.5) %>%
      RunUMAP(dims=1:15)

  dm <- DiffusionMap(Embeddings(seu_profiles, "pca")[,1:50], sigma = "global")
  dpt <- DPT(dm)
  seu_profiles$dpt_rank <- rank(dpt$dpt)
  seu_profiles$dpt <- dpt$dpt
  return(seu_profiles)
}


get_distances_by_feature <- function(feature, df, feature.col=NULL, measure='dtw') {
  t_start <- Sys.time()
  message(t_start)
  message(paste('Processing stain',feature,'...'))
  if(!is.null(feature.col)) {
    df <- df %>% rename_('stain'=feature.col)
  }
  df <- df %>% filter(stain==feature) %>% select(organoid, window, radial_position, intensity_scaled) %>%
    unite(col='id',c('organoid','window')) %>% pivot_wider(id_cols = id, values_from = intensity_scaled, names_from = "radial_position")
  M <- as.matrix(df %>% select(-id))
  rownames(M) <- df$id
  M_dist <- dist(M, method="TSDistances", distance=measure, diag=TRUE, upper=TRUE) %>% as.matrix()
  return(M_dist)
}

get_all_feature_distances <- function(df, measure="ccor"){
  message("Running dynamic time warping...")
  features <- df %>% distinct(stain) %>% pull()
  names(features) <- features
  map(features, get_distances_by_feature, df=df, measure=measure)
}

run_distatis <- function(M_dist_list){
  message('Running distatis...')
  feature_cube <- abind(M_dist_list, along=3)
  distatis(feature_cube)
}

run_destiny <- function(M_dist){
  M_dist <- as.dist(M_dist)
  dimnames(M_dist) <- NULL
  dm <- DiffusionMap(data.frame(id= 1:nrow(M_dist)), distance = M_dist)
  dpt <- DPT(dm)
  return(list(dm, dpt))
}

filter_matrix <- function(M, subset, df_wedges){
  subset <- df_wedges %>% filter(type %in% subset) %>% pull(profile_id)
  M[subset,subset]
}

get_louvain_clusters <- function(transitions) {
 	graph <- igraph::graph_from_adjacency_matrix(transitions, 'undirected', weighted = TRUE)
 	as.integer(unclass(igraph::membership(igraph::cluster_louvain(graph))))
}


map_contour_to_image <- function(df, point, annotation='predicted.id',
                         scale=NULL, cell_color=cols_ct_refine, return_gg=FALSE, pt_size=1) {

  dir <- '/links/groups/treutlein/DATA/imaging/charmel/shiny_input/'
  file <- str_c(point, '/channel_hoechst.tif')
  image <- image_read(str_c(dir, file)) %>% image_normalize()
  width <- image_info(image) %>% pull(width)
  height <- image_info(image) %>% pull(height)
  if (!is.null(scale)) {
    width_scale <- round(width * scale)
    image <- image %>% image_scale(as.character(width_scale))
  }

  df_plot <- df %>% filter(organoid == point) %>% mutate(Y = height - y, X = x - 1)

  if (!(annotation %in% c('seurat_clusters', 'predicted.id'))) {
    df_plot <- df_plot %>% select(c('X','Y',annotation)) %>% rename("feature" = annotation)
  }

  # set image as background in the plot
  bg <- rasterGrob(image, width = unit(1, 'npc'), height = unit(1, 'npc'), interpolate = TRUE)


  p1 <- ggplot(df_plot, aes(x = X, y = Y, col = feature)) +
      coord_fixed() +
      #guides(col=guide_legend(title=annotation)) +
      annotation_custom(grob = bg,
                        xmin = 0,
                        xmax = width,
                        ymin = 0,
                        ymax = height) +
      scale_x_continuous(limits = c(0, width),
                         expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, height),
                         expand = c(0, 0)) +
      geom_point(size = pt_size)

  if(return_gg){
    return(p1)
  } else {
   ggplotly(p1) %>%
    layout(images = list(source = plotly::raster2uri(as.raster(image)),
                         x = 0, y = 0,
                         sizex = width,
                         sizey = height,
                         xref = 'x', yref = 'y',
                         xanchor = 'left', yanchor = 'bottom',
                         sizing = 'fill', layer = 'below'),
           xaxis = list(showgrid = F, color = 'transparent'),
           yaxis = list(showgrid = F, color = 'transparent'),
           title = list(text=annotation,
                       font = list(
               size = 18,
               color = 'white')),
           legend = list(
             font = list(
               size = 12,
               color = 'white'),
             bgcolor = 'black',
             bordercolor = 'grey',
             borderwidth = 2,
             itemsizing = 'constant'),
           plot_bgcolor = 'black',
           paper_bgcolor = 'black')
  }

}



load_cluster_annotations <- function(sort_arrange=FALSE, version=2){
  if(version==1){
    df <- tibble(`Horizontal cells` = 1,
                            `Photoreceptor nuclei` = 2,
                            `Rods non-nuclear` = 3,
                            `Neurons` = 4,
                            `Glass slide` = 5,
                            `EphB2` = 6,
                            `Notch activated RPGs / MG in primary` = 7,
                            `Neuronal subset of GC` = 8,
                            `Rods bipolar cells` = 9,
                            `Outer membrane / adherent junctions` = 10,
                            `Pax6+ GC layer and inner part of inner nuclear layer` = 11,
                            `GC in adult retina` = 12,
                            `P-bodies / NUPC+ and GLRA2+` = 13,
                            `Intermediate filaments / Astrocytes / Müller cells (GFAP+)` = 14,
                            `Misc` = 15,
                            `Müller cell bodies (CRALBP+) / end feet` = 16,
                            `Mitochondira + VGLUT / outer plexiform layer in primary` = 17,
                            `Retinal GC (TUBB3+)` = 18,
                            `B-actin` = 19,
                            `Proliferating cells` = 20,
                            `Alpha-tubulin` = 21,
                            `Serotonin` = 22,
                            `Long wave cones and IER5+` = 23,
                            `Lower EphB2` = 24,
                            `Inner and outer plexiform layers` = 25,
                            `Dynamin 1, innersegments of photoreceptors` = 26,
                            `NeuroD1+ in RBPs, amacrine and cone bipolar cells` = 27,
                            `Inner plexiform layer and short wave cones` = 28,
                            `Nuclei of bipolar cells and Müller glia (Chx10+)` = 29,
                            `Yap1+` = 30,
                            `RGC (Brn3a+)` = 31,
                            `Collagen` = 32) %>%
       gather(annotation, cluster)
    if(sort_arrange){
      df %>% mutate(cluster=as.character(cluster-1)) %>% arrange(cluster)
    } else {
      df
    }
  }
  if(version==2){
    read.csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/labelannot_fusion_znormalised.txt',
             sep='\t') %>% as_tibble() %>% rename(cluster=Cluster, annotation=IDandLabel) %>%
      mutate(cluster=as.character(cluster))
  }


}

filter_selected_profile_clusters <- function(df_cluster, type='intensity'){

  if(type=='intensity'){
   exclude.list <- list(week39=c(5,6,9,11,12,13),
                       week24=NULL,
                       week18=c(8,9),
                       week12=c(4,5,8,9),
                       week6=1)
  }

  if (type=='cluster'){
    exclude.list <- list(week39=c(5,6,9,11,12,13),
                       week24=NULL,
                       week18=c(8,9),
                       week12=c(4,5,8,9),
                       week6=1)
  }
  if (type =='manual_pw'){
    exclude.list <- list(week39 = c(2,4,6,7,8,9,10,12,13),
                         week24 = 1:13,
                         week18 = c(2,3,4,5,7,8,9,11),
                         week12 = c(1,4,6,9,10,11),
                         week6 = c(5,6,7,8,10))
  }
  if (type =='manual_pw_include_w24'){
    exclude.list <- list(week39 = c(2,4,6,7,8,9,10,12,13),
                         week24 = c(2,4:6,8,10:13),
                         week18 = c(2,3,4,5,7,8,9,11),
                         week12 = c(1,4,6,9,10,11),
                         week6 = c(5,6,7,8,10))
  }
  if (type =='manual_pw_include_w24_v2'){
    exclude.list <- list(week39 = c(4,5,7,11:14,16:18,20:22,25:29),
                         week24 = c(1,4,6:9,11:13,15),
                         week18 = c(1,2,5,7:9,13:21,23),
                         week12 = c(1,2,4,8,11,12,14,19,20),
                         week6 = c(2:4,11:13,16))
  }
  exclude.list %>% enframe() %>%
    unnest() %>%
    mutate(cluster_type=str_c(name, value, sep='_')) %>%
    pull(cluster_type) -> excluded_clusters

  df_cluster %>% mutate(cluster_type=str_c(type,louvain,sep='_')) %>%
    filter(!(cluster_type %in% excluded_clusters))
}

filter_wedges_by_bg_count <- function(df, max_frac_bg){
  df_bg <- read_csv('/nas/groups/treutlein/DATA/imaging/charmel/wedges_from_clustering/df_bg_counts.csv') %>%
  mutate(id=str_remove(wedge, '.png'),
         bg_frac = bg_count / (1000*100))

  df <- df %>% left_join(df_bg, by='id')

  df <- df %>% filter(bg_frac <= max_frac_bg)

  return(df)
}

