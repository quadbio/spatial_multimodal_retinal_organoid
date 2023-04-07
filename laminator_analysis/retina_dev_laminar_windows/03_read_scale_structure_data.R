source('laminator_analysis/utils.R')

# Define functions
load_meta_files <- function(dir = 'data/processed/4i/laminator/results_protein') {
  # set directories and get file paths
  directories <- list.files(dir, full.names = 1)
  names(directories) <- list.files(dir)
  file_paths <- map(directories, list.files, pattern = ".csv", full.names = 1)
  meta_paths <- map(file_paths, 2)
  df_meta <- map_df(meta_paths, read_csv, col_types = cols(), .id = 'organoid')
  df_conditions <- read_csv('data/raw/4i/df_conditions.csv')
  df_meta <- df_meta %>% left_join(df_conditions %>% select(-order), by = 'organoid')
  return(df_meta)
}

read_intensity_profiles <- function(dir='data/processed/4i/laminator/results_protein',
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
      mutate(file=str_remove(file, 'data/raw/4i/images/'),
             stain=as.integer(stain)) %>%
      separate(file, into = c('sample','file'), sep = '/') %>% select(stain, file) %>%
      distinct(stain, file)

    df_channels  <- tibble(`color channel`=c('red','green','far_red'),
                       channel = c('0','1','2'))

    df_stains <- read_csv('data/raw/4i/df_meta.csv')

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
    message("Adding age information...")
    df_conditions <- read_csv('data/raw/4i/df_conditions.csv')

    df <- df %>% left_join(df_conditions %>% select(-order), by='organoid') %>%
            rename(age=type, gene=`gene ID`) %>%
            select("organoid","radial_position","window","intensity","gene","age", "x", "y") %>%
            filter(age %in% timepoints)

    return(df)
}

read_cluster_profiles <- function (dir='data/processed/4i/laminator/results_mtu',
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
        mutate(file=str_remove(file, 'data/processed/4i/MTU_results'),
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
    df_conditions <- read_csv('data/raw/4i/df_conditions.csv')

    df <- df %>% left_join(df_conditions %>% select(-order), by='organoid') %>%
            rename(age=type) %>%
            select("organoid","radial_position","stain","window","intensity","age","x","y") %>%
            filter(age %in% timepoints)

    return(df)
}

# protein stains
df_intensity <- read_intensity_profiles(timepoints = c("week6","week12","week18","week24","week39", "adult")) %>%
  unite('id',c('organoid','window'),remove=FALSE) %>% z_scoring_age_sample()

message('Loading meta files and filtering windows that originate from cut areas in the images...')
df_meta <- load_meta_files() %>%
  unite('id',c('organoid','window'),remove=FALSE) %>%
  group_by(organoid, contour) %>%
  filter(n()>12,
         !(organoid=="72"&contour==1),
         !(organoid =="72"&(window <10))) %>%
  detect_cut_areas() %>%
  mutate(selected=ifelse(organoid %in% c("43","6","39"), FALSE, selected))

df_meta <- bind_rows(df_meta %>% filter(organoid != "48"),
                     df_meta %>% filter(organoid=="48") %>% mutate(selected = !(window %in% 17:167)))

selected_wedges <- df_meta %>% filter(!selected) %>% pull(id)

message('Filtering intensity profiles for selected windows...')
df_intensity <- df_intensity %>%
  filter(id %in% selected_wedges) %>%
  ungroup() %>% rename(stain=gene) %>%
  select(id, organoid, window, radial_position, stain, intensity_scaled)

dir <- 'data/processed/4i/laminator/dataset_protein/'

saveRDS(df_meta,paste0(dir,'df_meta.rds'))

features <- df_intensity %>% distinct(stain) %>% pull()
for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:',features[feature]))
  df_intensity %>% filter(stain==features[feature]) %>% write.csv(paste0(dir,features[feature],'.csv'))
}

# mtu clusters
df_intensity <- read_cluster_profiles(timepoints = c("week6","week12","week18","week24","week39","adult")) %>%
  unite('id',c('organoid','window'), remove = FALSE) %>%
  scale_cluster_profiles()

selected_wedges <- df_meta %>% filter(!selected) %>% pull(id)

message('Filtering intensity profiles for selected windows...')
df_intensity <- df_intensity %>%
  filter(id %in% selected_wedges) %>%
  ungroup() %>%
  mutate(stain=str_c("cluster_",stain)) %>%
  select(id, organoid, window, radial_position, stain, intensity_scaled)

dir <- 'data/processed/4i/laminator/dataset_mtu/'

features <- df_intensity %>% distinct(stain) %>% pull()

for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:', features[feature]))
  df_intensity %>%
    filter(stain == features[feature]) %>%
    write.csv(paste0(dir, features[feature], '.csv'))
}




