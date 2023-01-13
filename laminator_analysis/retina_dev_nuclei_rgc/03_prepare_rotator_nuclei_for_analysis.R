source('laminator_analysis/utils.R')

# load 4i nuclei datasets
list_4i_nuclei <- readRDS('data/processed/4i/seq_integration/integrated_seurat_4i_scRNA_scATAC.rds')
names(list_4i_nuclei) <- c('adult', 'week12', 'week18', 'week24', 'week39', 'week6')

# protein stains
read_radial_profiles <- function(dir = 'data/processed/4i/laminator/results_nuclei',
                                 subset = NULL, df_subset_positions = NULL) {
  message(paste0("Input directory: ", dir))
  message("Loading meta files and marker informations")
  # set directories and get file paths
  directories <- list.files(dir, full.names = 1)
  names(directories) <- list.files(dir)
  if (!is.null(subset)) {
    directories <- directories[subset] %>% na.omit.list()
  }
  file_paths <- map(directories, list.files, pattern = ".csv", full.names = 1)
  files <- map(file_paths, 1)
  image_paths <- map(file_paths, 2)

  # link position in numpy stack via file paths to antibody stain
  df_image_files <- map_df(image_paths, read_csv, col_types = cols(), .id = 'organoid') %>%
    rename(stain = index) %>%
    mutate(file = str_remove(file, 'data/raw/4i/images/'),
           stain = as.integer(stain)) %>%
    separate(file, into = c('sample', 'file'), sep = '/') %>%
    select(stain, file) %>%
    distinct(stain, file)

  df_channels <- tibble(`color channel` = c('red', 'green', 'far_red'),
                        channel = c('0', '1', '2'))

  df_stains <- read_csv('data/raw/4i/df_meta.csv')

  df_stains <- left_join(df_stains, df_channels, by = 'color channel') %>%
    mutate(file = str_c('channel', channel, 'cycle', sep = '_') %>% str_c(cycle, '.tif')) %>%
    select(file, channel, AB, `gene ID`) %>%
    bind_rows(tibble(file = 'channel_hoechst.tif', AB = "Hoechst", `gene ID` = "Hoechst")) %>%
    left_join(df_image_files, by = 'file')

  # select stains to be excluded
  exclude <- df_stains %>%
    mutate(exclude = str_detect(`gene ID`, 'excluded')) %>%
    filter(exclude == 1) %>%
    pull(stain)

  # read intensity profiles and join with marker names and meta information
  message("Loading intensity profiles...")
  df <- map_df(files, read_csv, col_types = cols(), .id = "organoid")

  message('Excluding stains that did not work...')
  df <- df %>% select(-as.character(exclude))
  colnames(df) <- map_chr(colnames(df), rename_cols, df_stains = df_stains)

  if (!is.null(df_subset_positions)) {
    message('Subsetting for selected nuclei...')
    df <- df %>%
      left_join(df_subset_positions %>% select(cell_type, x, y), by = c('x', 'y')) %>%
      na.omit()
  }

  # adding age information
  message("Adding age information...")
  df_conditions <- read_csv('data/raw/4i/df_conditions.csv')

  df <- df %>%
    left_join(df_conditions %>% select(-order), by = 'organoid') %>%
    rename(age = type)

  message('Creating cell ids...')
  df <- df %>% mutate(id = df %>%
    group_by(x, y, organoid, cell_type) %>%
    group_indices())

  features <- df_stains %>%
    mutate(exclude = str_detect(`gene ID`, 'excluded')) %>%
    filter(exclude == 0) %>%
    pull("gene ID")
  message('Reshaping dataframe...')
  df <- df %>% pivot_longer(cols = features, names_to = 'gene', values_to = 'intensity')
  return(df)
}

df_subset <- map(list_4i_nuclei, function(x) {
  x$seu_4i@meta.data %>%
    mutate(organoid = as.character(organoid)) %>%
    as_tibble() }) %>%
  bind_rows() %>%
  filter(age == 'week12', cell_type == 'RGC') %>%
  mutate(x = round(X), y = round(Y))

subset <- df_subset %>%
  distinct(organoid) %>%
  pull()

df_intensity <- read_radial_profiles(subset = subset, df_subset_positions = df_subset) %>%
  z_scoring_age_sample()

dir <- 'data/processed/4i/laminator/dataset_nuclei/'

features <- df_intensity %>%
  distinct(gene) %>%
  pull()

for (feature in seq_along(features)) {
  message(paste('Writing dataset for stain:', features[feature]))
  df_intensity %>%
    filter(gene == features[feature]) %>%
    write.csv(paste0(dir, features[feature], '.csv'))
}

# for MTUs
read_radial_cluster_profiles <- function(dir = 'data/processed/4i/laminator/results_nuclei_mtu',
                                         subset = NULL, df_subset_positions = NULL) {
  message(paste0("Input directory: ", dir))
  message("Loading meta files and marker informations")
  # set directories and get file paths
  directories <- list.files(dir, full.names = 1)
  names(directories) <- list.files(dir)
  if (!is.null(subset)) {
    directories <- directories[subset] %>% na.omit.list()
  }
  file_paths <- map(directories, list.files, pattern = ".csv", full.names = 1)
  files <- map(file_paths, 1)
  image_paths <- map(file_paths, 2)

  df_image_files <- map_df(image_paths, read.csv, .id = 'organoid') %>%
    as_tibble() %>%
    rename(stain = index) %>%
    mutate(file = str_remove(file, 'data/processed/4i/MTU_results'),
           stain = as.integer(stain)) %>%
    separate(file, into = c('sample', 'folder', 'file'), sep = '/') %>%
    distinct(organoid, stain, file) %>%
    mutate(file = str_remove(file, 'cluster') %>% str_remove('.png'))

  # read intensity profiles and join with marker names and meta information
  message("Loading cluster profiles...")
  df <- map_df(files, read_csv, col_types = cols(), .id = "organoid")

  if (!is.null(df_subset_positions)) {
    message('Subsetting for selected nuclei...')
    df <- df %>%
      left_join(df_subset_positions %>% select(cell_type, x, y), by = c('x', 'y')) %>%
      na.omit()
  }

  # add cluster informations
  df <- df %>%
    pivot_longer(cols = as.character(0:31), values_to = 'intensity', names_to = 'stain') %>%
    mutate(stain = as.numeric(stain)) %>%
    left_join(df_image_files, by = c('organoid', 'stain')) %>%
    mutate(stain = as.numeric(file))

  # adding age information
  message("Adding age information...")
  df_conditions <- read_csv('data/raw/4i/df_conditions.csv')
  df <- df %>%
    left_join(df_conditions %>% select(-order), by = 'organoid') %>%
    rename(age = type)

  message('Creating cell ids...')
  df <- df %>%
    mutate(id = df %>%
      group_by(x, y, organoid, cell_type) %>%
      group_indices()) %>%
    select(-file)
  return(df)
}

dir <- 'data/processed/4i/laminator/dataset_nuclei_mtu/'

df_intensity <- read_radial_cluster_profiles(subset = subset, df_subset_positions = df_subset) %>%
  scale_cluster_profiles()

df_intensity <- df_intensity %>%
  ungroup() %>%
  mutate(stain = str_c("cluster_", stain)) %>%
  group_by(stain) %>%
  mutate(intensity = rescale(intensity)) %>%
  ungroup()


features <- df_intensity %>%
  distinct(stain) %>%
  pull()

for (feature in seq_along(features)) {
  message(paste('Writing dataset for stain:', features[feature]))
  df_intensity %>%
    filter(stain == features[feature]) %>%
    write.csv(paste0(dir, features[feature], '.csv'))
}
