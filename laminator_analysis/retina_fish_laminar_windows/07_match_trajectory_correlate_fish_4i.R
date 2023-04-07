source('laminator_analysis/retina_fish_laminar_windows/utils_compare.R')

df_contours <- read_csv('data/processed/fish/laminator/analysis_results/df_contours.csv')
dir_dataset <- 'data/processed/4i/laminator/dataset_fish'

message('Reading and filtering meta files...')
# Load meta files
df_meta <- read_rds(str_c(dir_dataset, '/df_meta.rds')) %>%
  mutate(type = ifelse(str_detect(organoid, 'A|B'), 'week32', 'week13'))
# pixelsize 4i: 0.1625 um per pixel
# pixelsize fish: 0.138 um per pixel
# scaling factor: 0.1625/0.138 = 1.177536
scaling_factor <- 0.1625 / 0.138
df_meta <- df_meta %>% mutate(slice_width = slice_width * scaling_factor,
                              window_size = window_size * scaling_factor)

df_meta <- df_meta %>%
  group_by(organoid) %>%
  filter(window %in% selected_windows[[unique(organoid)]]) %>%
  ungroup()

message('Loading 4i nuclei dataset...')
list_4i_nuclei <- readRDS('data/processed/4i/seq_integration/integrated_seurat_4i_scRNA_scATAC.rds')
names(list_4i_nuclei) <- c('adult', 'week12', 'week18', 'week24', 'week39', 'week6')

message('Loading fish dataset...')
files <- list.files('data/processed/fish/seq_integration/', full.names = TRUE)
names(files) <- list.files('data/processed/fish/seq_integration/')
list_fish_nuclei <- map(files, readRDS)
names(list_fish_nuclei) <- c('week13', 'week32')


timepoints <- names(list_fish_nuclei) %>% set_names()
results <- vector(mode = 'list', length = length(timepoints)) %>% set_names(timepoints)
for (i in seq_along(timepoints)) {

  message(paste('Processing timepoint:', timepoints[i]))
  seu_fish <- list_fish_nuclei[[timepoints[i]]]$seu_fish

  message('Assigning fish nuclei to fish laminar window data...')
  features <- seu_fish@assays$integrated@data %>%
    rownames() %>%
    set_names()
  fish_windows <- df_meta %>%
    filter(type == timepoints[i]) %>%
    distinct(id) %>%
    pull() %>%
    set_names()
  future::plan('multisession', workers = 25, gc = TRUE)
  df_fish_windows <- future_map(fish_windows, map_cells_and_distance_to_window,
                                df_meta = df_meta,
                                df_spots = seu_fish@meta.data %>% as_tibble(rownames = 'id'),
                                df_contours = df_contours) %>%
    bind_rows(.id = 'window') %>%
    select(window, id, organoid, cell_type, distance_to_contour) %>%
    mutate(distance_to_contour = distance_to_contour / scaling_factor) %>%
    filter(distance_to_contour <= 1000)

  df_fish_windows <- df_fish_windows %>% left_join(seu_fish@assays$integrated@scale.data %>%
                                                     as.matrix() %>%
                                                     t() %>%
                                                     as_tibble(rownames = 'id'), by = 'id')

  message('Extracting imputed transcriptome metaclusters from 4i data...')
  df_meta_cluster_expression <- map(list_4i_nuclei, magrittr::extract, 'meta_cluster_expression_scaled') %>%
    unlist(recursive = FALSE) %>%
    bind_rows(.id = 'age') %>%
    mutate(age = str_remove(age, '.meta_cluster_expression_scaled') %>% str_remove('_')) %>%
    select(c('age', 'meta_cluster', features))

  message('Extracting 4i window data and merge with transcriptome data...')
  df_4i_windows <- map(list_4i_nuclei, extract_celltypes_by_dpt_rank) %>%
    bind_rows() %>%
    as_tibble() %>%
    filter(distance_to_contour <= 1000) %>%
    left_join(df_meta_cluster_expression, by = c('age', 'meta_cluster')) %>%
    filter(organoid != 72) # remove adult sample

  # Correlation of averaged cell type expression
  df_4i_average <- df_4i_windows %>%
    group_by(age, cell_type) %>%
    summarise_at(features, mean) %>%
    ungroup() %>%
    mutate(type = str_c(age, cell_type, sep = '_'))
  M_4i_average <- df_4i_average %>%
    select(-type, -age, -cell_type) %>%
    as.matrix()
  rownames(M_4i_average) <- df_4i_average$type
  df_fish_average <- df_fish_windows %>%
    group_by(cell_type) %>%
    summarise_at(features, mean) %>%
    ungroup()
  M_fish_average <- df_fish_average %>%
    select(-cell_type) %>%
    as.matrix()
  rownames(M_fish_average) <- df_fish_average$cell_type
  cor_results <- cor(t(rbind(M_fish_average, M_4i_average)))
  df_cors <- cor_results[rownames(M_fish_average), rownames(M_4i_average)] %>%
    t() %>%
    as_tibble(rownames = 'type') %>%
    separate(type, c('age', 'cell_type'), sep = '_') %>%
    pivot_longer(cols = rownames(M_fish_average) %>% set_names(NULL), names_to = 'cell_type_fish', values_to = 'cor') %>%
    filter(cell_type == cell_type_fish)
  df_4i_windows %>%
    group_by(age) %>%
    mutate(n_total = n()) %>%
    group_by(age, cell_type) %>%
    mutate(n = n(), weight = n / n_total) %>%
    select(n, age, cell_type, n_total, weight) %>%
    distinct() -> df_frac
  df_fish_windows %>%
    mutate(n_total = n()) %>%
    group_by(cell_type) %>%
    mutate(n = n(), weight = n / n_total) %>%
    select(n, cell_type, n_total, weight) %>%
    distinct() -> df_frac


  df_cors_sum <- df_cors %>%
    left_join(df_frac) %>%
    group_by(age) %>%
    summarise(mean_weighted = weighted.mean(cor, weight))

  message('Matching fish windows with 4i windows...')
  list_results <- future_map(fish_windows[1:100], map_fish_windows_to_trajectory,
                             df_fish_windows = df_fish_windows,
                             df_4i_windows = df_4i_windows,
                             features = features)
  future::plan('multisession', workers = 1, gc = TRUE)

  message('Adding all available modalities to matched nuclei and running correlation analysis...')
  df_matched_nuclei <- map(list_results, magrittr::extract, 'df') %>%
    unlist(recursive = FALSE) %>%
    bind_rows()

  df_4i_intensities <- df_matched_nuclei %>%
    mutate(id_4i = as.numeric(id_4i)) %>%
    left_join(df_4i_windows %>%
                mutate(id_4i = ID) %>%
                select(ID, organoid, age, dpt_rank, id_4i), by = c('id_4i', 'dpt_rank')) %>%
    rowwise() %>%
    mutate(expression = extract_expression_profile(age, ID, organoid)) %>%
    mutate(names = str_c(fish_window, id_fish, sep = '_')) %>%
    unnest(expression)

  # matrix measured fish spots
  M_fish <- df_fish_windows %>%
    select(features) %>%
    as.matrix()
  rownames(M_fish) <- df_fish_windows %>%
    mutate(names = str_c(window, id, sep = '_')) %>%
    pull(names)

  # matrix predicted rna transcripts 4i nuclei
  df_matched_transcriptome <- df_matched_nuclei %>%
    mutate(id_4i = as.numeric(id_4i)) %>%
    left_join(df_4i_windows %>%
                select(c(dpt_rank, ID, features)) %>%
                rename(id_4i = ID), by = c('dpt_rank', 'id_4i')) %>%
    mutate(names = str_c(fish_window, id_fish, sep = '_'))

  M_4i <- df_matched_transcriptome %>%
    select(features) %>%
    as.matrix()
  rownames(M_4i) <- df_matched_transcriptome$names

  # expression values
  features_overlap <- features[features %in% unique(df_4i_intensities$gene)]
  M_fish_protein <- M_fish[, features_overlap]

  df_cor_fish_sc_rna <- list_fish_nuclei[[timepoints[i]]]$seu_fish@meta.data %>%
    as_tibble(rownames = 'id') %>%
    left_join(list_fish_nuclei[[timepoints[i]]]$
                meta_cluster_expression_scaled, by = 'meta_cluster') %>%
    select(c('id', features)) %>%
    pivot_longer(cols = features, names_to = 'feature', values_to = 'expression_predicted') %>%
    left_join(list_fish_nuclei[[timepoints[i]]]$seu_fish@assays$integrated@scale.data %>%
                t() %>%
                as_tibble(rownames = 'id') %>%
                pivot_longer(cols = features, names_to = 'feature', values_to = 'expression_fish'), by = c('id', 'feature')) %>%
    filter(id %in% df_matched_nuclei$id_fish) %>%
    group_by(feature) %>%
    summarise(cor = cor(expression_fish, expression_predicted)) %>%
    ungroup()

  df_compare <- M_4i %>%
    as_tibble(rownames = 'id') %>%
    pivot_longer(names_to = 'feature', values_to = 'expression_4i', cols = features) %>%
    left_join(M_fish %>%
                as_tibble(rownames = 'id') %>%
                pivot_longer(names_to = 'feature', values_to = 'expression_fish', cols = features), by = c('id', 'feature')) %>%
    group_by(feature) %>%
    summarise(cor = cor(expression_4i, expression_fish)) %>%
    mutate(type = 'predicted_RNA_4i_nuclei-fish_rna') %>%
    ungroup() %>%
    bind_rows(df_cor_fish_sc_rna %>% mutate(type = 'predicted_RNA_fish_nuclei-fish_rna'))

  df_compare_2 <- df_4i_intensities %>%
    select(names, value, gene) %>%
    rename(id = names, feature = gene, expression_4i = value) %>%
    filter(feature %in% features_overlap) %>%
    left_join(M_fish_protein %>%
                as_tibble(rownames = 'id') %>%
                pivot_longer(names_to = 'feature', values_to = 'expression_fish', cols = features_overlap),
              by = c('id', 'feature')) %>%
    group_by(feature) %>%
    summarise(cor = cor(expression_4i, expression_fish)) %>%
    ungroup() %>%
    mutate(type = 'protein_4i_nuclei-fish_rna') %>%
    bind_rows(M_4i %>%
                as_tibble(rownames = 'id') %>%
                pivot_longer(names_to = 'feature', values_to = 'expression_4i', cols = features) %>%
                left_join(M_fish %>%
                            as_tibble(rownames = 'id') %>%
                            pivot_longer(names_to = 'feature', values_to = 'expression_fish', cols = features), by = c('id', 'feature')) %>%
                group_by(feature) %>%
                summarise(cor = cor(expression_4i, expression_fish)) %>%
                ungroup() %>%
                filter(feature %in% features_overlap) %>%
                mutate(type = 'predicted_RNA_4i_nuclei-fish_rna')) %>%
    bind_rows(df_cor_fish_sc_rna %>%
                filter(feature %in% features_overlap) %>%
                mutate(type = 'predicted_RNA_fish_nuclei-fish_rna')) %>%
    bind_rows(df_matched_transcriptome %>%
                select(c(id_fish, features_overlap)) %>%
                pivot_longer(names_to = 'feature', values_to = 'expression_4i', cols = features_overlap) %>%
                left_join(df_4i_intensities %>%
                            select(id_fish, value, gene) %>%
                            rename(feature = gene, expression_4i_protein = value) %>%
                            filter(feature %in% features_overlap), by = c('feature', 'id_fish')) %>%
                group_by(feature) %>%
                summarise(cor = cor(expression_4i, expression_4i_protein)) %>%
                mutate(type = 'predicted_RNA_4i_nuclei-protein') %>%
                ungroup())


  results[[timepoints[i]]] <- list(df_cor_fish_4i_rna = df_compare, df_cor_fish_4i_protein_rna = df_compare_2, matched_trajectory = list_results,
                                   df_fish_windows = df_fish_windows)
  message(paste0('Finished processing of timepoint: ', timepoints[i]))
}
results$df_4i_windows <- df_4i_windows
message('Processed all timepoints. Saving results...')
saveRDS(results, 'data/processed/fish/laminator/analysis_results/matched_trajectory_correlation.rds')
message('Done.')







