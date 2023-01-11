# set working directory
setwd('integration_4i_sc_rna_atac')

# set result directory
dir_results <- 'data/processed/4i/seq_integration/'
dir.create(dir_results)

# load functions
source('utils.R')

# load meta files
df_conditions <- read_csv('data/raw/4i/df_conditions.csv')
df_meta <- read_csv('data/raw/4i/df_meta.csv')

# load scRNA-seq dataset
message('Loading scRNA-seq data...')
seu_rna <- readRDS('data/processed/rna-atac/metacells_RNA_ATAC.seurat.rds')

message('Rescaling data for all RNA features...')
seu_rna <- seu_rna %>% ScaleData(features=.@assays$RNA@data %>% rownames())

# rename intermediate cell_types
seu_rna@meta.data <- seu_rna@meta.data %>% mutate(major_ct=as.character(major_ct)) %>%
  mutate(major_ct=ifelse(is.na(major_ct),'Intermediate', major_ct))

# import 4i nuclei into seurat objects
message('Loading 4i nuclei datasets...')
points <- df_conditions %>% pull(point) %>% as.numeric()
names(points) <- as.character(points)

# setup of the future
options(future.globals.maxSize= 891289600*20)
future::plan('multisession', workers = 20, gc=TRUE)

message('Loading 4i nuclei...')
seu_objects_4i <- readRDS('data/processed/4i/4i_nuclei.rds')

# run integration pipeline by timepoints
message('Running integration pipeline for each timepoint...')

message('Processing week 6...')
run_integration_pipeline(week_4i = 6, week_rna = 6,
                        seu_rna = seu_rna, seu_objects_4i = seu_objects_4i,
                        dir_results = dir_results)

message('Processing week 12...')
run_integration_pipeline(week_4i = 12, week_rna = c(11,12,13),
                        seu_rna = seu_rna, seu_objects_4i = seu_objects_4i,
                        dir_results = dir_results)

message('Processing week 18...')
run_integration_pipeline(week_4i = 18, week_rna = 18,
                        seu_rna = seu_rna, seu_objects_4i = seu_objects_4i,
                        dir_results = dir_results)

message('Processing week 24...')
run_integration_pipeline(week_4i = 24, week_rna = 24,
                        seu_rna = seu_rna, seu_objects_4i = seu_objects_4i,
                        dir_results = dir_results)

message('Processing week 39...')
run_integration_pipeline(week_4i = 39, week_rna = c(36,38,40),
                        seu_rna = seu_rna, seu_objects_4i = seu_objects_4i,
                        dir_results = dir_results)

message('Processing adult sample...')
run_integration_pipeline(week_4i = 'adult', week_rna = c(36,38,40,46),
                         seu_rna = seu_rna, multiple_seurat_list = FALSE,
                         seu_objects_4i = import_4i(point='72', log_transform = TRUE, dim = 10),
                         dir_results = dir_results)

message('Dataset integration and mapping completed.')

