library(tidyverse)
library(Seurat)
library(Signac)
library(Matrix)
library(harmony)
library(furrr)

generate_convex_hull <- function(cell_id, df){
  M <- df %>% filter(cell==cell_id) %>% select(x,y) %>% as.matrix()
  chull_indexes <- chull(M)
  df_hull <- M[chull_indexes,] %>%
    as_tibble() %>%
    mutate(cell=cell_id)
  return(df_hull)
}

process_fish <- function(sample, dim=15, log_transform=FALSE,
                         dir='data/processed/fish/baysor'){
  dir_sample <- str_c(dir,sample, 'output', sep='/')

  df <- read_csv(str_c(dir_sample,'/segmentation.csv'))

  df_stats <- read_csv(str_c(dir_sample,'/segmentation_cell_stats.csv')) %>% mutate(cell=as.character(cell))

  cells <- df %>% filter(cell != 0) %>% distinct(cell) %>% pull()
  names(cells) <- cells

  df_hulls <- map(cells, generate_convex_hull, df=df) %>%
    bind_rows(.id='cell') %>%
    mutate(cell=cell) %>%
    group_by(cell) %>%
    mutate(centroid_x=mean(x),
           centroid_y=mean(y)) %>%
    select(cell,x,y,centroid_x,centroid_y) %>% ungroup()

  df_counts <- read_tsv(str_c(dir_sample,'/segmentation_counts.tsv'))

  M <- df_counts %>% select(-gene) %>% as.matrix()
  rownames(M) <- df_counts %>% pull(gene)
  if (log_transform){
    M <- log(M+1)
  }
  seu_fish <- CreateSeuratObject(M)

  seu_fish <- seu_fish %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>%
        FindNeighbors(dims = 1:dim) %>% FindClusters(resolution = 0.5) %>%
        RunUMAP(dims=1:dim)
  seu_fish@meta.data <- seu_fish@meta.data %>% mutate(cell=rownames(.), organoid=sample) %>%
    left_join(df_stats, by='cell')
  return(list(seu=seu_fish, df_hulls=df_hulls))
}
# define function for CCA integration of 4I datasets
integrate_fish <- function(seu_objects, dim=30){
  anchors_integration <- FindIntegrationAnchors(object.list = seu_objects, dims = 1:dim,  reduction = 'cca')
  seu_integrated <- IntegrateData(anchorset = anchors_integration, dims = 1:dim) %>% ScaleData() %>%
    RunPCA() %>% FindNeighbors(dims = 1:dim) %>% FindClusters(resolution = 0.5) %>%
    RunUMAP(dims=1:dim)
  return(seu_integrated)
}

create_meta_clusters <- function(seu_rna_sub, resolution, dim, integrate.by=NULL) {
  # run integration
  if(!is.null(integrate.by)){
    seu_rna_sub <- seu_rna_sub %>% RunHarmony(group.by.vars = integrate.by)
    seu_rna_sub <- seu_rna_sub %>%
    FindNeighbors(dims = 1:dim, reduction="harmony") %>% FindClusters(resolution = resolution) %>%
    RunUMAP(dims=1:dim, reduction="harmony")
  } else {
  # rerun clustering and UMAP
  seu_rna_sub <- seu_rna_sub %>%
    FindNeighbors(dims = 1:dim) %>% FindClusters(resolution = resolution) %>%
    RunUMAP(dims=1:dim)
  }
  seu_rna_sub@meta.data <- seu_rna_sub@meta.data %>% rename(meta_cluster = seurat_clusters)
  return(seu_rna_sub)
}

average_peaks <- function(cluster, seu_rna, seu_subset) {
  cells <- seu_subset@meta.data %>% filter(meta_cluster == cluster) %>% rownames()
  average_peaks <- apply(seu_rna@assays$ATAC_stage@data[,cells], 1, mean, na.rm=1)
}

average_peaks_in_subset <- function(seu_subset){
  clusters <- seu_subset@meta.data %>% distinct(meta_cluster) %>% pull(meta_cluster)
  names(clusters) <- as.character(clusters)
  peaks.list <-  map(clusters, average_peaks, seu_rna= seu_rna, seu_subset=seu_subset)
  M <- do.call(rbind, peaks.list)
  return(M)
}

average_enrichments <- function(cluster, seu_rna, seu_subset) {
  cells <- seu_subset@meta.data %>% filter(meta_cluster == cluster) %>% rownames()
  average_enrichments <- apply(seu_rna@assays$chromvar@data[,cells], 1, mean, na.rm=1)
  return(average_enrichments)
}

average_enrichments_in_subset <- function(seu_subset){
  clusters <- seu_subset@meta.data %>% distinct(meta_cluster) %>% pull(meta_cluster)
  names(clusters) <- as.character(clusters)
  peaks.list <-  map(clusters, average_enrichments, seu_rna= seu_rna, seu_subset=seu_subset)
  M <- do.call(rbind, peaks.list)
  return(M)
}
# define label transfer function
predict_annotations <- function(seu_ref, seu_query, dim=30){
  # run cca and anchoring
  features <- seu_query@assays$RNA[[]] %>% rownames()
  anchors <- FindTransferAnchors(reference = seu_ref, query = seu_query, features = features,
                                 dims = 1:dim, reduction = 'cca')

  # add cell type labels
  predictions <- TransferData(anchorset = anchors, refdata = seu_ref@meta.data['major_ct'] %>% pull(),
                              dims = 1:dim, weight.reduction = 'cca')
  seu_query <- AddMetaData(seu_query, metadata = predictions)
  seu_query@meta.data <- seu_query@meta.data %>% rename(cell_type=predicted.id)

  # add meta cluster labels
  predictions <- TransferData(anchorset = anchors, refdata = seu_ref@meta.data['meta_cluster'] %>% pull(),
                              dims = 1:dim, weight.reduction = 'cca')
  seu_query <- AddMetaData(seu_query, metadata = predictions)
  seu_query@meta.data <- seu_query@meta.data %>% rename(meta_cluster=predicted.id)

  return(list(seu_query=seu_query, anchors=anchors))
}

# define pipeline function
run_integration_pipeline <- function(fish_samples, week_rna, seu_rna, week_fish, seu_objects_fish,
                                     dir_results=NULL, sample_rna=NULL, sample_fish=NULL){

  # filter scRNA-seq data and edit metadata
  message('Subsetting scRNA/ATAC-seq data...')
  seu_rna_sub <- seu_rna[,seu_rna$age_weeks %in% week_rna]
  if(!is.null(sample_rna)){
    seu_rna_sub <- seu_rna_sub[, sample(colnames(seu_rna_sub), size =sample_rna, replace=F)]
  }

  # run meta clustering
  message('Running high-resolution clustering...')
  seu_rna_sub <- create_meta_clusters(seu_rna_sub, resolution=30, dim=15, integrate.by = 'sample')

  # calculate average expression
  message('Averaging gene expression by meta clusters...')
  average_expression_scaled <- AverageExpression(seu_rna_sub, group.by='meta_cluster',
                                          slot = 'scale.data', assays='RNA')[[1]] %>% t() %>%
    as_tibble(rownames = 'meta_cluster')
  average_expression <- AverageExpression(seu_rna_sub, group.by='meta_cluster',
                                          slot = 'data', assays='RNA')[[1]] %>% t() %>%
    as_tibble(rownames = 'meta_cluster')
  # calculate average peaks
  message('Averaging region peaks by meta clusters...')
  average_peaks <- average_peaks_in_subset(seu_rna_sub) %>% as_tibble(rownames = 'meta_cluster')

  # calculate average motif enrichments
  message('Averaging motif enrichment scores by meta clusters...')
  average_motifs <- average_enrichments_in_subset(seu_rna_sub) %>% as_tibble(rownames = 'meta_cluster')

  # select fish datasets for selected week
  message('Selecting fish datasets for specified timepoint...')
  seu_objects <- map(seu_objects_fish[fish_samples],1)
  df_hulls <- map(seu_objects_fish[fish_samples],2) %>% bind_rows(.id='organoid')

  # integrate 4i datasets
  message('Integrating fish datasets...')
  seu_integrated <- integrate_fish(seu_objects, dim = 10)
  seu_integrated@meta.data <- seu_integrated@meta.data %>% mutate(age_weeks=week_fish)
  if(!is.null(sample_fish)){
    seu_integrated <- seu_integrated[,sample(colnames(seu_integrated), size = sample_4i, replace=F)]
  }

  # predicting cell type annotation from scRNAseq onto 4i
  message('Predicting cell type labels and meta clusters...')
  seu_integrated_query <- predict_annotations(seu_ref = seu_rna_sub, seu_query = seu_integrated, dim = 10)

  # creating result list
  results <- list(seu_fish = seu_integrated_query$seu_query,
                  seu_rna_atac = seu_rna_sub,
                  anchors = seu_integrated_query$anchors,
                  meta_cluster_expression = average_expression,
                  meta_cluster_expression_scaled = average_expression_scaled,
                  meta_cluster_peaks = average_peaks,
                  meta_cluster_motifs = average_motifs,
                  convex_hulls=df_hulls)

  # save results
  message('Saving results...')
  if(!is.null(dir_results)){
    saveRDS(results, paste0(dir_results,'week_',week_fish,'.rds'))
    message('Done. Saved files to: ', paste0(dir_results,week_fish,'.rds'))
  } else {
    return(results)
  }

}

# set result directory
dir_results <- 'data/processed/fish/seq_integration/'

# load scRNA-seq dataset
message('Loading scRNA-seq data...')
seu_rna <- readRDS('data/processed/rna-atac/metacells_RNA_ATAC.seurat.rds')

message('Rescaling data for all RNA features...')
seu_rna <- seu_rna %>% ScaleData(features=.@assays$RNA@data %>% rownames())

# rename intermediate cell_types
seu_rna@meta.data <- seu_rna@meta.data %>% mutate(major_ct=as.character(major_ct)) %>%
  mutate(major_ct=ifelse(is.na(major_ct),'Intermediate', major_ct))

# read in FISH data
message('Reading in fish data, processing convex hulls, converting seurat...')
dir_input <- 'data/processed/fish/baysor'
files <- list.files(dir_input)
names(files) <- files
list_results <- map(files, process_fish, dir=dir_input)

samples_week32 <- files[str_detect(files,'A|B')]
samples_week13 <- files[str_detect(files,'C|D')]

options(future.globals.maxSize= 891289600*20)
future::plan('multisession', workers = length(samples_week13), gc=TRUE)
run_integration_pipeline(fish_samples=samples_week13,
                         week_rna = c(11,12,13,15),
                         seu_rna=seu_rna,
                         week_fish=13,
                         seu_objects_fish=list_results,
                         dir_results=dir_results)

future::plan('multisession', workers = length(samples_week32), gc=TRUE)
run_integration_pipeline(fish_samples=samples_week32,
                         week_rna=c(30,31,32,34),
                         seu_rna=seu_rna,
                         week_fish=32,
                         seu_objects_fish=list_results,
                         dir_results=dir_results)

