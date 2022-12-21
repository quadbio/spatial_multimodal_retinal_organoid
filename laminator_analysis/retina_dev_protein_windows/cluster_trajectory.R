library(tidyverse)
library(scales)
library(furrr)
library(Matrix)
library(Seurat)
library(presto)
library(RColorBrewer)
# load 4i nuclei datasets
message('Loading datasets...')
list_4i_nuclei <- readRDS('/local2/USERS/charmel/integrated_seurat_4i_scRNA_scATAC_rev.rds')
names(list_4i_nuclei) <- c('adult','week12','week18','week24','week39','week6')

# get space-time density matrices
extract_celltypes_by_dpt_rank <- function(result.list){
  dpt_cell_list <- result.list$seu_4i@assays$DptRankWindow@data %>% apply(1, function(x){names(x[x==1])})
  df_metadata_dpt <- map(dpt_cell_list, function(x){result.list$seu_4i@meta.data[x,]}) %>%
    bind_rows(.id='dpt_rank') %>% mutate(dpt_rank = as.numeric(dpt_rank), organoid = as.numeric(organoid))
  return(df_metadata_dpt)
}
get_spatial_feature_density <- function(feature, df_dpt, df_features, count_amplification=10, grid_size=c(771.25,250)){
  df <- df_dpt %>% left_join(df_features[,c('age','meta_cluster', feature)], by=c('age','meta_cluster')) %>%
    rename(value=feature) %>%
    mutate(value=rescale(value)) %>%
    mutate(density = round(value * count_amplification)) %>%
    select(distance_to_contour, dpt_rank, density) %>%
    uncount(density)
  if(nrow(df)<=10) {
    return(NULL)
  } else {
    bandwidths <- c(MASS::bandwidth.nrd(df$dpt_rank), MASS::bandwidth.nrd(df$distance_to_contour))
    if(sum(bandwidths > 0)==2){
      MASS::kde2d(df$dpt_rank, df$distance_to_contour, n=grid_size)
    } else {
      return(NULL)
    }
  }
}

get_spatial_celltype_density <- function(ct, df_dpt, grid_size = c(771,250)){
  df <- df_dpt %>% filter(cell_type == ct)
  MASS::kde2d(df$dpt_rank, df$distance_to_contour, n=grid_size)
}

spatial_map_to_tidy <- function(x){
  if(is.null(x$z)){
    return(tibble())
  }  else {
    reshape2::melt(x$z, c("x", "y"), value.name = "z") %>% as_tibble()
  }
}

df_dpt <- map(list_4i_nuclei, extract_celltypes_by_dpt_rank)  %>%
  bind_rows() %>% as_tibble() %>% filter(distance_to_contour <= 1000)

df_meta_cluster_expression <- map(list_4i_nuclei, magrittr::extract, 'meta_cluster_expression') %>% unlist(recursive = FALSE) %>%
  bind_rows(.id='age') %>% mutate(age=str_remove(age, '.meta_cluster_expression') %>% str_remove('_'))

features <- map(list_4i_nuclei,function(x){rownames(x$seu_rna_atac@assays$RNA@data)}) %>% unlist() %>% unique()
names(features) <- features
rm(list_4i_nuclei)
message('Setting up the future...')
options(future.globals.maxSize= 891289600*2)
future::plan('multisession', workers = 30, gc=TRUE)
message('Calculating the future of spatial densities...')
density_maps <- future_map(features, get_spatial_feature_density, df_dpt=df_dpt, df_features=df_meta_cluster_expression, grid_size=c(308.5,100))
cell_types <- df_dpt %>% distinct(cell_type) %>% pull()
names(cell_types) <- cell_types
density_maps_ct <- map(cell_types, get_spatial_celltype_density, df_dpt=df_dpt, grid_size=c(308.5,100))
future::plan('sequential')

message('Pivoting density maps...')
df_density_maps <- map(density_maps, spatial_map_to_tidy) %>% bind_rows(.id = 'feature') %>%
  pivot_wider(names_from = feature, values_from = z, id_cols = c('x','y')) %>%
  mutate(id=str_c(x,y,sep='_')) %>%
  select(c('id','x','y',features))

df_density_maps_ct <- map(density_maps_ct, spatial_map_to_tidy) %>% bind_rows(.id = 'cell_type') %>%
  pivot_wider(names_from = cell_type, values_from = z, id_cols = c('x','y')) %>%
  mutate(id=str_c(x,y,sep='_')) %>%
  select(c('id','x','y',cell_types))

# abusing Seurat for 'pixelwise' clustering
message('Building Seurat object...')
M <- df_density_maps %>% select(features) %>% t() %>% as.matrix(sparse=TRUE)
colnames(M) <- df_density_maps$id
M <- as(M,'dgCMatrix')
seu <- CreateSeuratObject(M)
seu@meta.data <- seu@meta.data %>% mutate(x=df_density_maps$x,
                                          y=df_density_maps$y)

# add cell type assay
M <- df_density_maps_ct %>% select(cell_types) %>% t() %>% as.matrix(sparse=TRUE)
colnames(M) <- df_density_maps_ct$id
seu[['CT']]<- CreateAssayObject(data = M)



message('Clustering space and time...')
seu <- seu %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.3)

message('Saving results...')
saveRDS(list(seu=seu,df_dpt=df_dpt,df_expression=df_meta_cluster_expression,density_maps=density_maps, density_maps_ct=density_maps_ct),'/local2/USERS/charmel/spatial_densities.rds')
message('Done.')

