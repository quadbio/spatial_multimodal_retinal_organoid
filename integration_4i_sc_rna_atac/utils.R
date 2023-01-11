# load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(Matrix)
library(harmony)
library(furrr)
library(pracma)
library(secr)

# define function for CCA integration of 4I datasets
integrate_4i <- function(seu_objects, dim=30){
  anchors_integration <- FindIntegrationAnchors(object.list = seu_objects, dims = 1:dim,  reduction = 'cca')
  seu_integrated <- IntegrateData(anchorset = anchors_integration, dims = 1:dim) %>% ScaleData() %>%
    RunPCA() %>% FindNeighbors(dims = 1:dim) %>% FindClusters(resolution = 0.5) %>%
    RunUMAP(dims=1:dim)
  return(seu_integrated)
}

run_cca <- function(seu_ref, seu_query, dim=30) {
  features <- seu_query@assays$integrated@scale.data %>% rownames()
  cca_res <- RunCCA(seu_ref, seu_query, features=features, num.cc = dim, assay1 = 'RNA', assay2 = seu_query %>% DefaultAssay())
  return(cca_res)
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
run_integration_pipeline <- function(week_4i, week_rna, seu_rna, seu_objects_4i,
                                     multiple_seurat_list=TRUE,
                                     dir_results=NULL, sample_rna=NULL, sample_4i=NULL){
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

  # filter 4i samples for selected week
  if(multiple_seurat_list){
    message('Selecting 4i datasets for specified timepoint...')
    week_4i <- str_c('week', week_4i)
    points_week <- df_conditions %>% filter(type == week_4i) %>% pull(point) %>% as.numeric()
    names(points_week) <- as.character(points_week)
    points_week <- points_week[points_week %in% names(seu_objects_4i)]

    # select 4i datasets for selected week
    seu_objects <- seu_objects_4i[as.character(points_week)]

    # integrate 4i datasets
    message('Integrating 4i datasets...')
    seu_integrated <- integrate_4i(seu_objects, dim = 10)
    seu_integrated@meta.data <- seu_integrated@meta.data %>% mutate(age=week_4i)
  } else {
    message('Just one seurat object passed. Skipping integration...')
    seu_integrated <- seu_objects_4i
    seu_integrated@meta.data <- seu_integrated@meta.data %>% mutate(age='adult')
  }
  if(!is.null(sample_4i)){
    seu_integrated <- seu_integrated[, sample(colnames(seu_integrated), size = sample_4i, replace=F)]
  }

  # predicting cell type annotation from scRNAseq onto 4i
  message('Predicting cell type labels and meta clusters...')
  seu_integrated_query <- predict_annotations(seu_ref = seu_rna_sub, seu_query = seu_integrated, dim = 10)

  # creating result list
  results <- list(seu_4i = seu_integrated_query$seu_query,
                  seu_rna_atac = seu_rna_sub,
                  anchors = seu_integrated_query$anchors,
                  meta_cluster_expression = average_expression,
                  meta_cluster_expression_scaled = average_expression_scaled,
                  meta_cluster_peaks = average_peaks,
                  meta_cluster_motifs = average_motifs)

  # save results
  if(!is.null(dir_results)){
    saveRDS(results, paste0(dir_results,week_4i,'.rds'))
  } else {
    return(results)
  }
}

get_wedge_coordinates <- function(x, y, angle, width, height){
  width <- width/2
  coords <- cbind(c(x+width, x+width, x-width, x-width),
                  c(y, y-height,y-height, y))
  coords_rotated <- rotate(coords, degrees=(angle-180)*-1, centrexy=c(x,y))
  colnames(coords_rotated) <- c('x','y')
  return(coords_rotated)
}

map_nuclei_to_wedge <- function(wedge_id, df_meta, df_nuclei){
  df_meta <- df_meta %>% filter(id == wedge_id)
  organoid_id <- df_meta %>% pull(organoid) %>% as.numeric()
  df_nuclei <- df_nuclei %>% filter(organoid == organoid_id)
  coords_wedge <- get_wedge_coordinates(df_meta$x, df_meta$y, df_meta$angle, df_meta$slice_width, df_meta$window_size)
  df_nuclei$in_wedge <- inpolygon(df_nuclei %>% pull(X),
                                  df_nuclei %>% pull(Y),
                                  coords_wedge[,1],
                                  coords_wedge[,2])
  return(df_nuclei %>% filter(in_wedge))
}

euclidean <- function(a, b) {
  sqrt(sum((a - b)^2))
}

get_distance_to_contour <- function(ids,x,y, df_contours){
  M <- df_contours  %>% select(x,y) %>% as.matrix()
  positions <- cbind(x,y)
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

