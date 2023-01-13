if (! require(AnnotationHub)){
  BiocManager::install("AnnotationHub")
  library(AnnotationHub)
}
if (! require(GenomeInfoDb)){
  BiocManager::install("GenomeInfoDb")
  library(GenomeInfoDb)
}
if (! require(Seurat)){
  install.packages("Seurat")
  library(Seurat)
}
if (! require(Signac)){
  install.packages("Signac")
  library(Signac)
}
if (! require(anndata)){
  install.packages("anndata")
  library(anndata)
}


# retrieve Ensembl annotation
retrieve_ensbd <- function(species = "EnsDb.Hsapiens", id = NULL, version = NULL, seqlevelsStyle = "UCSC", force = FALSE){
  ah <- AnnotationHub()
  ensdbs <- query(ah, c("EnsDb.Hsapiens"))
  
  if (!is.null(id) & sum(ensdbs$ah_id == id)==0) id <- NULL
  if (!is.null(version) & length(grep(paste0(" ",version," EnsDb"), ensdbs$title))==0) version <- NULL
  
  if (is.null(id) & is.null(version)){
    idx <- length(ensdbs$ah_id)
    id <- ensdbs$ah_id[idx]
    message(paste0("retrieving the newest annotation: ID:",id," Title:\"",ensdbs$title[idx],"\""))
  } else if(is.null(id)){
    idx <- grep(paste0(" ",version," EnsDb"), ensdbs$title)
    id <- ensdbs$ah_id[idx]
    message(paste0("retrieving the requested annotation: ID:",id," Title:\"",ensdbs$title[idx],"\""))
  } else{
    message(paste0("retrieving the requested annotation: ID:",id," Title:\"",ensdbs$title[which(ensdbs$ah_id==id)],"\""))
  }
  
  ensdb <- ensdbs[[id, force = force]]
  seqlevelsStyle(ensdb) <- seqlevelsStyle
  
  return(ensdb)
}

# create a Seurat object for scATAC-seq data
create_atac_seurat <- function(path, sample = NULL, assay = "peaks", project= "SeuratProject", annotations = NULL, annot_ens_version = NULL, force_ens_dl = TRUE, verbose=T){
  if (verbose)
    message(" % reading count matrix from the h5 file...")
  counts <- Read10X_h5(paste0(path, "/outs/filtered_peak_bc_matrix.h5"))
  
  if (verbose)
    message(" % reading metadata...")
  meta <- data.frame(read.csv(paste0(path, "/outs/singlecell.csv"), header=T, row.names=1)[colnames(counts),],
                     row.names = colnames(counts))
  if (!is.null(sample))
    meta <- data.frame(sample = sample, meta, row.names = rownames(meta))
  
  if (verbose)
    message(" % creating chromatin assay...")
  atac_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = paste0(path, '/outs/fragments.tsv.gz'),
    min.cells = 1
  )
  
  if (verbose)
    message(" % creating seurat object...")
  seurat <- CreateSeuratObject(counts = atac_assay, assay = assay, meta.data = meta, project = project)
  
  if (!is.null(annotations)){
    if (verbose)
      message(" % adding provided annotation...")
    Annotation(seurat) <- annotations
  } else if (!is.null(annot_ens_version)){
    if (verbose)
      message(" % adding annotation...")
    ensdb <- retrieve_ensbd(version = annot_ens_version, force = force_ens_dl)
    annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
    genome(annotations) <- "hg38"
    Annotation(seurat) <- annotations
  }
  
  if (verbose)
    message("done.")
  return(seurat)
}

# convert data to h5ad format
data_to_h5ad.Seurat <- function(object,
                                assay = DefaultAssay(object),
                                savefile = NULL,
                                verbose = F)
{
  if (verbose)
    message("start to create the anndata object...")
  adata <- anndata::AnnData(X = t(object[[assay]]@data),
                            obs = object@meta.data,
                            var = object[[assay]]@meta.features,
                            layers = list(count = t(object[[assay]]@counts)),
                            obsm = setNames(lapply(names(object@reductions), function(x) Embeddings(object, x)), names(object@reductions))
                            )
  if (verbose)
    message("done.")
  
  if (!is.null(savefile)){
    if (verbose)
      message(paste0("saving the anndata object to file: ",savefile,"..."))
    adata$write_h5ad(savefile)
    if (verbose)
      message("done.")
  }
  
  return(adata)
}

data_to_h5ad.default <- function(object,
                                 vars = NULL,
                                 obs = NULL,
                                 obsm = list(),
                                 layers = list(),
                                 savefile = NULL,
                                 verbose=F)
{
  if (verbose)
    message("start to create the anndata object...")
  if (is.null(vars))
    vars <- data.frame(id = colnames(object), row.names = colnames(object))
  if (is.null(obs))
    obs <- data.frame(id = rownames(object), row.names = rownames(object))
  adata <- anndata::AnnData(X = object,
                            obs = obs,
                            var = vars,
                            obsm = obsm,
                            layers = layers
                            )
  if (verbose)
    message("done.")
  
  if (!is.null(savefile)){
    if (verbose)
      message(paste0("saving the anndata object to file: ",savefile,"..."))
    adata$write_h5ad(savefile)
    if (verbose)
      message("done.")
  }
  
  return(adata)
}

data_to_h5ad <- function(object, ...) {
  UseMethod(generic = 'data_to_h5ad', object = object)
}

# summarize data to the group level (by average / sum / majority)
summarize_data_to_groups.Matrix <- function(object, groups, group_cols = F, use_mean = T)
{
  groups <- factor(groups)
  if (group_cols){
    if (use_mean){
      summfunc <- Matrix::rowMeans
    } else{
      summfunc <- Matrix::rowSums
    }
  } else{
    if (use_mean){
      summfunc <- Matrix::colMeans
    } else{
      summfunc <- Matrix::colSums
    }
  }
  
  summ_mat <- sapply(levels(groups), function(x){
    if (group_cols) obj <- object[,which(groups==x)]
    if (!group_cols) obj <- object[which(groups==x),]
    summfunc(obj)
  })
  if (!group_cols) summ_mat <- t(summ_mat)
  summ_mat <- Matrix(summ_mat, sparse = T)
  return(summ_mat)
}

summarize_data_to_groups.data.frame <- function(object, groups)
{
  groups <- factor(groups)
  ident_mat_groups <- sparseMatrix(i = which(!is.na(groups)),
                                   j = as.numeric(groups)[!is.na(groups)],
                                   x = 1,
                                   dims = c(length(groups), length(levels(groups))), dimnames = list(names(groups), levels(groups)))
  ident_mat_groups_norm <- colNorm(ident_mat_groups)
  
  df_summ <- do.call(cbind.data.frame, lapply(1:ncol(object), function(i){
    dat <- object[,i]
    if (is.numeric(dat)){
      if (sum(! is.na(dat)) == 0)
        return(setNames(rep(NA, length(levels(groups))), levels(groups)))
      dat <- matrix(dat, nrow = 1, dimnames = list(colnames(object)[i], rownames(object)))
      dat[is.na(dat)] <- 0
      summ_dat <- dat %*% ident_mat_groups_norm
      return(setNames(as.numeric(summ_dat), colnames(ident_mat_groups_norm)))
    } else{
      if (sum(! is.na(dat))==0){
        if (is.null(levels(dat)))
          return(setNames(rep(NA, length(levels(groups))), levels(groups)))
        return(setNames(factor(rep(NA, length(levels(groups))), levels = levels(dat)), levels(groups)))
      }
      dat <- factor(dat)
      lab_mat_groups <- sparseMatrix(i = which(!is.na(dat)),
                                     j = as.numeric(dat)[which(!is.na(dat))],
                                     x = 1,
                                     dims = c(length(dat), length(levels(dat))), dimnames = list(rownames(object), levels(dat)))
      summ_mat <- t(lab_mat_groups) %*% ident_mat_groups
      summ_dat <- factor(setNames(rownames(summ_mat)[apply(summ_mat, 2, which.max)], levels(groups)), levels = levels(dat))
      return(summ_dat)
    }
  }))
  colnames(df_summ) <- colnames(object)
  
  return(df_summ)
}

summarize_data_to_groups <- function(object, ...) {
  UseMethod(generic = 'summarize_data_to_groups', object = object)
}

# do feature plot with base R
plotFeature <- function (coord,
                         values = NULL,
                         emphasize = NULL,
                         col = NULL,
                         colorPal = NULL,
                         col_alpha = NULL,
                         balance_col = TRUE,
                         mask_col = "#efefef50",
                         value_ceiling = NULL,
                         value_floor = NULL,
                         pt_border = FALSE,
                         col_border = "#303030",
                         lwd_border = 0.2,
                         main = NA,
                         axis = FALSE,
                         xlab = "Dim-1", 
                         ylab = "Dim-2",
                         cex = 1,
                         cex.main = 1,
                         cex.lab = 1,
                         cex.axis = 1,
                         random_order = TRUE,
                         sort_by_value = FALSE,
                         edges = NULL,
                         lwd.edges = 0.2,
                         col.edges = "#bdbdbd50",
                         do_label = FALSE,
                         cex.label = 1,
                         label_min_cell = 10,
                         label_round = FALSE,
                         label_round_cex = cex.label * 3,
                         label_round_col = "#303030",
                         label_round_lwd = 0.5,
                         label_round_bg = "#fefefe50",
                         do_legend = F,
                         legend_pos = "right",
                         legend_cex = 1,
                         legend_pt_cex = legend_cex * 2,
                         ...) 
{
  if (! is.null(col)){ do_legend = F; do_label = F }
  if (is.numeric(values)){ do_legend = F; do_label = F }
  
  if (is.null(col)) {
    if (is.numeric(values)) {
      values <- c(value_ceiling, value_floor, values)
      
      if (is.null(colorPal)) {
        colorPal <- grDevices::colorRampPalette(c("darkgreen", 
                                                  "yellow", "red"))
      } else if (is.character(colorPal)) {
        colorPal <- colorRampPalette(colorPal)
      }
      
      if (balance_col & min(values, na.rm=T) < 0 & max(values, na.rm=T) > 0){
        values <- c(-max(abs(values)), max(abs(values)), values)
        cellColor <- adjustcolor(colorPal(30), alpha = 0.8)[as.numeric(cut(values, 
                                                                           breaks = 30, right = F, include.lowest = T))]
        values <- values[-(1:2)]
        cellColor <- cellColor[-(1:2)]
      } else{
        cellColor <- adjustcolor(colorPal(30), alpha = 0.8)[as.numeric(cut(values, 
                                                                           breaks = 30, right = F, include.lowest = T))]
        if (min(values, na.rm = T) == 0) 
          cellColor[values == 0] <- "#bdbdbd30"
      }
      
      values <- values[(length(c(value_ceiling, value_floor))+1):length(values)]
      cellColor <- cellColor[(length(c(value_ceiling, value_floor))+1):length(cellColor)]
    }
    else {
      if (is.character(values)) values <- as.factor(values)
      
      cols <- NULL
      if (!is.null(names(colorPal)) & sum(values %in% names(colorPal))>0){
        cols <- colorPal
      } else{
        if (is.null(colorPal)) colorPal <- scales::hue_pal()
        if (is.character(colorPal)) colorPal <- colorRampPalette(colorPal)
        if (is.function(colorPal)) cols <- colorPal(length(levels(values)))
        if (is.null(names(cols))) cols <- setNames(cols, levels(values))
      }
      cellColor <- cols[as.character(values)]
    }
    
    if (!is.null(col_alpha))
      cellColor <- adjustcolor(cellColor, col_alpha)
  }
  else {
    if (length(col) == 1) 
      col <- rep(col, nrow(coord))
    cellColor <- col
  }
  col_border <- rep(col_border, length(cellColor))
  col_border[is.na(cellColor)] <- NA
  
  if (axis){
    plot(coord, type = "n", main = main, xlab = xlab, ylab = ylab, 
         cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, 
         ...)
  } else{
    plot(coord, type = "n", main = main, bty = "n", xlab = NA, ylab = NA, xaxt = "n", yaxt = "n",
         cex.main = cex.main, ...)
  }
  if (!is.null(edges) & (is.matrix(edges) | is.data.frame(edges))) {
    for (i in 1:nrow(edges)) lines(coord[as.numeric(edges[i, ]), 1], coord[as.numeric(edges[i, ]), 2],
                                   lwd = lwd.edges, 
                                   col = col.edges)
  }
  
  if (sum(is.na(cellColor))>0){
    idx <- which(!is.na(cellColor))
    if (is.null(emphasize))
      emphasize <- idx
    emphasize <- intersect(emphasize, idx)
  }
  
  if (is.null(emphasize)) {
    idx_order <- seq(nrow(coord))
    if (random_order){
      idx_order <- sample(idx_order)
    } else if (sort_by_value){
      idx_order <- order(values)
    }
    
    if (pt_border){
      points(coord[idx_order,], col = col_border[idx_order], bg = cellColor[idx_order], pch = 21, cex = cex, lwd = lwd_border)
    } else{
      points(coord[idx_order,], col = cellColor[idx_order], pch = 16, cex = cex)
    }
  }
  else {
    points(coord, col = mask_col, pch = 16, cex = cex)
    idx_order <- emphasize
    if (length(idx_order) > 0){
      if (random_order){
        idx_order <- sample(idx_order)
      } else if (sort_by_value){
        idx_order <- idx_order[order(values[idx_order])]
      }
      
      if (pt_border){
        points(coord[idx_order, ], col = col_border[idx_order], bg = cellColor[idx_order], pch = 21, cex = cex, lwd = lwd_border)
      } else{
        points(coord[idx_order, ], col = cellColor[idx_order], pch = 16, cex = cex)
      }
    }
  }
  
  if (is.factor(values) & do_label){
    labs <- table(values)
    labs <- names(which(labs >= label_min_cell))
    labs_coord <- t(sapply(labs, function(x) colMeans(as.data.frame(coord)[which(values == x),], na.rm = T)))
    if (label_round) points(labs_coord, pch = 21, cex = label_round_cex, lwd = label_round_lwd, col = label_round_col, bg = label_round_bg)
    text(labs_coord, labels = labs, cex = cex.label)
  }
  if (do_legend){
    legend(legend_pos, legend = names(cols), pch = 16, col = cols, pt.cex = legend_pt_cex, cex = legend_cex, bty="n")
  }
}

