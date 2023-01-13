# this script is to apply the RNA-ATAC integration pipeline to each of the RNA-ATAC paired samples

library(Matrix)
library(dplyr)
library(Seurat)
library(Signac)
library(reticulate)

source("ext/util.r")

# read paired sample list
paired.sample.list <- read.table("data/bimodal/paired_samples.csv",
                               sep=",", head=T, stringsAsFactors = F)

# ATAC
seurat.atac.list <- setNames(lapply(paired.sample.list[,1], function(x){
  seurat.atac <- create_atac_seurat(paste0("data/scATACseq/processed/",x), sample = x, annot_ens_version = 93)
  seurat.atac <- NucleosomeSignal(object = seurat.atac)
  seurat.atac <- TSSEnrichment(object = seurat.atac, fast = FALSE)
  seurat.atac$high.tss <- ifelse(seurat.atac$TSS.enrichment > 2, 'High', 'Low')
  seurat.atac$pct_reads_in_peaks <- seurat.atac$peak_region_fragments / seurat.atac$passed_filters * 100
  seurat.atac$blacklist_ratio <- seurat.atac$blacklist_region_fragments / seurat.atac$peak_region_fragments

  seurat.atac <- subset(seurat.atac,
                        subset = peak_region_fragments > 200 &
                                 peak_region_fragments < 10000 &
                                 pct_reads_in_peaks > 20 &
                                 blacklist_ratio < 0.025 &
                                 nucleosome_signal < 3 &
                                 TSS.enrichment > 2) %>%
    RunTFIDF()
  det_rates <- rowMeans(seurat.atac@assays$peaks@counts > 0)
  VariableFeatures(seurat.atac) <- names(which(det_rates>0.005))
  seurat.atac <- RunSVD(seurat.atac) %>%
    RunUMAP(reduction="lsi", dims = 2:30)
    FindNeighbors(reduction="lsi", dims=2:30) %>%
    FindClusters(resolution = 0.5)
  seurat.atac[['RNA']] <- CreateAssayObject(GeneActivity(seurat.atac))
  seurat.atac <- NormalizeData(seurat.atac, assay = "RNA")
  return(seurat.atac)
}), paired.sample.list[,1])

saveRDS(seurat.atac.list, file="data/scATACseq/seurat_objs/samples_list.rds")

# RNA
blacklist <- read.table("ext/blacklist.txt")[,1]
seurat.rna.list <- SplitObject(readRDS("data/scRNAseq/seurat_objs/merged.rds"), "sample")

# calculate the correlations between RNA and ATAC clusters for each paired RNA-ATAC data
source_python("ext/matching.py")
grp.idx.mat.by.pair <- apply(paired.sample.list, 1, function(x){
  rna.sample <- x[2]
  atac.sample <- x[1]
  seurat.rna <- seurat.rna.list[[rna.sample]]
  seurat.atac <- seurat.atac.list[[atac.sample]]
  
  levels(seurat.rna$RNA_snn_res.0.5) <- paste0("RNA_", levels(seurat.rna$RNA_snn_res.0.5))
  levels(seurat.atac$peaks_snn_res.0.5) <- paste0("ATAC_", levels(seurat.atac$peaks_snn_res.0.5))
  
  ## correlate RNA and ATAC clusters (expression vs. gene activity)
  selected.hvg <- intersect(VariableFeatures(seurat.rna), rownames(seurat.atac[['RNA']]))
  rna.hvg.expr <- AverageExpression(seurat.rna, group.by = "RNA_snn_res.0.5")$RNA[selected.hvg,]
  atac.hvg.expr <- AverageExpression(seurat.atac, group.by = "peaks_snn_res.0.5")$RNA[selected.hvg,]
  scc.cl2cl <- cor(atac.hvg.expr, rna.hvg.expr, method="spearman")
  # get the most similar RNA cluster for each ATAC cluster
  max.rna.cl <- colnames(scc.cl2cl)[apply(scc.cl2cl, 1, which.max)]
  names(max.rna.cl) <- rownames(scc.cl2cl)
  # get the most similar ATAC cluster for each RNA cluster
  max.atac.cl <- rownames(scc.cl2cl)[apply(scc.cl2cl, 2, which.max)]
  names(max.atac.cl) <- colnames(scc.cl2cl)
  
  ## pairs of matched RNA and ATAC clusters
  pairs <- unique(rbind(cbind(names(max.rna.cl), max.rna.cl),cbind(max.atac.cl, names(max.atac.cl))))
  graph <- graph_from_edgelist(pairs, directed=FALSE)
  grouped.cells <- components(graph)
  group.vec <- grouped.cells$membership
  group.list <- lapply(sort(unique(group.vec)), function(i){
    names(group.vec)[which(group.vec==i)]
  })
  
  ## CCA integration
  DefaultAssay(seurat.atac) <- "RNA"
  VariableFeatures(seurat.atac) <- selected.hvg
  seurat.atac <- ScaleData(seurat.atac)
  anchor_features <- SelectIntegrationFeatures(
    object.list = list('RNA'=seurat.rna, 'ATAC'=seurat.atac),
    nfeatures = length(selected.hvg),
    assay = c("RNA","RNA")
  )
  
  seurat.rna <- RenameCells(seurat.rna, new.names=paste(colnames(seurat.rna), "RNA", sep="_"))
  seurat.atac <- RenameCells(seurat.atac, new.names=paste(colnames(seurat.atac), "ATAC", sep="_"))
  seurat.cca <- RunCCA(object1=seurat.rna, object2=seurat.atac, assay1="RNA", assay2="RNA", features = anchor_features, num.cc = 20)
  seurat.cca$modality <- sapply(colnames(seurat.cca), function(x){
    vec <- strsplit(x, "_")[[1]]
    vec[length(vec)]
  })
  
  # run MCMF within each component
  rna.full.dat <- data.frame(Embeddings(seurat.cca, "cca")[which(seurat.cca$modality == "RNA"),])
  atac.full.dat <- data.frame(Embeddings(seurat.cca, "cca")[which(seurat.cca$modality == "ATAC"),])
  grp.idx.list <- lapply(group.list, function(cls){
    rna.cl <- grep("RNA", cls, value=T)
    atac.cl <- grep("ATAC", cls, value=T)
    rna.cells <- colnames(seurat.rna)[which(seurat.rna$RNA_snn_res.0.5%in%rna.cl)]
    atac.cells <- colnames(seurat.atac)[which(seurat.atac$peaks_snn_res.0.5%in%atac.cl)]
    rna.dat <- rna.full.dat[rna.cells,]
    atac.dat <- atac.full.dat[atac.cells,]
    cost_graph <- get_cost_knn_graph(
      source = rna.dat,
      target = atac.dat,
      knn_k = 10,
      knn_n_jobs = 10,
      null_cost_percentile = 99,
      capacity_method = "uniform"
    )
    grp.idx <- setNames(data.frame(do.call(cbind, mcmf(cost_graph))+1), c("RNA","ATAC"))
    grp.cell <- data.frame("RNA"=rownames(rna.dat)[grp.idx[,1]], "ATAC"=rownames(atac.dat)[grp.idx[,2]])
  })
  grp.idx.mat <- do.call(rbind, grp.idx.list)
  return(grp.idx.mat)
}

rna_atac_pairs <- setNames(data.frame(do.call(rbind, grp.idx.mat.by.pair),
                                      rep(paired.sample.list[,2], sapply(grp.idx.mat.by.pair, nrow)),
                                      rep(paired.sample.list[,1], sapply(grp.idx.mat.by.pair, nrow))),
                           c("RNA_cell_barcode", "ATAC_cell_barcode", "RNA_lib_idx", "ATAC_lib_idx"))

rna_atac_pairs$RNA_cell_ident <- paste0(sapply(strsplit(rna_atac_pairs$RNA_cell_barcode,"-"),"[",1), "_", rna_atac_pairs$RNA_lib_idx)
rna_atac_pairs$RNA_cell_ident[is.na(rna_atac_pairs$RNA_cell_barcode)] <- NA
rna_atac_pairs$ATAC_cell_ident <- paste0(sapply(strsplit(rna_atac_pairs$ATAC_cell_barcode,"-"),"[",1), "_", rna_atac_pairs$ATAC_lib_idx)
rna_atac_pairs$ATAC_cell_ident[is.na(rna_atac_pairs$ATAC_cell_barcode)] <- NA
rna_atac_pairs <- data.frame(dataset = "RNA_ATAC_paired",
                             rna_atac_pairs[which(!is.na(rna_atac_pairs$RNA_cell_ident) & !is.na(rna_atac_pairs$ATAC_cell_ident) & rna_atac_pairs$RNA_cell_ident %in% seurat_rna_all$cell_ident),c(3,5,4,6)])

saveRDS(rna_atac_pairs, file="data/bimodal/paired_cells_rna_atac_paired.rds")
