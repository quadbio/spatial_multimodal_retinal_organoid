library(Seurat)
library(Signac)
library(dplyr)
library(reticulate)
source("ext/util.r")
source_python("ext/matching.py")

# benchmark with multiome: matching of cells
blacklist <- read.table("ext/genes_blacklist.txt")[,1]
srts <- readRDS("data/multiome/seurat_objs/sample_QCed_list.rds")

## clustering with mono-/bi-modal information
srts <- lapply(srts, function(srt){
  # RNA
  DefaultAssay(srt) <- "RNA"
  srt <- NormalizeData(srt) %>%
    FindVariableFeatures(nfeatures = 3000)
  VariableFeatures(srt) <- setdiff(VariableFeatures(srt), blacklist)
  srt <- CellCycleScoring(srt, g2m.features=cc.genes.updated.2019$g2m.genes, s.features=cc.genes.updated.2019$s.genes) %>%
    ScaleData(vars.to.regress = c("G2M.Score", "S.Score")) %>%
    RunPCA(npcs = 20, verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.5)
  
  # ATAC
  DefaultAssay(srt) <- "ATAC"
  srt <- RunTFIDF(srt) %>%
    FindTopFeatures(min.cutoff = ncol(srt) * 0.005)
  det_rates <- rowMeans(srt[['ATAC']]@counts > 0)
  VariableFeatures(srt) <- names(which(det_rates>0.005))
  srt <- RunSVD(srt) %>%
    FindNeighbors(reduction = "lsi", dims = 2:30) %>%
    FindClusters(resolution = 0.5)
  srt[['ACT']] <- CreateAssayObject(GeneActivity(srt))
  
  # bimodal
  srt <- FindMultiModalNeighbors(srt,
                                 reduction.list = list("pca","lsi"),
                                 dims.list = list(1:20, 2:30),
                                 verbose = T) %>%
    FindClusters(graph.name = "wsnn", resolution = 0.5)
  return(srt)
})

## running the integration pipeline
linked_idx_multiome <- lapply(srts, function(srt){
  ## RNA assay
  DefaultAssay(srt) <- "RNA"
  rna <- DietSeurat(srt, assays = "RNA")
  rna$modality <- "RNA"
  
  ## ATAC assay
  DefaultAssay(srt) <- "ATAC"
  atac <- DietSeurat(srt, assays = "ATAC")
  atac[['RNA']] <- srt[['ACT']]
  atac$modality <- "ATAC"
  
  ## integration: CCA
  rna <- RenameCells(rna, new.names = paste(colnames(rna), "RNA", sep="_"))
  atac <- RenameCells(atac, new.names = paste(colnames(atac), "ATAC", sep="_"))
  anchor_features <- SelectIntegrationFeatures(
    object.list = list('RNA'=rna, 'ATAC'=atac),
    nfeatures = 2000,
    assay = rep('RNA', 2)
  )
  cca <- RunCCA(object1=rna, object2=atac, assay1="RNA", assay2="RNA", features = anchor_features, num.cc = 20)
  
  ## anchor RNA and ATAC clusters based on gene expression vs. gene activities
  selected.hvg <- intersect(VariableFeatures(rna), rownames(atac[['RNA']]))
  rna.expr <- AverageExpression(rna, features = selected.hvg, assays="RNA")$RNA
  colnames(rna.expr) <- paste0("RNA_", colnames(rna.expr))
  atac.expr <- AverageExpression(atac, features = selected.hvg, assays="RNA")$RNA
  colnames(atac.expr) <- paste0("ATAC_", colnames(atac.expr))
  scc.cl2cl <- cor(atac.expr, rna.expr, method="spearman")
  max.rna.cl <- colnames(scc.cl2cl)[apply(scc.cl2cl, 1, which.max)]
  max.atac.cl <- rownames(scc.cl2cl)[apply(scc.cl2cl, 2, which.max)]
  names(max.rna.cl) <- rownames(scc.cl2cl)
  names(max.atac.cl) <- colnames(scc.cl2cl)
  pairs <- unique(rbind(cbind(names(max.rna.cl), max.rna.cl),cbind(max.atac.cl, names(max.atac.cl))))
  graph <- graph_from_edgelist(pairs, directed=FALSE)
  grouped.cells <- components(graph)
  group.vec <- grouped.cells$membership
  group.list <- lapply(sort(unique(group.vec)), function(i){
    names(group.vec)[which(group.vec==i)]
  })
  
  ## MCMF on anchored clusters
  rna.full.dat <- data.frame(Embeddings(cca, "cca")[which(cca$modality == "RNA"),])
  atac.full.dat <- data.frame(Embeddings(cca, "cca")[which(cca$modality == "ATAC"),])
  grp.idx.mat <- Reduce(rbind, lapply(group.list, function(x){
      rna.cl <- sub("RNA_", "", grep("RNA", x, value=T))
      atac.cl <- sub("ATAC_", "", grep("ATAC", x, value=T))
      rna.cells <- colnames(rna)[which(rna$RNA_snn_res.0.5%in%rna.cl)]
      atac.cells <- colnames(atac)[which(atac$ATAC_snn_res.0.5%in%atac.cl)]
      rna.dat <- rna.full.dat[rna.cells,]
      atac.dat <- atac.full.dat[atac.cells,]
      cost_graph <- get_cost_knn_graph(source = rna.dat,
                                       target = atac.dat,
                                       knn_k = 10,
                                       knn_n_jobs = 10,
                                       null_cost_percentile = 99,
                                       capacity_method = "uniform")
      grp.idx <- setNames(data.frame(do.call(cbind, mcmf(cost_graph))+1), c("RNA","ATAC"))
      grp.cell <- data.frame(RNA = rownames(rna.dat)[grp.idx[,1]],
                             ATAC = rownames(atac.dat)[grp.idx[,2]])
      return(grp.cell)
  }))
  
  return(grp.idx.mat)
})

## compare the cluster label of a cell with its paired cell
### check proportions of paired cells being in the same clusters
linked_idx_multiome <- setNames(lapply(names(linked_idx_multiome), function(sample){
  linked_idx <- linked_idx_multiome[[sample]]
  srt <- srts[[sample]]
  
  linked_idx$RNA <- gsub("_RNA$","", linked_idx$RNA)
  linked_idx$ATAC <- gsub("_ATAC$","", linked_idx$ATAC)
  linked_idx$same_RNA_cl <- srt$RNA_snn_res.0.5[linked_idx$RNA] == srt$RNA_snn_res.0.5[linked_idx$ATAC]
  linked_idx$same_ATAC_cl <- srt$ATAC_snn_res.0.5[linked_idx$RNA] == srt$ATAC_snn_res.0.5[linked_idx$ATAC]
  linked_idx$same_wnn_cl <- srt$wsnn_res.0.5[linked_idx$RNA] == srt$wsnn_res.0.5[linked_idx$ATAC]
  
  linked_idx_permut <- lapply(1:100, function(i){
    set.seed(i)
    linked_idx$ATAC <- sample(linked_idx$ATAC)
    linked_idx$same_RNA_cl <- srt$RNA_snn_res.0.5[linked_idx$RNA] == srt$RNA_snn_res.0.5[linked_idx$ATAC]
    linked_idx$same_ATAC_cl <- srt$ATAC_snn_res.0.5[linked_idx$RNA] == srt$ATAC_snn_res.0.5[linked_idx$ATAC]
    linked_idx$same_wnn_cl <- srt$wsnn_res.0.5[linked_idx$RNA] == srt$wsnn_res.0.5[linked_idx$ATAC]
    return(linked_idx)
  })
  
  res <- list(obs = linked_idx,
              permut = linked_idx_permut)
  return(res)
}), names(linked_idx_multiome))

for(sample in names(linked_idx_multiome)){
  linked_idx <- linked_idx_multiome[[sample]]$obs
  linked_idx_permut <- linked_idx_multiome[[sample]]$permut
  
  layout(matrix(1:2,nrow=1)); par(mar=c(6,5,1,1))
  barplot(cbind(c(sum(linked_idx$same_RNA_cl+linked_idx$same_ATAC_cl==2, na.rm=T),sum(linked_idx$same_RNA_cl+linked_idx$same_ATAC_cl==1, na.rm=T),sum(!linked_idx$same_RNA_cl & !linked_idx$same_ATAC_cl, na.rm=T))/sum(!is.na(linked_idx$same_RNA_cl) & !is.na(linked_idx$same_ATAC_cl)),
                rowMeans(sapply(linked_idx_permut, function(linked_idx) c(sum(linked_idx$same_RNA_cl & linked_idx$same_ATAC_cl, na.rm=T),sum(linked_idx$same_RNA_cl+linked_idx$same_ATAC_cl==1, na.rm=T),sum(!linked_idx$same_RNA_cl & !linked_idx$same_ATAC_cl, na.rm=T))/sum(!is.na(linked_idx$same_RNA_cl) & !is.na(linked_idx$same_ATAC_cl)) ))),
          col = c("#303030","#909090","#dedede"), border=NA, names.arg = c("Real","Permutation"), las=2, ylab = "Proportion of pairs", cex.lab=1.5)
  barplot(cbind(c(sum(linked_idx$same_wnn_cl, na.rm=T),sum(!linked_idx$same_wnn_cl, na.rm=T))/sum(!is.na(linked_idx$same_wnn_cl)),
                rowMeans(sapply(linked_idx_permut, function(linked_idx) c(sum(linked_idx$same_wnn_cl, na.rm=T),sum(!linked_idx$same_wnn_cl, na.rm=T))/sum(!is.na(linked_idx$same_wnn_cl)) ))),
          col = c("#303030","#dedede"), border=NA, names.arg = c("Real","Permutation"), las=2, ylab = "Proportion of pairs", cex.lab=1.5)

}

### calculate jaccard index
jaccard_mats <- lapply(names(linked_idx_multiome), function(sample){
  srt <- srts[[sample]]
  linked_idx <- linked_idx_multiome[[sample]]$obs
  
  sapply(levels(srt$RNA_snn_res.0.5), function(cl)
    sapply(levels(srt$RNA_snn_res.0.5), function(cl2)
      length(intersect(which(srt$RNA_snn_res.0.5[linked_idx$RNA]==cl),
                       which(srt$RNA_snn_res.0.5[linked_idx$ATAC]==cl2))) /
      length(intersect(which(!is.na(linked_idx$RNA) & !is.na(linked_idx$ATAC)),
                       union(which(srt$RNA_snn_res.0.5[linked_idx$RNA]==cl),
                             which(srt$RNA_snn_res.0.5[linked_idx$ATAC]==cl2))))))
})

for(sample in names(linked_idx_multiome)){
  avg_expr_cl <- sapply(levels(srts[[sample]]$RNA_snn_res.0.5), function(cl)
                   rowMeans(srts[[sample]][['RNA']]@data[,which(srts[[sample]]$RNA_snn_res.0.5 == cl)]))
  hcl_cl <- hclust(as.dist(1 - cor(avg_expr_cl[VariableFeatures(srts[[sample]], assay="RNA"),])))

  gplots::heatmap.2(jaccard_mats[[sample]], Rowv=as.dendrogram(hcl_cl), Colv=as.dendrogram(hcl_cl), scale="none", trace="none", col=greyscale_colscheme(30), key=T)
}
#### compare jaccard index of w15 and w36
boxplot(diag(jaccard_mats[[1]]), diag(jaccard_mats[[2]]), frame=F, names=c("Week-15", "Week-36"), cex.names = 1.2, ylab="Jaccard index", cex.lab=1.2, ylim=c(0,1))






# benchmark by comparing multiome ATAC and RNA-paired ATAC profiles in different cell types
seurat <- readRDS("data/bimodal/onlyRetina_integrated.rds")
seurat_with_atac <- subset(seurat, subset = paired_ATAC)

cols_major_ct <- c("RPC" = "#85C1E9",
                   "MG" = "#2980B9",
                   "RPE" = "#2C3E50",
                   "Rods" = "#CB4335",
                   "Cones" = "#E67E22",
                   "BC" = "#D4AC0D",
                   "AC" = "#1D8348",
                   "HC" = "#17A589",
                   "RGC" = "#117A65")
avg_acc_ct <- setNames(lapply(sort(unique(seurat_with_atac$orig.ident)), function(ds) sapply(levels(seurat$major_ct), function(ct)
  rowMeans(seurat_with_atac@assays$ATAC@data[,which(seurat_with_atac$orig.ident==ds & seurat_with_atac$major_ct == ct)]))), sort(unique(seurat_with_atac$orig.ident)))
corr_acc_ct <- cor(do.call(cbind, avg_acc_ct))
hcl_acc_ct <- hclust(as.dist(1-corr_acc_ct), method="ward.D2")
corr_acc_ct_mods <- cor(avg_acc_ct[[2]], avg_acc_ct[[1]])

gplots::heatmap.2(corr_acc_ct, Rowv=as.dendrogram(hcl_acc_ct), Colv=as.dendrogram(hcl_acc_ct), scale="none", trace="none", key=F, keysize=0.5, ColSideColors = rep(cols_major_ct[levels(seurat$major_ct)],2), RowSideColors = rep(c("#303030","#cdcdcd"),each=length(levels(seurat$major_ct))), col=bluered_colscheme(30), margins = c(7,7))
gplots::heatmap.2(corr_acc_ct_mods, Rowv=NA, Colv=NA, dendrogram="none", scale="col", trace="none", key=F, keysize=0.5, ColSideColors = cols_major_ct[levels(seurat$major_ct)], RowSideColors = cols_major_ct[levels(seurat$major_ct)], col=bluered_colscheme(30), margins = c(8,8))





# benchmark by estimating similalrities between RNA and paired ATAC gene activity scores
seurat <- readRDS("data/bimodal/onlyRetina_integrated.rds")
seurat_ATAC <- readRDS("data/scATACseq/seurat_objs/onlyRetina_incl_multiome.rds")
seurat_ATAC[['RNA']] <- CreateAssayObject(GeneActivity(seurat_ATAC))

## RNA-ATAC pairs
### pairs: paired RNA-ATAC + multiome
rna_atac_pairs <- readRDS(file="data/bimodal/paired_cells_rna_atac.rds")
rna_atac_pairs <- rna_atac_pairs[which(rna_atac_pairs$RNA_cell_ident %in% seurat$cell_ident[seurat$paired_ATAC] & rna_atac_pairs$ATAC_cell_ident %in% seurat_ATAC$cell_ident),]
rna_atac_pairs$RNA_cell_idx <- pbapply::pbsapply(rna_atac_pairs$RNA_cell_ident, function(cell) which(seurat$cell_ident == cell))
rna_atac_pairs$ATAC_cell_idx <- pbapply::pbsapply(rna_atac_pairs$ATAC_cell_ident, function(cell) which(seurat_ATAC$cell_ident == cell))

### pairs: multiome by constrained CCA-MCMF
linked_idx_obs_multiome <- lapply(names(linked_idx_multiome), function(sample){
  linked_idx <- linked_idx_multiome[[sample]]$obs %>%
    filter(!is.na(RNA) & !is.na(ATAC)) %>%
    mutate(dataset = "multiome_MCMF",
           RNA_lib_idx = sample,
           ATAC_lib_idx = sample) %>%
    mutate(RNA_cell_ident = paste0(gsub("\\-1$","",RNA),"_",sample)) %>%
    mutate(ATAC_cell_ident = paste0(gsub("\\-1$","",ATAC),"_",sample)) %>%
    select(dataset, RNA_lib_idx, RNA_cell_ident, ATAC_lib_idx, ATAC_cell_ident) %>%
    filter(RNA_cell_ident %in% seurat$cell_ident & ATAC_cell_ident %in% seurat_ATAC$cell_ident)
})
linked_idx_multiome_MCMF <- Reduce(rbind, linked_idx_obs_multiome) %>%
  mutate(RNA_cell_idx = setNames(1:ncol(seurat), seurat$cell_ident)[RNA_cell_ident]) %>%
  mutate(ATAC_cell_idx = setNames(1:ncol(seurat_ATAC), seurat_ATAC$cell_ident)[ATAC_cell_ident]) %>%
  mutate(highres_cluster = setNames(seurat$highres_cl_lab, seurat$cell_ident)[RNA_cell_ident])

### pairs: paired RNA-ATAC + multiome + multiome by constrained CCA-MCMF
rna_atac_pairs_add <- rbind(rna_atac_pairs[,colnames(linked_idx_multiome_MCMF)],
                            linked_idx_multiome_MCMF)

## correlation across genes per pair, on highres cluster level
groups2check <- rna_atac_pairs_add %>%
  group_by(dataset, RNA_lib_idx, highres_cluster) %>%
  summarise(num_cells = n_distinct(RNA_cell_ident)) %>%
  filter(num_cells > 20)
groups2check$cell_idx <- apply(groups2check, 1, function(x) which(rna_atac_pairs_add$RNA_lib_idx == x[2] & rna_atac_pairs_add$highres_cluster == x[3]) )

expr_rna <- as.matrix(seurat[['RNA']]@data[highvar,rna_atac_pairs_add$RNA_cell_idx])
expr_geneactivity <- as.matrix(seurat_ATAC[['RNA']]@data[highvar,rna_atac_pairs_add$ATAC_cell_idx])

corr_rna2geneact_highres <- pbapply::pbsapply(groups2check$cell_idx, function(idx){
  rna <- rowMeans(expr_rna[,idx])
  act <- rowMeans(expr_geneactivity[,idx])
  cor(rna, act)
})
groups2check$corr_rna2geneact <- corr_rna2geneact_highres

corr_rna2geneact_highres_permut <- pbapply::pblapply(1:nrow(groups2check), function(i){
  idx1 <- groups2check$cell_idx[[i]]
  rna <- rowMeans(expr_rna[,idx1])
  corrs <- sapply(1:50, function(n){
    set.seed(n)
    idx2 <- sample(which(rna_atac_pairs_add$RNA_lib_idx == groups2check$RNA_lib_idx[i]), length(idx1))
    act <- rowMeans(expr_geneactivity[,idx2])
    cor(rna, act)
  })
  return(corrs)
})
groups2check$corr_rna2geneact_permut <- corr_rna2geneact_highres_permut
groups2check$corr_rna2geneact_permut_mean <- sapply(corr_rna2geneact_highres_permut, mean)
groups2check$corr_rna2geneact_over_permut <- groups2check$corr_rna2geneact - groups2check$corr_rna2geneact_permut_mean

#### correlation in relative to permutations
sample_ages <- seurat@meta.data[,c("orig.ident","sample","age_weeks")] %>%
  distinct() %>%
  filter(orig.ident != "hRO_Cowan") %>%
  select(sample, age_weeks) %>%
  tibble::deframe()
sample_ages <- c(sample_ages, setNames(sample_ages[c("w15_retina_org", "w36_retina_org")], c("w15_retina_org_MCMF","w36_retina_org_MCMF")))
sample_datasets <- seurat@meta.data[,c("orig.ident","sample")] %>%
  distinct() %>%
  filter(orig.ident != "hRO_Cowan") %>%
  select(sample, orig.ident) %>%
  tibble::deframe()
sample_datasets <- c(sample_datasets, setNames(rep("hRO_timecourse_multiome_MCMF",2), c("w15_retina_org_MCMF","w36_retina_org_MCMF")))
ordered_samples <- sort(unique(groups2check$RNA_lib_idx))[order(sample_datasets[sort(unique(groups2check$RNA_lib_idx))], sample_ages[sort(unique(groups2check$RNA_lib_idx))])]

par(mar=c(8,5,1,1), cex=0.8)
boxplot(c(tapply(groups2check$corr_rna2geneact_over_permut, groups2check$RNA_lib_idx, list)[ordered_samples],list(NA),list(NA),
          tapply(groups2check$corr_rna2geneact_over_permut, groups2check$dataset, list)),
        frame = F, xlab = "Sample age (weeks)", ylab = "Higher correlation than permut. (genes per state)", cex.lab = 1.2,
	las=2, names = c(sample_ages[ordered_samples],"","","Multiome","Multiome-MCMF","Paired"),
        col = c(rep(c("#303030","#909090","#cdcdcd"), c(2, 2, length(ordered_samples)-4)),NA,NA,"#303030","#909090","#cdcdcd"), border = "#000000",
        outline = F, pch = 16, cex = 1)
abline(h = 0, lty=2)

