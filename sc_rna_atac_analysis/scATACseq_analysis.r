library(Matrix)
library(Seurat)
library(Signac)
library(simspec)
library(harmony)
library(dplyr)

source("util.r")

# read data
seurat_rna_all <- readRDS("data/scRNAseq/seurat_objs/all_processed.rds")
seurat_rna_noMC <- readRDS("data/scRNAseq/seurat_objs/noMC_processed.rds")
seurat_rna_retina <- readRDS("data/bimodal/onlyRetina_integrated.rds")

seurat_rna_all$cell_ident <- paste0(sapply(strsplit(colnames(seurat_rna_all),"-"),"[",1), "_", seurat_rna_all$sample)
seurat_rna_all@meta.data$ct <- NA
seurat_rna_all$ct[colnames(seurat_rna_retina)] <- as.character(seurat_rna_retina$annot_ct_refine)
seurat_rna_all$ct[which(seurat_rna_all$RNA_css_snn_res.0.1 == "6")] <- "Mesenchyme"
seurat_rna_all$ct[is.na(seurat_rna_all$ct)] <- seurat_rna_noMC$annot_ct[colnames(seurat_rna_all)[which(is.na(seurat_rna_all$ct))]]
seurat_rna_all$ct <- factor(seurat_rna_all$ct, levels=c(levels(seurat_rna_retina$annot_ct_refine), c("Ast","Dien.","Tel.","Mesenchyme")))

seurat_multiome_atac <- readRDS("data/multiome/seurat_objs/atac_only.rds")
seurat_atac_integrated <- readRDS("data/scATACseq/seurat_objs/aggr_requant.rds")

# make a merged scATAC-seq seurat object including all the ATAC-seq cells
seurat_atac <- merge(seurat_atac_integrated, seurat_multiome_atac)

# apply ATAC-only analysis
seurat_atac <- RunTFIDF(seurat_atac)
det_rates <- rowMeans(seurat_atac@assays$peaks@counts > 0)
VariableFeatures(seurat_atac) <- names(which(det_rates>0.005))
seurat_atac <- RunSVD(seurat_atac) %>%
  RunUMAP(reduction="lsi", dims = 2:30)

### CSS integration
seurat_atac <- cluster_sim_spectrum(seurat_atac, label_tag = "sample", use_dr = "lsi", dims = 2:30, corr_method = "pearson", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 30, reduction.name = "css_pca", reduction.key = "CSSPCA_") %>%
  RunUMAP(reduction = "css_pca", dims = 1:20, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

saveRDS(seurat_atac, file="data/scATACseq/seurat_objs/all_incl_multiome.rds")


# transfer the RNA-based annotation
## update RNA-ATAC pairing by adding the multiome cells
rna_atac_pairs <- readRDS("data/bimodal/paired_cells_rna_atac_paired.rds")
multiome_pairs <- data.frame(dataset = "multiome",
                             RNA_lib_idx = seurat_rna_all$sample[seurat_rna_all$orig.ident == "hRO_timecourse_multiome"],
                             RNA_cell_ident = seurat_rna_all$cell_ident[seurat_rna_all$orig.ident == "hRO_timecourse_multiome"],
                             ATAC_lib_idx = seurat_rna_all$sample[seurat_rna_all$orig.ident == "hRO_timecourse_multiome"],
                             ATAC_cell_ident = seurat_rna_all$cell_ident[seurat_rna_all$orig.ident == "hRO_timecourse_multiome"], row.names = NULL)

rna_atac_pairs <- rbind(rna_atac_pairs, multiome_pairs)
rownames(rna_atac_pairs) <- NULL
saveRDS(rna_atac_pairs, file="data/bimodal/paired_cells_rna_atac.rds")

df_ct_atac <- setNames(data.frame(sapply(levels(seurat_rna_all$ct), function(ct){
  rna_cells <- seurat_rna_all$cell_ident[which(seurat_rna_all$ct==ct)]
  atac_cells <- rna_atac_pairs$ATAC_cell_ident[rna_atac_pairs$RNA_cell_ident %in% rna_cells]
  seurat_atac$cell_ident %in% atac_cells
}), row.names = colnames(seurat_atac)), levels(seurat_rna_all$ct))
seurat_atac$paired_ct <- "Unpaired"
seurat_atac$paired_ct[which(rowSums(df_ct_atac)==1)] <- colnames(df_ct_atac)[apply(df_ct_atac[which(rowSums(df_ct_atac)==1),],1,which)]
seurat_atac$paired_ct[which(rowSums(df_ct_atac)>1)] <- "Multiple"
seurat_atac$paired_ct <- factor(seurat_atac$paired_ct, levels=c(levels(seurat_rna_all$ct),"Multiple","Unpaired"))

# transferred-information-based cleanup: unpaired-enriched clusters
seurat_atac <- FindClusters(seurat_atac, graph.name = "peaks_lsi_snn", resolution = 0.1) %>%
  FindClusters(graph.name = "peaks_lsi_snn", resolution = 0.3) %>%
  FindClusters(graph.name = "peaks_lsi_snn", resolution = 0.5)

cl_freq <- table(seurat_atac$paired_ct, seurat_atac$peaks_lsi_snn_res.0.5)
cl_freq <- cl_freq[setdiff(rownames(cl_freq),c("Multiple")),]
p_enrichment <- apply(apply(cl_freq, 1, function(x){
  sapply(1:ncol(cl_freq), function(i){
    mat <- matrix(c(x[i],sum(x[-i]),sum(cl_freq[,i])-x[i], sum(cl_freq[,-i])-sum(x[-i])),c(2,2))
    p <- fisher.test(mat, alternative="g")$p.value
  })
}),2,p.adjust, method="bonferroni")
cl_unpaired <- colnames(cl_freq)[which(p_enrichment[,which(colnames(p_enrichment)=="Unpaired")]<0.01 & rowSums(p_enrichment[,which(colnames(p_enrichment)!="Unpaired")]<0.1)==0)]
cl_remained <- setdiff(levels(seurat_atac$peaks_lsi_snn_res.0.5), cl_unpaired)

seurat_atac_cleaned <- subset(seurat_atac, subset = peaks_lsi_snn_res.0.5 %in% cl_remained)

det_rates <- rowMeans(seurat_atac_cleaned@assays$peaks@counts > 0)
VariableFeatures(seurat_atac_cleaned) <- names(which(det_rates>0.005))
seurat_atac_cleaned <- RunSVD(seurat_atac_cleaned) %>%
  RunUMAP(reduction="lsi", dims = 2:30)

seurat_atac_cleaned <- cluster_sim_spectrum(seurat_atac_cleaned, label_tag = "sample", use_dr = "lsi", dims = 2:30, corr_method = "pearson", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 30, reduction.name = "css_pca", reduction.key = "CSSPCA_") %>%
  RunUMAP(reduction = "css_pca", dims = 1:20, reduction.name = "umap_css", reduction.key = "UMAPCSS_") %>%

# network propagation to seal annotation for cells with unpaired/multiple labels
seurat_atac <- seurat_atac_cleaned

seurat_atac <- FindNeighbors(seurat_atac, reduction = "lsi", dims = 2:30)
seurat_atac[['peaks_lsi_nn']] <- seurat_atac[['peaks_nn']]
seurat_atac[['peaks_lsi_snn']] <- seurat_atac[['peaks_snn']]
seurat_atac <- FindNeighbors(seurat_atac, reduction = "css_pca", dims = 1:20)
seurat_atac[['peaks_css_nn']] <- seurat_atac[['peaks_nn']]
seurat_atac[['peaks_css_snn']] <- seurat_atac[['peaks_snn']]
seurat_atac <- FindNeighbors(seurat_atac, reduction = "harmony", dims = 1:ncol(Embeddings(seurat_atac,"harmony")))
seurat_atac[['peaks_harmony_nn']] <- seurat_atac[['peaks_nn']]
seurat_atac[['peaks_harmony_snn']] <- seurat_atac[['peaks_snn']]

avg_snn_atac <- (as(seurat_atac[['peaks_lsi_snn']], "Matrix") + as(seurat_atac[['peaks_css_snn']], "Matrix") + as(seurat_atac[['peaks_harmony_snn']], "Matrix"))/3
mat_propagate <- summary(avg_snn_atac)
mat_propagate <- sparseMatrix(i = mat_propagate$i[mat_propagate$i %in% which(seurat_atac$paired_ct %in% c("Unpaired","Multiple"))],
                              j = mat_propagate$j[mat_propagate$i %in% which(seurat_atac$paired_ct %in% c("Unpaired","Multiple"))],
                              x = mat_propagate$x[mat_propagate$i %in% which(seurat_atac$paired_ct %in% c("Unpaired","Multiple"))],
                              dims = dim(avg_snn_atac), dimnames = dimnames(avg_snn_atac))
diag(mat_propagate) <- 1
mat_propagate <- rowNorm(mat_propagate)
mat_ct_ident <- Matrix(as.matrix(df_ct_atac[colnames(seurat_atac),]) + 0, sparse = T)

propagated_ct_mat <- mat_ct_ident
for(i in 1:100)
  propagated_ct_mat <- mat_propagate %*% propagated_ct_mat

seurat_atac$propag_ct <- colnames(propagated_ct_mat)[apply(propagated_ct_mat,1,which.max)]

# exclude non-retinal cell clusters
seurat_atac <- subset(seurat_atac, cells = colnames(seurat_atac)[-which(seurat_atac$paired_ct %in% c("Unpaired","Ast","Dien.","Tel.","Mesenchyme") & seurat_atac$propag_ct %in% c("Ast","Dien.","Tel.","Mesenchyme"))])

seurat_atac <- RunTFIDF(seurat_atac)
det_rates <- rowMeans(seurat_atac@assays$peaks@counts > 0)
VariableFeatures(seurat_atac) <- names(which(det_rates>0.005))
seurat_atac <- RunSVD(seurat_atac) %>%
  RunUMAP(reduction="lsi", dims = 2:30)

### CSS integration
seurat_atac <- cluster_sim_spectrum(seurat_atac, label_tag = "sample", use_dr = "lsi", dims = 2:30, corr_method = "pearson", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 30, reduction.name = "css_pca", reduction.key = "CSSPCA_") %>%
  RunUMAP(reduction = "css_pca", dims = 1:20, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

### Harmony integration
seurat_atac <- RunHarmony(seurat_atac, group.by.vars = "sample", reduction="lsi", dims = 2:30, assay.use = "peaks", max.iter.cluster = 50, max.iter.harmony = 50, project.dim = F)
seurat_atac <- RunUMAP(seurat_atac, reduction="harmony", dims = 1:ncol(Embeddings(seurat_atac,"harmony")), reduction.name = "umap_harmony", reduction.key = "UMAPHARMONY_")

### network propagation for scRNA-seq cell types and niches
seurat_atac <- FindNeighbors(seurat_atac, reduction = "lsi", dims = 2:30)
seurat_atac[['peaks_lsi_nn']] <- seurat_atac[['peaks_nn']]
seurat_atac[['peaks_lsi_snn']] <- seurat_atac[['peaks_snn']]
seurat_atac <- FindNeighbors(seurat_atac, reduction = "css_pca", dims = 1:20)
seurat_atac[['peaks_css_nn']] <- seurat_atac[['peaks_nn']]
seurat_atac[['peaks_css_snn']] <- seurat_atac[['peaks_snn']]
seurat_atac <- FindNeighbors(seurat_atac, reduction = "harmony", dims = 1:ncol(Embeddings(seurat_atac,"harmony")))
seurat_atac[['peaks_harmony_nn']] <- seurat_atac[['peaks_nn']]
seurat_atac[['peaks_harmony_snn']] <- seurat_atac[['peaks_snn']]

avg_snn_atac <- (as(seurat_atac[['peaks_lsi_snn']], "Matrix") + as(seurat_atac[['peaks_css_snn']], "Matrix") + as(seurat_atac[['peaks_harmony_snn']], "Matrix"))/3

###### propagation for retina cell types
df_ct_atac <- setNames(data.frame(sapply(levels(seurat_rna_retina$annot_ct_refine), function(ct){
  rna_cells <- seurat_rna_retina$cell_ident[which(seurat_rna_retina$annot_ct_refine==ct)]
  atac_cells <- rna_atac_pairs$ATAC_cell_ident[rna_atac_pairs$RNA_cell_ident %in% rna_cells]
  seurat_atac$cell_ident %in% atac_cells
}), row.names = colnames(seurat_atac)), levels(seurat_rna_retina$annot_ct_refine))
seurat_atac$paired_ct_ret <- "Unpaired"
seurat_atac$paired_ct_ret[which(rowSums(df_ct_atac)==1)] <- colnames(df_ct_atac)[apply(df_ct_atac[which(rowSums(df_ct_atac)==1),],1,which)]
seurat_atac$paired_ct_ret[which(rowSums(df_ct_atac)>1)] <- "Multiple"
seurat_atac$paired_ct_ret <- factor(seurat_atac$paired_ct_ret, levels=c(levels(seurat_rna_retina$annot_ct_refine),"Multiple","Unpaired"))

mat_propagate <- summary(avg_snn_atac)
mat_propagate <- sparseMatrix(i = mat_propagate$i[mat_propagate$i %in% which(seurat_atac$paired_ct_ret %in% c("Unpaired","Multiple"))],
                              j = mat_propagate$j[mat_propagate$i %in% which(seurat_atac$paired_ct_ret %in% c("Unpaired","Multiple"))],
                              x = mat_propagate$x[mat_propagate$i %in% which(seurat_atac$paired_ct_ret %in% c("Unpaired","Multiple"))],
                              dims = dim(avg_snn_atac), dimnames = dimnames(avg_snn_atac))
diag(mat_propagate) <- 1
mat_propagate <- rowNorm(mat_propagate)
mat_ct_ident <- Matrix::Matrix(as.matrix(df_ct_atac[colnames(seurat_atac),]) + 0, sparse = T)

propagated_ct_mat <- mat_ct_ident
for(i in 1:100)
  propagated_ct_mat <- mat_propagate %*% propagated_ct_mat
seurat_atac$propag_ct_ret <- colnames(propagated_ct_mat)[apply(propagated_ct_mat,1,which.max)]
seurat_atac[['propag_prob_ct_ret']] <- CreateAssayObject(data = t(propagated_ct_mat))

##### convert to major cell type representations
seurat_atac$major_ct <- factor(seurat_atac$propag_ct_ret, levels = intersect(levels(seurat_atac$paired_ct_ret), unique(seurat_atac$propag_ct_ret)))
levels(seurat_atac$major_ct)[grep("^BC-", levels(seurat_atac$major_ct))] <- "BC"
levels(seurat_atac$major_ct)[grep("RPC$", levels(seurat_atac$major_ct))] <- "RPC"
seurat_atac$major_ct <- droplevels(seurat_atac$major_ct, exclude = setdiff(levels(seurat_atac$major_ct), c("RPC","MG","RPE","RGC","AC","HC","BC","Rods","Cones")))


##### propagation for correlation to well-defined cell types (transcriptomeic level)
rna_atac_pairs_retina <- rna_atac_pairs[rna_atac_pairs$RNA_cell_ident %in% seurat_rna_retina$cell_ident & rna_atac_pairs$ATAC_cell_ident %in% seurat_atac$cell_ident,]
mat_rna_atac_retina <- sparseMatrix(i = setNames(1:ncol(seurat_rna_retina),seurat_rna_retina$cell_ident)[rna_atac_pairs_retina$RNA_cell_ident],
                                    j = setNames(1:ncol(seurat_atac),seurat_atac$cell_ident)[rna_atac_pairs_retina$ATAC_cell_ident],
                                    x = 1, dims = c(ncol(seurat_rna_retina),ncol(seurat_atac)), dimnames = list(colnames(seurat_rna_retina),colnames(seurat_atac)))
mat_rna_atac_retina_norm <- colNorm(mat_rna_atac_retina)

corr_rna <- t(seurat_rna_retina@assays$corr_ct@data)
corr_atac_paired <- t(corr_rna) %*% mat_rna_atac_retina_norm

mat_propagate <- summary(avg_snn_atac)
mat_propagate <- sparseMatrix(i = mat_propagate$i[mat_propagate$i %in% which(colSums(mat_rna_atac_retina_norm)==0)],
                              j = mat_propagate$j[mat_propagate$i %in% which(colSums(mat_rna_atac_retina_norm)==0)],
                              x = mat_propagate$x[mat_propagate$i %in% which(colSums(mat_rna_atac_retina_norm)==0)],
                              dims = dim(avg_snn_atac), dimnames = dimnames(avg_snn_atac))
diag(mat_propagate) <- 1
mat_propagate <- rowNorm(mat_propagate)
propagated_corr_mat <- t(corr_atac_paired)
for(i in 1:100)
  propagated_corr_mat <- mat_propagate %*% propagated_corr_mat
seurat_atac[['propag_corr_rna']] <- CreateAssayObject(data = as.matrix(t(propagated_corr_mat)))

saveRDS(seurat_atac, file="data/scATACseq/seurat_objs/onlyRetina_incl_multiome.rds")



# motif enrichment
library(BSgenome.Hsapiens.UCSC.hg38)
load("ext/integrated_motifs.rdata")

seurat_atac <- AddMotifs(seurat_atac, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm_motifs)
seurat_atac <- RunChromVAR(seurat_atac, genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(seurat_atac, file="data/scATACseq/seurat_objs/onlyRetina_incl_multiome.rds")

DefaultAssay(seurat_rna_retina) <- "ATAC"
seurat_rna_retina <- AddMotifs(seurat_rna_retina, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm_motifs)
seurat_rna_retina <- RunChromVAR(seurat_rna_retina, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(seurat_rna_retina) <- "RNA"
saveRDS(seurat_rna_retina, file="data/bimodal/onlyRetina_integrated.rds")



# compare transcriptomic and epigenomic differences between differentiated cell types as well as RPC
idx_atac <- setNames(lapply(list("Rods","Cones",c("BC-ON","BC-OFF"),"AC","HC","RGC","MG","RPE",c("RPC","Prolif. RPC")), function(ct) return(which(seurat_atac$propag_ct %in% ct)) ),
                     c("Rods","Cones","BC","AC","HC","RGC","MG","RPE","RPC"))
idx_rna <- setNames(lapply(list("Rods","Cones",c("BC-ON","BC-OFF"),"AC","HC","RGC","MG","RPE",c("RPC","Prolif. RPC")), function(ct) return(which(seurat_rna_retina$annot_ct_refine %in% ct)) ),
                    c("Rods","Cones","BC","AC","HC","RGC","MG","RPE","RPC"))

avg_expr_ct <- sapply(idx_rna, function(idx) rowMeans(seurat_rna_retina@assays$RNA@data[,idx]))
avg_atac_ct <- sapply(idx_atac, function(idx) rowMeans(seurat_atac@assays$peaks@data[,idx]))

layout(matrix(1:4,nrow=2, byrow=T)); par(mar=c(2,5,4,1))
plot(hclust(as.dist(1-cor(avg_expr_ct[VariableFeatures(seurat_rna_retina),])), method="ward.D2"), xlab=NA, sub=NA, main="RNA Dendrogram (expr/highvar)")
plot(hclust(as.dist(1-cor(avg_atac_ct[VariableFeatures(seurat_atac),])), method="ward.D2"), xlab=NA, sub=NA, main="ATAC Dendrogram (TFIDF/0.5%)")
plot(hclust(as.dist(1-cor(avg_expr_ct)), method="ward.D2"), xlab=NA, sub=NA, main="RNA Dendrogram (expr/all)")
plot(hclust(as.dist(1-cor(avg_atac_ct)), method="ward.D2"), xlab=NA, sub=NA, main="ATAC Dendrogram (TFIDF/all)")




