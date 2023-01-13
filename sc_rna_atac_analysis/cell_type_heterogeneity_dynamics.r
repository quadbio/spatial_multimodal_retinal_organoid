library(Matrix)
library(dplyr)
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(simspec)


library(igraph)
library(presto)
library(glmnet)
library(reticulate)
library(gplots)
library(ggplot2)

source("ext/util.r")


# load data
seurat <- readRDS("data/bimodal/onlyRetina_integrated.rds")

# cell type annotation
seurat$annot_ct <- setNames(c("RPC","Rods","MG","Cones","Prolif. RPC","RGC","Early BC","Early IN","AC/HC","Ast","PR-PC/Doublet","PC-MG","BC","Tel.","Early PR","Dien.","RPE","Rods"), levels(seurat$RNA_css_snn_res.0.5))[as.character(seurat$RNA_css_snn_res.0.5)]
seurat$annot_ct_refine <- as.character(seurat$annot_ct)

## Early PR/BC
idx_PR_train <- unlist(lapply(c("Rods","Cones","BC","RPC"), function(x){
  set.seed(30)
  sample(which(seurat$annot_ct == x), 1000)
}))
library(doParallel)
registerDoParallel(10)
model_PR <- cv.glmnet(as.matrix(t(seurat@assays$RNA@data[VariableFeatures(seurat),idx_PR_train])),
                      droplevels(seurat$annot_ct[idx_PR_train]),
                      family = "multinomial",
                      parallel = T)
stopImplicitCluster()

pred_PR <- predict(model_PR, as.matrix(t(seurat@assays$RNA@data[VariableFeatures(seurat),which(seurat$annot_ct %in% c("Early PR/BC","PR-PC/Doublet"))])), type="response")[,,1]
seurat$annot_ct_refine[intersect(colnames(seurat)[seurat$annot_ct == "Early PR/BC"], rownames(pred_PR)[apply(pred_PR,1,which.max)==which(colnames(pred_PR)=="BC") & apply(pred_PR,1,max)>0.5])] <- "Immature BC"
seurat$annot_ct_refine[intersect(colnames(seurat)[seurat$annot_ct == "Early PR/BC"], rownames(pred_PR)[apply(pred_PR,1,which.max)==which(colnames(pred_PR)=="Rods") & apply(pred_PR,1,max)>0.5])] <- "Immature rods"
seurat$annot_ct_refine[intersect(colnames(seurat)[seurat$annot_ct == "Early PR/BC"], rownames(pred_PR)[apply(pred_PR,1,which.max)==which(colnames(pred_PR)=="Cones") & apply(pred_PR,1,max)>0.5])] <- "Immature cones"
seurat$annot_ct_refine[seurat$annot_ct_refine == "Early PR/BC"] <- "PR/BC-PC"

## PR-PC/Doublet
idx_stage_train <- unlist(lapply(list(which(seurat$annot_ct=="PR-PC/Doublet" & seurat$age_weeks<15),
                                      which(seurat$annot_ct=="PR-PC/Doublet" & seurat$age_weeks>30)), function(x){
                                        set.seed(30)
                                        sample(x, 1000)
                                      }))
library(doParallel)
registerDoParallel(10)
model_stage <- cv.glmnet(as.matrix(t(seurat@assays$RNA@data[VariableFeatures(seurat),idx_stage_train])),
                         ifelse(seurat$age_weeks[idx_stage_train]<15, "early", "late"),
                         family = "binomial",
                         parallel = T)
stopImplicitCluster()
pred_stage <- predict(model_stage, as.matrix(t(seurat@assays$RNA@data[VariableFeatures(seurat),which(seurat$annot_ct == "PR-PC/Doublet")])), type="response")
seurat$annot_ct_refine[rownames(pred_stage)[pred_stage[,1]<0.5]] <- "Early PR/BC-IPC"
seurat$annot_ct_refine[rownames(pred_stage)[pred_stage[,1]>0.5]] <- "Late PR/BC-IPC"

## AC/HC
seurat_IN <- subset(seurat, subset = annot_ct == "AC/HC") %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 10, verbose=F) %>%
  RunUMAP(dims = 1:10)
seurat_IN_mnn <- RunFastMNN(SplitObject(seurat_IN, "sample"))
seurat_IN[['mnn']] <- CreateDimReducObject(Embeddings(seurat_IN_mnn,"mnn")[colnames(seurat_IN),], key="MNN_", assay = "RNA")
seurat_IN <- RunUMAP(seurat_IN, reduction="mnn", dims = 1:10, reduction.name="umap_mnn", reduction.key="UMAPMNN_")
seurat_IN <- FindNeighbors(seurat_IN, reduction="mnn", dims=1:10) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 0.8) %>%
  FindClusters(resolution = 1)
seurat_IN$annot_ct_refine[seurat_IN$RNA_snn_res.0.5 %in% c("6","3")] <- "AC"
seurat_IN$annot_ct_refine[seurat_IN$RNA_snn_res.0.5 %in% c("4","7","5","2")] <- "HC"

seurat$annot_ct_refine[colnames(seurat_IN)] <- seurat_IN$annot_ct_refine
seurat$annot_ct_refine[seurat$annot_ct_refine == "Early IN"] <- "Immature IN"

## RGC
seurat_RGC <- subset(seurat, subset = annot_ct == "RGC") %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 10, verbose=F) %>%
  RunUMAP(dims = 1:10)
seurat_RGC_mnn <- RunFastMNN(SplitObject(seurat_RGC, "sample"))
seurat_RGC[['mnn']] <- CreateDimReducObject(Embeddings(seurat_RGC_mnn,"mnn")[colnames(seurat_RGC),], key="MNN_", assay = "RNA")
seurat_RGC <- RunUMAP(seurat_RGC, reduction="mnn", dims = 1:10, reduction.name="umap_mnn", reduction.key="UMAPMNN_")
seurat_RGC <- FindNeighbors(seurat_RGC, reduction="mnn", dims=1:10) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.3) %>%
  FindClusters(resolution = 0.5)
seurat_RGC$annot_ct_refine <- NA
seurat_RGC$annot_ct_refine[seurat_RGC$RNA_snn_res.0.3 %in% c(2,4)] <- "RGC"
seurat_RGC$annot_ct_refine[! seurat_RGC$RNA_snn_res.0.3 %in% c(2,4)] <- "Immature IN"

seurat$annot_ct_refine[colnames(seurat_RGC)] <- seurat_RGC$annot_ct_refine

## BC
seurat_BC <- subset(seurat, subset = annot_ct == "BC") %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 10, verbose=F) %>%
  RunUMAP(dims = 1:10)
seurat_BC_mnn <- RunFastMNN(SplitObject(seurat_BC, "sample"))
seurat_BC[['mnn']] <- CreateDimReducObject(Embeddings(seurat_BC_mnn,"mnn")[colnames(seurat_BC),], key="MNN_", assay = "RNA")
seurat_BC <- RunUMAP(seurat_BC, reduction="mnn", dims = 1:10, reduction.name="umap_mnn", reduction.key="UMAPMNN_")

seurat_BC <- FindNeighbors(seurat_BC, reduction="mnn", dims=1:10) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 0.8) %>%
  FindClusters(resolution = 1)
seurat_BC$annot_ct_refine[seurat_BC$RNA_snn_res.0.1 %in% c("1","2")] <- "BC-ON"
seurat_BC$annot_ct_refine[seurat_BC$RNA_snn_res.0.1 %in% c("0")] <- "BC-OFF"

seurat$annot_ct_refine[colnames(seurat_BC)] <- seurat_BC$annot_ct_refine



# focus on the well-defined cell types: RPC, MG, RPE, RGC, AC, HC, BC, Rods, Cones
seurat$major_ct <- seurat$annot_ct_refine
levels(seurat$major_ct)[grep("^BC-", levels(seurat$major_ct))] <- "BC"
levels(seurat$major_ct)[grep("RPC$", levels(seurat$major_ct))] <- "RPC"
seurat$major_ct <- droplevels(seurat$major_ct, exclude = setdiff(levels(seurat$major_ct), c("RPC","MG","RPE","RGC","AC","HC","BC","Rods","Cones")))

seurat_ct <- DietSeurat(seurat, assays = "RNA") %>%
  subset(cells = colnames(seurat)[!is.na(seurat$major_ct)])

### expression features of cell types
DE_ct <- presto::wilcoxauc(seurat_ct, "major_ct")
DEG_ct <- lapply(split(DE_ct, DE_ct$group), function(x){
  x <- x[which(x$padj < 0.01 & x$auc > 0.7 & x$logFC > log(1.2) & x$pct_in - x$pct_out > 20 & x$pct_out < 20 & !x$feature %in% blacklist),]
  x <- x[order(x$pct_out),]
  return(x)
})
tfs <- read.table("~/Work/databases/animalTFDB/animalTFDB_v3/Homo_sapiens_TF.tsv", sep="\t", header=T)
DEG_tf_ct <- lapply(split(DE_ct, DE_ct$group), function(x){
  x <- x[which(x$padj < 0.01 & x$auc > 0.65 & x$logFC > log(1.2) & x$pct_in - x$pct_out > 20 & x$pct_out < 20 & !x$feature %in% blacklist & x$feature %in% tfs$Symbol),]
  x <- x[order(x$pct_out),]
  return(x)
})

### accessibility features of cell types
seurat_ct_withatac <- subset(seurat_ct, subset = paired_ATAC & orig.ident == "hRO_timecourse_scRNAseq")
DA_ct <- wilcoxauc(seurat_ct_withatac, "major_ct", assay = "data", seurat_assay = "ATAC_stage")

DAP_ct <- lapply(split(DA_ct,DA_ct$group), function(x){
  x <- x[which(x$padj < 0.01 & x$auc > 0.51 & x$pct_in / x$pct_out > 5),]
  x <- x[order(x$pct_in - x$pct_out, decreasing=T),]
  return(x)
})


# visualization of cell type annotation: major cell type centric
avg_major_ct <- sapply(levels(seurat$major_ct), function(ct) rowMeans(seurat@assays$RNA@data[,which(seurat$major_ct == ct)]))
corr2ct <- qlcMatrix::corSparse(seurat@assays$RNA@data[VariableFeatures(seurat),], avg_major_ct[VariableFeatures(seurat),])
dimnames(corr2ct) <- list(colnames(seurat), levels(seurat$major_ct))

seurat[['corr_ct']] <- CreateAssayObject(data = t(corr2ct))

### cor kernel probability to mix colors
lambda <- 50
prob <- t(apply(exp(corr2ct * lambda) * exp(-lambda), 1, function(x) x/sum(x)))
mat_cols_ct <- t(col2rgb(cols_major_ct))[colnames(prob),]
mat_cols <- prob %*% mat_cols_ct
cols <- (mat_cols * 1 + 222 * 1) / 2
cols <- apply(cols,1,function(x) rgb(x[1],x[2],x[3], maxColorValue = 255))
cols <- paste0(cols, as.hexmode(as.integer(apply(prob,1,max)*150)))
cols[!is.na(seurat$major_ct)] <- cols_major_ct[as.character(seurat$major_ct[!is.na(seurat$major_ct)])]

par(mar=c(1,1,1,1), cex=0.5)
plotFeature(Embeddings(seurat,"umap_css")[order(!is.na(seurat$major_ct)),], col = cols[order(!is.na(seurat$major_ct))], cex=0.3, pt_border = F, lwd_border=0.01, random_order = F)
legend("topleft", legend = names(cols_major_ct), pch = 16, col = cols_major_ct, pt.cex = 2, cex = 1, bty="n")











# scvelo and pseudotime analysis
matrices_dropest_RNA <- lapply(sort(unique(seurat$sample[which(seurat$orig.ident=="hRO_timecourse_scRNAseq")])), function(sample){
  message(paste0("start ", sample))
  cell_names <- setNames(colnames(seurat)[seurat$sample==sample], gsub("_.+$", "", colnames(seurat)[seurat$sample==sample]))
  matrices <- readRDS(paste0("data/scRNAseq/dropest/",sample,"/",sample,".matrices.rds"))
  matrices <- lapply(matrices[c("exon","intron")], function(mat){
    mat <- mat[,names(cell_names)]
    colnames(mat) <- cell_names[colnames(mat)]
    return(mat)
  })
  matrices <- lapply(matrices, function(x) x[intersect(rownames(matrices$exon), rownames(matrices$intron)),])
  return(matrices)
})
genes_shared <- names(which(table(unlist(lapply(matrices_dropest_RNA, function(x) rownames(x$exon))))==length(matrices_dropest_RNA)))
genes_scvelo <- intersect(genes_shared, VariableFeatures(seurat))

mat_exonic <- do.call(cbind, lapply(matrices_dropest_RNA, function(x) x$exon[genes_shared,]))
mat_intronic <- do.call(cbind, lapply(matrices_dropest_RNA, function(x) x$intron[genes_shared,]))
seurat_timecourse <- subset(seurat, cells = colnames(mat_exonic))
seurat_timecourse$paired_ATAC <- factor(seurat_timecourse$paired_ATAC)
h5ad_timecourse <- data_to_h5ad(object = t(seurat_timecourse@assays$RNA@data[genes_scvelo,]),
                                vars = seurat_timecourse@assays$RNA@meta.features[genes_scvelo,],
                                obs = seurat_timecourse@meta.data,
                                obsm = list(pca = Embeddings(seurat_timecourse,"pca"),
                                            css = Embeddings(seurat_timecourse,"css"),
                                            csspca = Embeddings(seurat_timecourse,"css_pca"),
                                            umap = Embeddings(seurat_timecourse,"umap_css")),
                                layers = list(spliced = t(mat_exonic[genes_scvelo,]),
                                              unspliced = t(mat_intronic[genes_scvelo,])),
                                savefile = "data/scRNAseq/h5ad/onlyRetina_timecourse.h5ad")

## run scvelo
scvelo <- import("scvelo")
anndata <- import("anndata")
pd <- import("pandas")
plt <- import("matplotlib")
plt$use("Agg", force = TRUE)

h5ad_timecourse <- scvelo$read("data/scRNAseq/h5ad/onlyRetina_timecourse.h5ad")
scvelo$pp$moments(h5ad_timecourse,
                  n_neighbors = as.integer(30),
                  use_rep = "csspca")
scvelo$tl$velocity(h5ad_timecourse,
                   groupby = "sample")
scvelo$tl$velocity_graph(h5ad_timecourse)
scvelo$tl$velocity_pseudotime(h5ad_timecourse)

scvelo$pl$velocity_embedding_stream(h5ad_timecourse,
                                    basis="umap",
                                    color="annot_ct_refine",
                                    dpi=120,
                                    figsize = c(8,8),
                                    save="scvelo_stream.png")
scvelo$pl$velocity_embedding_grid(h5ad_timecourse,
                                  density = 1,
                                  arrow_length = 3,
                                  basis="umap",
                                  color="annot_ct_refine",
                                  dpi=120,
                                  figsize = c(8,8),
                                  save="scvelo_grid.png")
h5ad_timecourse$write_h5ad("data/scRNAseq/h5ad/onlyRetina_timecourse.h5ad")

## embed scvelo results into the seurat object
h5ad_timecourse <- anndata::read_h5ad("data/scRNAseq/h5ad/onlyRetina_timecourse.h5ad")
seurat_timecourse$velocity_pseudotime <- h5ad_timecourse$obs$velocity_pseudotime
seurat_timecourse$velocity_self_transition <- h5ad_timecourse$obs$velocity_self_transition
seurat_timecourse$root_cells <- h5ad_timecourse$obs$root_cells
seurat_timecourse$end_points <- h5ad_timecourse$obs$end_points

## cellrank in python, with the h5ad with velocity included as the input
system("python run_cellrank.py")

## next, read in the results and embed into the seurat object
cellrank_fates <- read.table("data/scRNAseq/cellrank_lineages_absorption_prob.tsv")[colnames(seurat_timecourse),]
colnames(cellrank_fates) <- paste0("cellrank_prob_", colnames(cellrank_fates))
seurat_timecourse <- AddMetaData(seurat_timecourse, cellrank_fates)
saveRDS(seurat_timecourse, file="data/scRNAseq/seurat_objs/onlyRetina_timecourse.rds")

## transfer info to other cells (network propagation)
knn <- RANN::nn2(Embeddings(seurat, "css_pca"), k=50)

idx2est <- which(! colnames(seurat) %in% colnames(seurat_timecourse))
mat_propag <- sparseMatrix(i = c(rep(idx2est, ncol(knn$nn.idx)-1), setdiff(1:nrow(knn$nn.idx), idx2est)),
                           j = c(as.numeric(knn$nn.idx[idx2est,-1]), setdiff(1:nrow(knn$nn.idx), idx2est)),
                           x = 1,
                           dims = c(nrow(knn$nn.idx), nrow(knn$nn.idx)),
                           dimnames = list(colnames(seurat), colnames(seurat)))
mat_propag_norm <- t(colNorm(t(mat_propag)))

init_prob <- matrix(1/ncol(cellrank_fates), ncol = ncol(cellrank_fates), nrow = ncol(seurat))
dimnames(init_prob) <- list(colnames(seurat), colnames(cellrank_fates))
init_prob[rownames(cellrank_fates),] <- as.matrix(cellrank_fates)

init_pt_info <- matrix(0.5, ncol=4, nrow = ncol(seurat))
dimnames(init_pt_info) <- list(colnames(seurat), colnames(seurat_timecourse@meta.data)[31:34])
init_pt_info[colnames(seurat_timecourse),] <- as.matrix(seurat_timecourse@meta.data[,31:34])

propagated_prob <- Matrix(init_prob)
for(i in 1:100)
  propagated_prob <- mat_propag_norm %*% propagated_prob
propagated_prob_ranked <- apply(propagated_prob, 2, rank)
seurat <- AddMetaData(seurat, data.frame(as.matrix(propagated_prob)))

propagated_pt <- Matrix(init_pt_info)
for(i in 1:100)
  propagated_pt <- mat_propag_norm %*% propagated_pt
seurat <- AddMetaData(seurat, data.frame(as.matrix(propagated_pt)))

saveRDS(seurat, file="data/bimodal/onlyRetina_integrated.rds")








# cell state niche summary (cluster with res=20) of the atlas
## high-resolution clustering
cols_terminal <- c("Rods" = "#CB4335",
                   "Cones" = "#E67E22",
                   "BC" = "#D4AC0D",
                   "AC" = "#1D8348",
                   "HC" = "#17A589",
                   "RGC" = "#117A65",
                   "MG" = "#2980B9",
                   "RPE" = "#2C3E50")

seurat_rna <- DietSeurat(seurat, assays="RNA", dimreducs = "css_pca")
seurat_AC <- subset(seurat_rna, subset = annot_ct_refine == "AC")
seurat_HC <- subset(seurat_rna, subset = annot_ct_refine == "HC")
seurat_others <- subset(seurat_rna, subset = annot_ct_refine %in% setdiff(levels(seurat_rna$annot_ct_refine), c("AC","HC")))

seurat_AC <- FindNeighbors(seurat_AC, reduction = "css_pca", dims = 1:20) %>%
  FindClusters(resolution = 0.5)
seurat_HC <- FindNeighbors(seurat_HC, reduction = "css_pca", dims = 1:20) %>%
  FindClusters(resolution = 0.5)
seurat_others <- FindNeighbors(seurat_others, reduction = "css_pca", dims = 1:20) %>%
  FindClusters(resolution = 20)

levels(seurat_AC@active.ident) <- paste0(levels(seurat_AC@active.ident), "_A")
levels(seurat_HC@active.ident) <- paste0(levels(seurat_HC@active.ident), "_H")
levels(seurat_others@active.ident) <- paste0(levels(seurat_others@active.ident), "_O")
seurat$highres_cl_lab <- factor(setNames(c(as.character(seurat_AC@active.ident), as.character(seurat_HC@active.ident), as.character(seurat_others@active.ident)), c(colnames(seurat_AC),colnames(seurat_HC),colnames(seurat_others))),
                                levels=c(levels(seurat_others@active.ident),levels(seurat_AC@active.ident),levels(seurat_HC@active.ident)))[colnames(seurat)]

saveRDS(seurat, file="data/bimodal/onlyRetina_integrated.rds")

## cluster abstraction
avg_niche <- AverageExpression(seurat, group.by="highres_cl_lab", assays = c("RNA","ATAC"))
expr_niche <- avg_niche$RNA
acc_niche <- avg_niche$ATAC
meta_niche <- summarize_data_to_groups(seurat@meta.data, groups = seurat$highres_cl_lab)

## graph edges filtering
### raw edges: cellrank probabilities
annot_niche_ct <- setNames(meta_niche$annot_ct_refine, rownames(meta_niche))
levels(annot_niche_ct)[grep("^BC-",levels(annot_niche_ct))] <- "BC"

fate_prob <- meta_niche[,grep("cellrank_prob", colnames(meta_niche))]
fate_prob_rank <- apply(fate_prob, 2, rank)

knn_niche <- RANN::nn2(fate_prob, k = 20)
adj_knn_niche <- sparseMatrix(i = rep(1:nrow(knn_niche$nn.idx), ncol(knn_niche$nn.idx)),
                              j = as.numeric(knn_niche$nn.idx),
                              dims = c(nrow(meta_niche), nrow(meta_niche)), dimnames = list(rownames(meta_niche),rownames(meta_niche)))
adj_knn_niche <- adj_knn_niche | t(adj_knn_niche)
edges_niche <- summary(adj_knn_niche)

### filter 0: terminal states should be disconnected
edges_niche <- edges_niche[-which(annot_niche_ct[edges_niche$i] != annot_niche_ct[edges_niche$j] &
                                    annot_niche_ct[edges_niche$i] %in% names(cols_terminal) & 
                                    annot_niche_ct[edges_niche$j] %in% names(cols_terminal) &
                                    apply(fate_prob,1,max)[edges_niche$i] >= 0.5 &
                                    apply(fate_prob,1,max)[edges_niche$j] >= 0.5),]

### filter 1: pseudotime & terminal states being terminal
idx_max_pt_terminals <- sapply(names(cols_terminal), function(ct) which(annot_niche_ct==ct)[which.max(meta_niche$velocity_pseudotime[which(annot_niche_ct==ct)])])
df_edges <- data.frame(from = edges_niche$i, to = edges_niche$j)
df_edges <- df_edges[meta_niche$velocity_pseudotime[df_edges$from] < meta_niche$velocity_pseudotime[df_edges$to],]
df_edges <- df_edges[! df_edges$from %in% idx_max_pt_terminals,]
df_edges <- df_edges[which((! annot_niche_ct[df_edges$from] %in% names(cols_terminal)) |
                             apply(fate_prob,1,max)[df_edges$from] < 0.5 |
                             annot_niche_ct[df_edges$from] == annot_niche_ct[df_edges$to]),]

adj_niche <- sparseMatrix(i = df_edges$from,
                          j = df_edges$to,
                          dims = c(nrow(meta_niche), nrow(meta_niche)), dimnames = list(rownames(meta_niche), rownames(meta_niche)))

### filter 2: PAGA connectivity
meta <- seurat@meta.data
meta$paired_ATAC <- factor(meta$paired_ATAC)
h5ad <- data_to_h5ad(object = t(seurat@assays$RNA@data),
                     vars = seurat@assays$RNA@meta.features,
                     obs = meta,
                     obsm = list(pca = Embeddings(seurat,"pca"),
                                 css = Embeddings(seurat,"css"),
                                 csspca = Embeddings(seurat,"css_pca"),
                                 umap = Embeddings(seurat,"umap_css")),
                     savefile = "data/scRNAseq/h5ad/onlyRetina.h5ad")

sc <- import("scanpy")
sc$pp$neighbors(h5ad, n_neighbors=20L, use_rep='csspca')
sc$tl$paga(h5ad, groups='highres_cl_lab')
h5ad$write_h5ad("data/scRNAseq/h5ad/onlyRetina_RNA.h5ad")

niche_paga_conn <- h5ad$uns$paga$connectivities$tocsc()
dimnames(niche_paga_conn) <- list(levels(h5ad$obs$highres_cl_lab), levels(h5ad$obs$highres_cl_lab))

adj_niche <- adj_niche & (niche_paga_conn > 0.2)
df_edges <- summary(drop0(adj_niche))

### get the graph and generate force-directed layout
df_graph <- data.frame(from = rownames(meta_niche)[df_edges$i],
                       to = rownames(meta_niche)[df_edges$j])
graph_niche <- graph_from_data_frame(df_graph)
set.seed(30)
coord_niche <- layout_nicely(graph_niche)
rownames(coord_niche) <- names(V(graph_niche))
coord_niche <- coord_niche[rownames(meta_niche),]

### check the graph layout
layout(matrix(1:6,nrow=2, byrow=T)); par(mar=c(1,1,1,1))
plotFeature(coord_niche, meta_niche$annot_ct_refine,
            edges = df_graph[,1:2], lwd.edges = 0.3, col.edges = "#bdbdbd", colorPal = cols_ct_refine, cex=1.3, do_legend = T, legend_pos = "topleft", legend_cex = 0.7, pt_border = T)
plotFeature(coord_niche, meta_niche$transferred_ct_dev,
            edges = edges_niche[,1:2], col.edges = "#bdbdbd", colorPal = rev(setNames(prettyrainbow_colscheme(11), levels(meta_niche$transferred_ct_dev)[c(11,7,8,9,1,2,6,3,4,5,10)])), cex=1.3, do_legend = T, legend_pos = "topleft", legend_cex = 0.7, pt_border = T)
plotFeature(coord_niche, annot_niche_ct,
            edges = edges_niche[,1:2], col.edges = "#bdbdbd", colorPal = cols_terminal, cex=1.3, do_legend = T, legend_pos = "topleft", legend_cex = 0.7, pt_border = T, emphasize = which(annot_niche_ct %in% names(cols_terminal)))
plotFeature(coord_niche, meta_niche$velocity_pseudotime, edges = edges_niche[,1:2], col.edges = "#bdbdbd", colorPal = bluewhitered_colscheme, pt_border = T, cex = 1.3)
plotFeature(coord_niche, meta_niche$age_weeks, edges = edges_niche[,1:2], col.edges = "#bdbdbd", colorPal = prettyrainbow_colscheme, pt_border = T, cex = 1.3)

## random walk from root
### define roots: nodes with smallest pt that can reach out to most of the other nodes
graph_dist_niches <- distances(graph_niche, mode = "out")[rownames(meta_niche)[order(meta_niche$velocity_pseudotime)],rownames(meta_niche)[order(meta_niche$velocity_pseudotime)]]
reachable_niches <- sapply(1:nrow(graph_dist_niches), function(i) colSums(!is.infinite(graph_dist_niches[rep(1:i,2),]))>0)
num_roots <- min(which(colSums(reachable_niches) >= length(V(graph_niche))*0.975))
nodes_roots <- rownames(meta_niche)[order(meta_niche$velocity_pseudotime)[1:num_roots]]

roots_random_walk <- lapply(nodes_roots, function(v){
  t(sapply(1:100000, function(i){
    set.seed(i)
    path <- random_walk(graph_niche, start = v, steps = 300, mode = "out")
    return(rownames(meta_niche) %in% names(path))
  }))
})

roots_random_walk_all <- Matrix(do.call(rbind, roots_random_walk)[,order(meta_niche$velocity_pseudotime)]+0, sparse = T)
colnames(roots_random_walk_all) <- rownames(meta_niche)[order(meta_niche$velocity_pseudotime)]
df_root_random_walk <- summary(roots_random_walk_all)
idx_ends_random_walks <- order(meta_niche$velocity_pseudotime)[sapply(split(df_root_random_walk, df_root_random_walk$i), function(x) max(x$j))]
ct_ends_random_walk <- annot_niche_ct[idx_ends_random_walks]
levels(ct_ends_random_walk)[levels(ct_ends_random_walk) %in% c("BC-ON","BC-OFF")] <- "BC"
idx_ends_terminal <- which(ct_ends_random_walk %in% names(cols_terminal))

#### estimate multipotency
prop_ct_path <- sapply(levels(droplevels(ct_ends_random_walk[idx_ends_terminal])), function(ct){
  idx <- which(ct_ends_random_walk == ct)
  prop <- setNames(rep(0, nrow(meta_niche)), rownames(meta_niche))
  prop[order(meta_niche$velocity_pseudotime)] <- colMeans(roots_random_walk_all[idx,])
  return(prop)
})
prop_ct_path_norm <- rowNorm(prop_ct_path)

cutoff_prop <- 0.01
plotFeature(coord_niche, factor(rowSums(prop_ct_path_norm>cutoff_prop)),
            edges = edges_niche[,1:2], col.edges = "#bdbdbd", cex=1.3, pt_border = T, colorPal = prettyrainbow_colscheme, do_legend=T, legend_pos = "topleft", legend_cex = 0.7)
apply(prop_ct_path_norm>cutoff_prop, 2, function(x) apply(prop_ct_path_norm>cutoff_prop, 2, function(y) sum(x&y,na.rm=T)))
apply(prop_ct_path_norm[apply(prop_ct_path_norm>cutoff_prop,1,sum,na.rm=T)==1,]>cutoff_prop, 2, sum)


## save the results
meta_niche$init_ct <- annot_niche_ct
dat_niche <- list(expr = expr_niche,
                  peak = acc_niche,
                  peak_stages = acc2_niche,
                  meta = meta_niche,
                  meta_peaks = seurat@assays$ATAC@meta.features,
                  reductions = list(coord = coord_niche,
                                    prop_terminal_path = prop_ct_path,
                                    prop_terminal_path_norm = prop_ct_path_norm),
                  edge = list(raw = edges_niche[,1:2],
                              prunned = df_edges[,1:2]),
                  adjmat = list(raw = adj_knn_niche,
                                prunned = adj_niche),
                  rw = roots_random_walk_all)

annot_niche <- setNames(dat_niche$meta$annot_ct_refine, rownames(dat_niche$meta))
levels(annot_niche)[levels(annot_niche) %in% c("BC-ON","BC-OFF")] <- "BC"
annot_niche <- droplevels(annot_niche, exclude = setdiff(levels(annot_niche), names(cols_terminal)))

ct_prob <- summary(prop_ct_path_norm)
ct_prob <- ct_prob[ct_prob$x > 0.01,]
ct_prob <- ct_prob[ct_prob$i %in% which(is.na(annot_niche)),]
ct_prob <- rbind(ct_prob, data.frame(i = which(!is.na(annot_niche)),
                                     j = setNames(1:ncol(prop_ct_path_norm), colnames(prop_ct_path_norm))[as.character(annot_niche[!is.na(annot_niche)])],
                                     x = 1))
ct_prob <- sparseMatrix(i = ct_prob$i, j = ct_prob$j, x = ct_prob$x, dims = dim(prop_ct_path_norm), dimnames = dimnames(prop_ct_path_norm))
ct_prob <- rowNorm(ct_prob)

mat_cols_terminal <- t(col2rgb(cols_terminal))
mat_cols_niche <- ct_prob %*% mat_cols_terminal
mat_cols_niche[which(is.na(annot_niche)),] <- (mat_cols_niche[which(is.na(annot_niche)),] + 255) / 2
cols_niche <- apply(mat_cols_niche, 1, function(x) rgb(x[1],x[2],x[3],maxColorValue = 255))
dat_niche$meta$col_niche <- cols_niche[rownames(dat_niche$meta)]

saveRDS(dat_niche, file="data/scRNAseq/graph_abstraction.rds")


