library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(simspec)
library(dplyr)
library(Signac)
source("ext/util.r")

# 10x scRNA-seq data and QC
counts_samples <- setNames(lapply(rownames(meta_samples)[meta_samples$Type == "RNA"], function(sample)
  Read10X_h5(paste0("data/scRNAseq/processed/", sample, "/outs/filtered_feature_bc_matrix.h5"))),
  rownames(meta_samples)[meta_samples$Type == "RNA"])
seurat_samples <- lapply(names(counts_samples), function(sample){
  meta <- data.frame(sample = rep(sample,ncol(counts_samples[[sample]])), line = meta_samples[sample,"OrganoidLine"], age_weeks = meta_samples[sample,"WeekAge"], orig_cell_id = colnames(counts_samples[[sample]]), row.names = colnames(counts_samples[[sample]]))
  srt <- CreateSeuratObject(counts_samples[[sample]], project = "hRO_timecourse_scRNAseq", meta.data = meta)
  saveRDS(srt, file=paste0("data/scRNAseq/seurat_objs/",sample,".rds"))
})
seurat_rna <- merge(seurat_samples[[1]], seurat_samples[-1])
seurat_rna <- PercentageFeatureSet(seurat_rna, "^MT-", col.name="percent.mt")
seurat_rna <- PercentageFeatureSet(seurat_rna, features=read.table("ext/RPgenes.txt")[,1], col.name="percent.rp")
seurat_rna <- subset(seurat_rna, subset = nCount_RNA > 1500 & nCount_RNA < 20000 & (percent.mt < 20 | (percent.mt < 40 & sample == "GB2_scRNAseq")))

saveRDS(seurat_rna, file = "data/scRNAseq/seurat_objs/timecourse_merged.rds")

# 10x multiome data (scRNA-seq part)
counts_samples <- setNames(lapply(rownames(meta_samples)[meta_samples$Type == "multiome"], function(sample)
  Read10X_h5(paste0("data/multiome/processed/",sample,"/outs/filtered_feature_bc_matrix.h5"))),
  rownames(meta_samples)[meta_samples$Type == "multiome"])
seurat_samples <- setNames(lapply(names(counts_samples), function(sample){
  meta <- data.frame(sample = rep(sample,ncol(counts_samples[[sample]]$`Gene Expression`)), line = meta_samples[sample,"OrganoidLine"], age_weeks = meta_samples[sample,"WeekAge"], orig_cell_id = colnames(counts_samples[[sample]]$`Gene Expression`), row.names = colnames(counts_samples[[sample]]$`Gene Expression`))
  seurat <- CreateSeuratObject(counts_samples[[sample]]$`Gene Expression`, project = "hRO_timecourse_multiome", meta.data = meta)
  seurat$Log10_nCount_RNA <- log10(seurat$nCount_RNA)
  seurat$Log10_nFeature_RNA <- log10(seurat$nFeature_RNA) 
  seurat <- PercentageFeatureSet(seurat, "^MT-", col.name="percent.mt")
  seurat <- PercentageFeatureSet(seurat, features=read.table("ext/RPgenes.txt")[,1], col.name="percent.rp")
  return(seurat)
}), names(counts_samples))

ensdb <- retrieve_ensbd(version = 98) # cellranger-arc 2020-A uses the Ens98/GENCODE32 human annotation
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"
seurat_samples <- setNames(lapply(names(counts_samples), function(sample){
  grange.counts <- StringToGRanges(rownames(counts_samples[[sample]]$Peaks), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- counts_samples[[sample]]$Peaks[as.vector(grange.use), ]
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0("data/multiome/processed/",sample,"/outs/atac_fragments.tsv.gz"),
    min.cells = 1,
    annotation = annotations
  )
  
  seurat <- seurat_samples[[sample]]
  seurat[["ATAC"]] <- chrom_assay
  
  saveRDS(seurat, file=paste0("data/multiome/seurat_objs/",sample,".rds"))
  return(seurat)
}), names(counts_samples))

## QC
seurat_samples$w15_retina_org <- subset(seurat_samples$w15_retina_org, subset = nCount_RNA > 1000 & nCount_RNA < 10000 & nCount_ATAC > 3000 & nCount_ATAC < 30000 & percent.mt < 10)
seurat_samples$w36_retina_org <- subset(seurat_samples$w36_retina_org, subset = nCount_RNA > 500 & nCount_RNA < 20000 & nCount_ATAC > 500 & nCount_ATAC < 20000 & percent.mt < 30)
saveRDS(seurat_samples, file="data/multiome/seurat_objs/sample_QCed_list.rds")

seurat_samples$w15_retina_org[['ATAC']] <- NULL
seurat_samples$w36_retina_org[['ATAC']] <- NULL
seurat_multiome_rna <- merge(seurat_samples$w15_retina_org, seurat_samples$w36_retina_org)
saveRDS(seurat_multiome_rna, file="data/multiome/seurat_objs/rna_QCed.rds")

# Cowan et al. data set
counts_b7_w6_38 <- readRDS("data/Cowan/human_retinal_organoid_f49b7/week_6_38/umi_counts_human.rds")
rownames(counts_b7_w6_38) <- make.unique(readRDS("data/Cowan/human_retinal_organoid_f49b7/week_6_38/gene_anno_human.rds")$symbol)
meta_b7_w6_38 <- data.frame(sample = readRDS("data/Cowan/human_retinal_organoid_f49b7/week_6_38/ega_sample_alias.rds"),
                            line = "B7-YYSH9",
                            age_weeks = NA,
                            orig_cell_id = colnames(counts_b7_w6_38),
                            row.names = colnames(counts_b7_w6_38))
meta_b7_w6_38$age_weeks <- as.numeric(gsub("^W", "", sapply(strsplit(meta_b7_w6_38$sample, "-"), "[", 2)))
seurat_cowan_b7_w6_38 <- CreateSeuratObject(counts_b7_w6_38, project="hRO_Cowan", meta.data = meta_b7_w6_38)

counts_b7_w46 <- readRDS("data/Cowan/human_retinal_organoid_f49b7/week_46/umi_counts_human.rds")
rownames(counts_b7_w46) <- make.unique(readRDS("data/Cowan/human_retinal_organoid_f49b7/week_46/gene_anno_human.rds")$symbol)
meta_b7_w46 <- data.frame(sample = readRDS("data/Cowan/human_retinal_organoid_f49b7/week_46/ega_sample_alias.rds"),
                          line = "B7-YYSH9",
                          age_weeks = 46,
                          orig_cell_id = colnames(counts_b7_w46),
                          row.names = colnames(counts_b7_w46))
seurat_cowan_b7_w46 <- CreateSeuratObject(counts_b7_w46, project="hRO_Cowan", meta.data = meta_b7_w46)

counts_imr90_w30_38 <- readRDS("data/Cowan/human_retinal_organoid_imr90-4/week_30_38/umi_counts_human.rds")
rownames(counts_imr90_w30_38) <- make.unique(readRDS("data/Cowan/human_retinal_organoid_imr90-4/week_30_38/gene_anno_human.rds")$symbol)
meta_imr90_w30_38 <- data.frame(sample = readRDS("data/Cowan/human_retinal_organoid_imr90-4/week_30_38/ega_sample_alias.rds"),
                                line = "IMR90-PG4",
                                age_weeks = NA,
                                orig_cell_id = colnames(counts_imr90_w30_38),
                                row.names = colnames(counts_imr90_w30_38))
meta_imr90_w30_38$age_weeks <- as.numeric(gsub("^W", "", sapply(strsplit(meta_imr90_w30_38$sample, "-"), "[", 2)))
seurat_cowan_imr90_w30_38 <- CreateSeuratObject(counts_imr90_w30_38, project="hRO_Cowan", meta.data = meta_imr90_w30_38)

seurat_cowan <- merge(seurat_cowan_b7_w6_38, list(seurat_cowan_b7_w46, seurat_cowan_imr90_w30_38))
seurat_cowan <- PercentageFeatureSet(seurat_cowan, "^MT-", col.name="percent.mt")
seurat_cowan <- PercentageFeatureSet(seurat_cowan, features=read.table("ext/RPgenes.txt")[,1], col.name="percent.rp")

genes_cowan <- read.table("data/Cowan/annotation/meta_genes.tsv", sep="\t")
idx_sel_genes <- which(genes_cowan[,3] %in% c(""))
seurat_cowan <- subset(seurat_cowan, subset = nCount_RNA > 1000 & nCount_RNA < 20000 & percent.mt < 10)


# merge and integrate the three data sets
seurat <- merge(seurat_rna, list(seurat_multiome_rna, seurat_cowan))
seurat <- NormalizeData(seurat) %>% FindVariableFeatures()
blacklist <- c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat), value=T), read.table("ext/RPgenes.txt")[,1])
highvar_candidates <- setdiff(names(which(table(c(rownames(seurat_rna), rownames(seurat_multiome_rna), rownames(seurat_cowan)))==3)), blacklist)
VariableFeatures(seurat) <- highvar_candidates[order(seurat@assays$RNA@meta.features[highvar_candidates,"vst.variance.standardized"], decreasing=T)[1:3000]]

seurat_rna <- NormalizeData(seurat_rna) %>% FindVariableFeatures()
seurat_multiome_rna <- NormalizeData(seurat_multiome_rna) %>% FindVariableFeatures()
seurat_cowan <- NormalizeData(seurat_cowan) %>% FindVariableFeatures()
seurat@assays$RNA@meta.features$vst.variance.standardized.rna <- NA
seurat@assays$RNA@meta.features$vst.variance.standardized.multiome <- NA
seurat@assays$RNA@meta.features$vst.variance.standardized.cowan <- NA
seurat@assays$RNA@meta.features[rownames(seurat_rna),"vst.variance.standardized.rna"] <- seurat_rna@assays$RNA@meta.features$vst.variance.standardized
seurat@assays$RNA@meta.features[rownames(seurat_multiome_rna),"vst.variance.standardized.multiome"] <- seurat_multiome_rna@assays$RNA@meta.features$vst.variance.standardized
seurat@assays$RNA@meta.features[rownames(seurat_cowan),"vst.variance.standardized.cowan"] <- seurat_cowan@assays$RNA@meta.features$vst.variance.standardized
freq_highvar <- table(highvar_candidates[as.numeric(apply(seurat@assays$RNA@meta.features[highvar_candidates,c("vst.variance.standardized.rna","vst.variance.standardized.multiome","vst.variance.standardized.cowan")], 2, order, decreasing=T)[1:3000,])])
VariableFeatures(seurat) <- intersect(rownames(seurat)[order(seurat@assays$RNA@meta.features$vst.variance.standardized, decreasing=T)], names(which(freq_highvar > 1)))

seurat <- CellCycleScoring(seurat, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)
seurat <- ScaleData(seurat, vars.to.regress = c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
seurat <- cluster_sim_spectrum(seurat, label_tag = "sample", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 20, reduction.name = "css_pca", reduction.key = "CSSPCA_") %>%
  RunUMAP(reduction = "css_pca", dims = 1:20, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

seurat <- FindNeighbors(seurat, reduction = "css_pca", dims = 1:20)
seurat[['RNA_css_nn']] <- seurat[['RNA_nn']]
seurat[['RNA_css_snn']] <- seurat[['RNA_snn']]
seurat <- FindClusters(seurat, graph.name = "RNA_css_snn", resolution = 0.1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 0.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 0.8) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 2)

saveRDS(seurat, file="data/scRNAseq/seurat_objs/all_processed.rds")


# exclude mesenchymal cells
seurat <- subset(seurat, subset = RNA_css_snn_res.0.1 != 6, features = names(which(table(c(rownames(seurat_rna), rownames(seurat_multiome_rna), rownames(seurat_cowan)))==3)))
seurat <- NormalizeData(seurat) %>% FindVariableFeatures()
seurat_datasets <- SplitObject(seurat, "orig.ident")
seurat_datasets <- lapply(seurat_datasets, function(x) NormalizeData(x) %>% FindVariableFeatures())
seurat@assays$RNA@meta.features$vst.variance.standardized.rna <- seurat_datasets$hRO_timecourse_scRNAseq@assays$RNA@meta.features$vst.variance.standardized
seurat@assays$RNA@meta.features$vst.variance.standardized.multiome <- seurat_datasets$hRO_timecourse_multiome@assays$RNA@meta.features$vst.variance.standardized
seurat@assays$RNA@meta.features$vst.variance.standardized.cowan <- seurat_datasets$hRO_Cowan@assays$RNA@meta.features$vst.variance.standardized

highvar_candidates <- setdiff(rownames(seurat), blacklist)
freq_highvar <- table(highvar_candidates[as.numeric(apply(seurat@assays$RNA@meta.features[highvar_candidates,c("vst.variance.standardized.rna","vst.variance.standardized.multiome","vst.variance.standardized.cowan")], 2, order, decreasing=T)[1:3000,])])
VariableFeatures(seurat) <- intersect(rownames(seurat)[order(seurat@assays$RNA@meta.features$vst.variance.standardized, decreasing=T)], names(which(freq_highvar > 1)))

seurat <- ScaleData(seurat, vars.to.regress = c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
seurat <- cluster_sim_spectrum(seurat, label_tag = "sample", cluster_resolution = 1, use_fast_rank=F) %>%
  run_PCA(reduction = "css", npcs = 20, reduction.name = "css_pca", reduction.key = "CSSPCA_") %>%
  RunUMAP(reduction = "css_pca", dims = 1:20, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

seurat <- FindNeighbors(seurat, reduction = "css_pca", dims = 1:20)
seurat[['RNA_css_nn']] <- seurat[['RNA_nn']]
seurat[['RNA_css_snn']] <- seurat[['RNA_snn']]
seurat <- FindClusters(seurat, graph.name = "RNA_css_snn", resolution = 0.1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 0.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 0.8) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 2)

saveRDS(seurat, file="data/scRNAseq/seurat_objs/noMC_processed.rds")

# exclude brain cells
seurat <- subset(seurat_all, subset = RNA_css_snn_res.0.5 != "9" & RNA_css_snn_res.0.5 != "13" & RNA_css_snn_res.0.5 != "15")
seurat <- RunUMAP(seurat, reduction = "css_pca", dims = 1:20, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

saveRDS(seurat, file="data/scRNAseq/seurat_objs/onlyRetina_processed.rds")
