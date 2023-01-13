library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)

source("ext/util.r")

# load scRNA-seq data (output by script.integrate_all_scRNAseq.r)
seurat_rna <- readRDS(file="data/scRNAseq/seurat_objs/all_processed.rds")
seurat_rna$cell_ident <- paste0(sapply(strsplit(colnames(seurat_rna),"-"),"[",1), "_", seurat_rna$sample)

# load multiome with both RNA and ATAC
seurat_multiome_samples <- readRDS("data/multiome/seurat_objs/sample_QCed_list.rds")

# load aggregated scATAC-seq data
seurat_atac_aggr <- create_atac_seurat(path = "data/scATACseq/processed/aggr",
                                       project = "hRO_timecourse_ATAC", annot_ens_version = 93)
meta_samples_aggr <- read.csv("data/scATACseq/aggr/processed/outs/aggregation_csv.csv")
sample_info <- read.table("data/genomic_sample_info_full.tsv", sep="\t", header=T)
seurat_atac_aggr$sample <- meta_samples_aggr$library_id[as.numeric(sapply(strsplit(colnames(seurat_atac_aggr),"-"),"[",2))]
seurat_atac_aggr$cell_ident <- paste0(sapply(strsplit(colnames(seurat_atac_aggr),"-"),"[",1),"_",seurat_atac_aggr$sample)
seurat_atac_aggr$age_weeks <- sample_info[seurat_atac_aggr$sample,"WeekAge"]

saveRDS(seurat_atac_aggr, file="data/scATACseq/seurat_objs/aggr.rds")

# load RNA-ATAC cell pairing results
rna_atac_pairs <- readRDS("data/bimodal/paired_cells_rna_atac.rds")
rna_atac_pairs$RNA_cell_ident <- paste0(sapply(strsplit(rna_atac_pairs$RNA_cell_barcode,"-"),"[",1), "_", rna_atac_pairs$RNA_lib_idx)
rna_atac_pairs$RNA_cell_ident[is.na(rna_atac_pairs$RNA_cell_barcode)] <- NA
rna_atac_pairs$ATAC_cell_ident <- paste0(sapply(strsplit(rna_atac_pairs$ATAC_cell_barcode,"-"),"[",1), "_", rna_atac_pairs$ATAC_lib_idx)
rna_atac_pairs$ATAC_cell_ident[is.na(rna_atac_pairs$ATAC_cell_barcode)] <- NA
rna_atac_pairs <- data.frame(dataset = "RNA_ATAC_paired",
                             rna_atac_pairs[which(!is.na(rna_atac_pairs$RNA_cell_ident) & !is.na(rna_atac_pairs$ATAC_cell_ident) & rna_atac_pairs$RNA_cell_ident %in% seurat_rna$cell_ident),c(3,5,4,6)])

multiome_pairs <- data.frame(dataset = "multiome",
                             RNA_lib_idx = seurat_rna$sample[seurat_rna$orig.ident == "hRO_timecourse_multiome"],
                             RNA_cell_ident = seurat_rna$cell_ident[seurat_rna$orig.ident == "hRO_timecourse_multiome"],
                             ATAC_lib_idx = seurat_rna$sample[seurat_rna$orig.ident == "hRO_timecourse_multiome"],
                             ATAC_cell_ident = seurat_rna$cell_ident[seurat_rna$orig.ident == "hRO_timecourse_multiome"], row.names = NULL)

rna_atac_pairs <- rbind(rna_atac_pairs, multiome_pairs)
rownames(rna_atac_pairs) <- NULL

# remove low-quality integration result
rna_atac_pairs <- rna_atac_pairs[which(! rna_atac_pairs$ATAC_lib_idx %in% c("Lib203","Lib206","Lib87_C7_w7_B7","Lib102","Lib205","Lib91_G7_w10_409b2")),]

# requantify ATAC data: split samples based on time points and call peaks for each group separately, and then merge
## scATAC-seq
seurat_atac_aggr$age_group <- ceiling(seurat_atac_aggr$age_weeks/10)
seurat_atac_aggr$age_group[seurat_atac_aggr$age_group > 4] <- 4
peaks <- CallPeaks(seurat_atac_aggr, group.by="age_group")

annotations <- retrieve_ensbd(version = 93)
counts_atac_aggr <- FeatureMatrix(seurat_atac_aggr@assays$peaks@fragments,
                                  features = peaks,
                                  cells = colnames(seurat_atac_aggr))
assay_peaks_stage_atac <- CreateChromatinAssay(counts_atac_aggr, fragments = seurat_atac_aggr@assays$peaks@fragments)
Annotation(assay_peaks_stage_atac) <- annotations
seurat_atac_aggr[['peaks']] <- assay_peaks_stage_atac

saveRDS(seurat_atac_aggr, file="data/scATACseq/seurat_objs/aggr_requant.rds")

## scMultiome
counts_peaks_multiome <- lapply(seurat_multiome_samples, function(x)
  FeatureMatrix(x@assays$ATAC@fragments,
                features = seurat_atac_aggr@assays$peaks@ranges,
                cells = colnames(x)))
seurat_multiome_atac <- setNames(lapply(names(counts_peaks_multiome), function(sample){
  chr_assay <- CreateChromatinAssay(counts_peaks_multiome[[sample]],
                                    fragments = seurat_multiome_samples[[sample]]@assays$ATAC@fragments)
  srt <- CreateSeuratObject(chr_assay,
                            project="hRO_multiome_ATAC",
                            assay = "peaks"
                            meta.data = seurat_multiome_samples[[sample]]@meta.data)
  Annotation(srt) <- seurat_multiome_samples[[1]]@assays$peaks@annotation
  return(srt)
}), names(counts_peaks_multiome))
seurat_multiome_atac <- merge(seurat_multiome_atac[[1]], seurat_multiome_atac[[2]])
seurat_multiome_atac$cell_ident <- paste0(sapply(strsplit(colnames(seurat_multiome_atac),"-"), "[", 1), "_", seurat_multiome_atac$sample)

saveRDS(seurat_muiltiome_atac, file="data/multiome/seurat_objs/atac_only.rds")

## embed the new matrix into RNA (only retinal cells)
seurat <- readRDS("data/scRNAseq/seurat_objs/onlyRetina_processed.rds")
counts_atac <- cbind(seurat_atac_aggr@assays$peaks@counts, seurat_multiome_atac@assays$peaks@counts)
colnames(counts_atac) <- c(seurat_atac_aggr$cell_ident, seurat_multiome_atac$cell_ident)

idx_rna <- setNames(1:ncol(seurat_rna), seurat_rna$cell_ident)
idx_atac <- setNames(1:ncol(counts_atac), colnames(counts_atac))
mat_rna_atac_pairs <- sparseMatrix(i = idx_atac[rna_atac_pairs$ATAC_cell_ident],
                                   j = idx_rna[rna_atac_pairs$RNA_cell_ident],
                                   dims = c(ncol(counts_atac), ncol(seurat_rna)),
                                   dimnames = list(colnames(counts_atac), colnames(seurat_rna)))

counts_atac_rna <- counts_atac %*% mat_rna_atac_pairs
assay_atac_rna <- CreateChromatinAssay(counts_atac_rna, min.cells = 1, min.features = -1)
seurat_rna[['ATAC']] <- assay_atac_rna
seurat_rna <- RunTFIDF(seurat, assay = "ATAC")
saveRDS(seurat_rna, file="data/bimodal/onlyRetina_integrated.rds")
