library(reticulate)
use_condaenv("conda_3.7.9", required=TRUE)
Sys.unsetenv("LD_LIBRARY_PATH")

library(FlowSOM)
library('SamSPECTRAL')
library(dplyr)
library(yaml)


np <- import("numpy")


if (!file.exists(paste0(data.path, '/fSOM_output/fSOM.rds'))){
pixel_matrix <- np$load(paste0(data.path, 'pixel_matrices/scaled_fusion_sub250_z_normalised_pixel_matrix.npz'))
pixel_matrix = pixel_matrix$f[["arr_0"]]

colnames = read_yaml('metadata.yml')
params = read_yaml('params.yml')

colnames(pixel_matrix) = colnames

pixel_matrix = pixel_matrix[, -which(colnames(pixel_matrix) == 'excluded')]

ydim_SOM = 30
xdim_SOM = 30
dist_metric = 2
num_runs = 10
marker_cols = colnames(pixel_matrix)
data_FlowSOM <- flowCore::flowFrame(pixel_matrix)

set.seed(1234)

# run FlowSOM
fSOM = FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
fSOM = FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols, ydim = ydim_SOM, xdim = xdim_SOM, distf = dist_metric, rlen = num_runs)

dir.create(paste0(data.path, '/fSOM_output/'), recursive = TRUE)
saveRDS(fSOM, file = paste0(params$data_path, '/fSOM_output/fSOM.rds'))
som_codes = fSOM$map$codes
np$savez_compressed(paste0(params$data_path, '/fSOM_output/som_codes.npz'), som_codes)
} else {
  fSOM = readRDS(paste0(params$data_path, '/fSOM_output/fSOM.rds'))
  som_codes = fSOM$map$codes
}

library(Rphenograph)
n_clusters_detected = list()
neigh_to_check = c(2:50)
for (i in neigh_to_check) {
  phgr_node_labels <- Rphenograph(fSOM$map$codes, k = i)[[2]] %>% membership()
  n_clusters_detected[i] <- length(unique(phgr_node_labels))
}


library(Rphenograph)
n_clusters_detected = list()
neigh_to_check = c(2:50)
for (i in neigh_to_check) {
  phgr_node_labels <- Rphenograph(fSOM$map$codes, k = i)[[2]] %>% membership()
  n_clusters_detected[i] <- length(unique(phgr_node_labels))
}

detected <- kneepointDetection(vect = rev(unlist(n_clusters_detected)), PlotFlag = FALSE)
print(detected$MinIndex)
df = data.frame(neighs = rev(neigh_to_check), clusts = rev(unlist(n_clusters_detected)))
unlist(n_clusters_detected)[15]

png(file=paste0(params$data_path, '/fSOM_output/kneeplot.png'), width=600, height=350)
plot(rev(unlist(n_clusters_detected)), rev(neigh_to_check), xaxp = c(0, 200, 20))
abline(v = df[detected$MinIndex, 2], col = "blue", lwd = c(1, 3))
abline(h = df[detected$MinIndex, 1], col = 'red', lwd = c(1, 3))
dev.off()

print(df)
print(paste0('Choose k nearest neighbors to use for phenograph fom kneeplot in ', params$data_path, '/fSOM_output. Or choose from table printed above, or use default of 15 neighbors to get 31 MTUs'))






















































