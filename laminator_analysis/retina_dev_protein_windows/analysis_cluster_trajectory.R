# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 01.06.22

library(tidyverse)
library(scales)
library(furrr)
library(Matrix)
library(Seurat)
library(presto)
library(RColorBrewer)
library(paletteer)
library(patchwork)
library(grid)
library(ggplotify)

list_results <- readRDS('/local2/USERS/charmel/spatial_densities.rds')
list_results$seu <- list_results$seu %>% FindClusters(resolution = 0.35)
saveRDS(list_results, '/local2/USERS/charmel/spatial_densities_corrected.rds')
# plot cluster zones
qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
set.seed(42)
cluster_cols <- sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))),list_results$seu@meta.data %>% distinct(seurat_clusters)%>%nrow())
names(cluster_cols) <- list_results$seu@meta.data %>% distinct(seurat_clusters) %>% pull() %>% sort()
p1 <- ggplot(list_results$seu@meta.data, aes(x=x, y=y, fill=seurat_clusters)) + geom_raster() +
  theme_void() + coord_equal() + scale_fill_manual(values=cluster_cols) + scale_y_reverse() +
  geom_text(data = list_results$seu@meta.data  %>% group_by(seurat_clusters) %>% summarise(x=mean(x), y= mean(y)), aes(x=x, y=y, label=seurat_clusters))+
  theme(legend.position = 'none') +
  ggtitle('Pseudotemporal-spatial cluster zones')

# find cluster markers
df_de <- wilcoxauc(list_results$seu, group_by = 'seurat_clusters', slot='scale.data')
df_de <- AverageExpression(list_results$seu, group.by = 'seurat_clusters', assay='RNA', slot='scale.data')$RNA %>% wilcoxauc(colnames(.)) %>% as_tibble()
# correlation gene patterns vs gene patterns
list_results$seu@assays$RNA@scale.data %>% t() %>% cor() -> cor_results
de_genes<- df_de  %>% group_by(group) %>% arrange(-logFC) %>% slice(1:20) %>% distinct(feature) %>% pull(feature) %>% unique()

# heat maps of correlated genes
p2 <- as.ggplot(~heatmap(cor_results[de_genes,de_genes], scale='none', symm = 1, col=brewer.pal(11,"RdBu") %>% rev())) +
  ggtitle('Correlation transcript densities - top 20 logFC of each cluster')

# heat map for cluster vs genes
p3 <- as.ggplot(~AverageExpression(list_results$seu,
                  group.by = 'seurat_clusters',
                  assay='RNA', features = de_genes,
                  slot='scale.data')$RNA %>%
   t() %>% heatmap(scale='col',RowSideColors = cluster_cols, col=brewer.pal(11,"RdBu") %>% rev())) +
  coord_equal() +
  ggtitle('Average transcript density by cluster - top 20 logFC of each cluster')


# heat map celltypes
p4 <- as.ggplot(~AverageExpression(list_results$seu, group.by = 'seurat_clusters', assay='CT')$CT %>%
   t() %>% heatmap(scale='col',RowSideColors = cluster_cols,
                   col=brewer.pal(11,"RdBu") %>% rev())) +
  coord_equal() +
  ggtitle('Average celltype densities by cluster')

((p1 )/p4)|p3

p2


df <- list_results$seu@meta.data %>% mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

n <- max(df$x)
m <- max(df$y)
M <- matrix(0, nrow=m, ncol=n)
for (i in 1:nrow(df)){
    M[df[i,'y'], df[i,'x']] = df[i,"seurat_clusters"]
}

library(raster)
r <- raster(M)

get_adjacent_clusters <- function(cell,r){
  cell_cluster <- extract(r, cell)
  cells <- adjacent(r, cell, 8, pairs=FALSE)
  adj_cluster <- cbind(xyFromCell(r, cells), cluster_adj=extract(r, cells)) %>%
    as_tibble() %>% distinct(cluster_adj) %>%
    mutate(cluster_cell = cell_cluster)
  return(adj_cluster)
}
df_adj <- map(1:30900, get_adjacent_clusters, r=r) %>% bind_rows()


library(igraph)
M_adj <- df_adj %>% distinct() %>% arrange(cluster_adj, cluster_cell) %>% graph_from_data_frame() %>% as_adjacency_matrix() %>% as.matrix()
saveRDS(M_adj, '/local2/USERS/charmel/adj_matrix_cluster_zones.rds')