# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 02.02.22

# Title     : TODO
# Objective : TODO
# Created by: harmelc
# Created on: 11.10.21

setwd('/home/charmel/4i_organoid_pipeline/retina/Rotator/Analysis')
source('../utils.R')

# load 4i nuclei datasets
list_4i_nuclei <- readRDS('/local2/USERS/charmel/integrated_seurat_4i_scRNA_scATAC_v3.rds')
names(list_4i_nuclei) <- c('adult','week12','week18','week24','week39','week6')

df_subset <- map(list_4i_nuclei, function (x){x$seu_4i@meta.data %>% mutate(organoid=as.character(organoid)) %>% as_tibble()}) %>%
  bind_rows() %>% filter(age=='week12', cell_type=='RGC') %>% mutate(x=round(X), y=round(Y))
subset <- df_subset %>% distinct(organoid) %>% pull()


df_intensity <- read_radial_cluster_profiles(subset=subset, df_subset_positions=df_subset) %>%
  scale_cluster_profiles()

message('Filtering intensity profiles for selected windows...')
df_intensity <- df_intensity %>%
  ungroup() %>%
  mutate(stain=str_c("cluster_",stain)) %>%
  group_by(stain) %>%
  mutate(intensity=rescale(intensity)) %>% ungroup()



df_expression <- df_intensity %>%  left_join(df_profiles %>% select(x,y,louvain), by=c('x','y')) %>% group_by(louvain,stain) %>%
  summarise(intensity=mean(intensity)) %>% pivot_wider(names_from = 'stain', values_from = 'intensity') %>% arrange(louvain) %>% ungroup()


groups <- df_expression$louvain

wilcoxauc(df_expression %>% select(-c('louvain')) %>% as.matrix() %>% t(), groups) %>% as_tibble() -> df_de_neigh
df_markers_neigh <- top_markers(df_de_neigh, 10)

df_de_plot_neigh <- df_de_neigh %>%
  left_join(df_markers_neigh %>% select(-rank) %>% pivot_longer(cols = colnames(.), names_to = 'group', values_to = 'feature') %>% mutate(marker=TRUE), by=c('group','feature')) %>%
  mutate(color=ifelse(abs(logFC)>=1,'red','grey50'),
         color_new=ifelse(!is.na(marker),'blue',color))
ggplot(df_de_plot_neigh , aes(x=logFC, y=-log10(padj), col=color_new)) + geom_point(size=.25) + facet_wrap(~group, nrow = 2, scales = 'free') + theme_light() +
  geom_vline(xintercept = c(1,-1), linetype = "dashed", col='grey70') +
  geom_text_repel(data = df_de_plot_neigh %>% filter(color_new !='grey50'),aes(label=feature), col='black', size=2, box.padding = 0.1) +
  coord_cartesian(clip = "off") +
  scale_color_identity()





df_profiles <- readRDS('/local2/USERS/charmel/dpt_heterogeneity_w12_rgc.rds') %>% bind_rows(.id='max_radius') %>% distinct() %>% filter(max_radius=='40')

df_intensity %>% filter(!(stain %in% c('cluster_15','cluster_5','cluster_27','cluster_30'))) %>% left_join(df_profiles %>% select(x,y,louvain), by=c('x','y')) %>% group_by(louvain,stain) %>%
  summarise(intensity=mean(intensity)) %>% pivot_wider(names_from = 'stain', values_from = 'intensity') %>% arrange(louvain) %>% ungroup() -> df_heatmap

df_cluster_labels <- read_tsv( '/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/labelannot_fusion_znormalised.txt')
labels_cluster <- df_cluster_labels %>% mutate(cluster=str_c('cluster',Cluster,sep = '_')) %>% select(cluster, IDandLabel) %>% deframe()
M <- df_heatmap %>%select(-louvain) %>% as.matrix()
rownames(M) <- df_heatmap$louvain
colnames(M) <- as_tibble(colnames(M)) %>% mutate(names = ifelse(str_detect(value,'cluster*'), labels_cluster[value], value)) %>% pull(names)
heatmap(M, scale='none',col=RColorBrewer::brewer.pal(11,"RdBu") %>% rev(), margins = c(20,5))




df_intensity %>% filter(stain=='cluster_13') %>% left_join(df_profiles %>% select(x,y,louvain), by=c('x','y')) %>% group_by(louvain,x,y) %>%
  summarise(intensity=mean(intensity))  %>% arrange(louvain) %>% ungroup() %>%
  left_join(df_subset %>% select(x,y,distance_to_contour), by=c('x','y')) %>%
  mutate(louvain_grouped=ifelse(louvain %in% c('1','5'),'1,5',ifelse(louvain %in% c('7','8'),'7,8','2,3,4,6'))) %>%
  group_by(louvain_grouped) %>% summarise(cor=cor(intensity, distance_to_contour))


  ggplot(aes(x=distance_to_contour, y=intensity, col=louvain)) + geom_point() + facet_wrap(~louvain_grouped)







dir <- '/local2/USERS/charmel/rotator_dataset_nuclei_mtu/'

features <- df_intensity %>% distinct(stain) %>% pull()

for (feature in seq_along(features)){
  message(paste('Writing dataset for stain:', features[feature]))
  df_intensity %>%
    filter(stain == features[feature]) %>%
    write.csv(paste0(dir, features[feature], '.csv'))
}