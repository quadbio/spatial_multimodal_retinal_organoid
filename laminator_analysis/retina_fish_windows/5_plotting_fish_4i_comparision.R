library(tidyverse)
library(patchwork)
library(scales)
library(Seurat)
library(ggnewscale)
library(wesanderson)

files <- list.files('/links/groups/treutlein/DATA/imaging/charmel/seq_integration_fish_reimaged', full.names = TRUE)
names(files) <- list.files('/links/groups/treutlein/DATA/imaging/charmel/seq_integration_fish_reimaged')
list_fish_nuclei <- map(files, readRDS)
names(list_fish_nuclei) <- c('week13', 'week32')

df_transcript_abundance <- map(list_fish_nuclei, function(x){
  x$seu_fish@assays$RNA@data %>% rowSums(na.rm = FALSE) %>% enframe(name = 'feature', value = 'count')
}) %>% bind_rows(.id='age') %>% mutate(log_count=log10(count+1))



timepoints_all <- c("#8966A9","#707AB9","#578FC9","#3FA4D9","#3EAACC","#3EB1C0","#3EB8B4","#64BB84","#8BBE54","#B2C224","#C9C428","#E0C82C","#F8CB31",
"#F7BD2F","#F6AF2D","#F6A22B","#F28C29","#EF7727","#EC6325","#E6542B","#E14631","#DC3838")
names(timepoints_all) <- str_c('week',c(6, 7, 8, 10, 11, 12, 13, 15, 16, 18, 21, 24, 25,
28, 30, 31, 32, 34, 36, 38,  40, 46))

col_timepoint <- c(timepoints_all['week6'],
                   timepoints_all['week12'],
                   timepoints_all['week18'],
                   timepoints_all['week24'],
                   timepoints_all['week36'],
                   timepoints_all['week46']) %>% set_names(c('week6','week12','week18','week24','week39','adult'))

results <- readRDS('/local2/USERS/charmel/rotator_analysis_results_fish/matched_trajectory_correlation_reimaged_2.rds')
df_4i_windows <- results$df_4i_windows
results$df_4i_windows <- NULL

df_conf <- map(results, function(x){x$df_fish_windows %>%distinct(id) %>% nrow()}) %>%
  bind_rows() %>% pivot_longer(values_to='n', names_to='age', cols=colnames(.)) %>%
  rowwise() %>%
  mutate(q=0.01,
         conf=qt(q, n-2) / sqrt(n - 2 + qt(q, n-2)^2)) %>% ungroup()

df_plot_fish_4i_rna <- map(results,magrittr::extract,'df_cor_fish_4i_rna') %>%
  unlist(recursive = 0) %>%
  bind_rows(.id='age') %>%
  mutate(age = str_remove(age,'.df_cor_fish_4i_rna'))

p1 <- ggplot(NULL) +
  geom_bar(data=df_plot_fish_4i_rna,aes(y=reorder(feature, cor), x=cor, fill=type, group=type),
           stat='identity', position = 'dodge') +
  ggtitle('Correlation FISH - predicted scRNA-seq / 4i transcriptome')+
  xlab('Correlation') +
  ylab('') +
  geom_vline(data=df_conf, aes(xintercept=conf),linetype="dotted") +
  geom_vline(data=df_conf, aes(xintercept=abs(conf)),linetype="dotted") +
  scale_fill_manual(values = c('grey70','black') %>% set_names(unique(df_plot_fish_4i_rna$type))) +
  facet_wrap(~age) +
  new_scale_fill() +
  geom_tile(data=df_transcript_abundance %>% filter(feature %in% unique(df_plot_fish_4i_rna$feature)),
            aes(y=feature,x=0.8, width=0.025, fill=log_count), color='black') +
  scale_fill_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys") +
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6), limits=c(NA,0.85)) +
  theme_minimal()

ggsave(plot=p1,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/correlation_matched_nuclei_by_timepoint.pdf'),
       dpi = 250, unit='mm', width = 400.05, height = 215.22*1.5)

df_plot_scatter <- df_plot_fish_4i_rna %>%
  pivot_wider(id_cols = c('age','feature'), values_from = cor, names_from = type) %>%
  left_join(df_transcript_abundance, by = c('age','feature')) %>%
  mutate(age = factor(age, levels = c('week13','week32'), labels = c('Week 13','Week 32')),
         feature = ifelse(feature %in% c("PDE6B","VSX2","RBP3","SLC17A7","VSX2","NR2F1","VEGFA","ONECUT1","OTX2","PPIB","STAT3"),feature, ''))
library(ggrepel)
p <- ggplot(df_plot_scatter, aes(x=`predicted_RNA_4i_nuclei-fish_rna`, y=`predicted_RNA_fish_nuclei-fish_rna`, label=feature)) +
  geom_smooth(method = 'lm', col='grey10', size=0.5) +
  geom_point(aes(col=log_count), size=1.5, alpha=1) +
  facet_wrap(~age, scales='free') +
  scale_x_continuous()+
  scale_y_continuous()+
  theme_classic() +
  theme(aspect.ratio = 1,
        strip.background = element_blank()) +
  geom_vline(data = df_conf %>% mutate(age = factor(age, levels = c('week13','week32'), labels = c('Week 13','Week 32'))), aes(xintercept=abs(conf)), linetype='dotted') +
  geom_hline(data = df_conf %>% mutate(age = factor(age, levels = c('week13','week32'), labels = c('Week 13','Week 32'))), aes(yintercept=abs(conf)), linetype='dotted') +

  scale_color_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys") +
  xlim(c(NA,0.75)) +
  ylim(c(NA,0.75)) +
  xlab('Imputed RNA 4i - RNA FISH') +
  ylab('Imputed RNA FISH - RNA FISH') +
  ggtitle('Correlation matched 4i - FISH nuclei') +
  geom_label_repel( box.padding = 0.75, max.overlaps = Inf, size=1.5)#+
  #guides(color=guide_legend('log10(x+1)'))

ggsave(plot=p,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/scatter_matched_nuclei_by_timepoint.pdf'),
       dpi = 250, unit='mm', width = 400.05/2, height = 215.22/2)



df_plot_fish_4i_protein_rna <- map(results,magrittr::extract,'df_cor_fish_4i_protein_rna') %>%
  unlist(recursive = 0) %>%
  bind_rows(.id='age') %>%
  mutate(age = str_remove(age,'.df_cor_fish_4i_protein_rna'))

types <-  c('predicted_RNA_4i_nuclei-fish_rna','predicted_RNA_4i_nuclei-protein','predicted_RNA_fish_nuclei-fish_rna','protein_4i_nuclei-fish_rna')
label_types <-  c('RNA 4i (Predicted) - RNA FISH (Measured)',
                                                      'RNA 4i (Predicted) - Protein 4i (Measured)',
                                                      'RNA FISH (Predicted) - RNA FISH (Measured)',
                                                      'Protein 4i (Measured) - RNA FISH (Measured)')
p <- df_plot_fish_4i_protein_rna %>% filter(feature =='VSX2') %>%
  mutate(type=factor(type, levels = types, labels = label_types)) %>%
  ggplot(aes(x=age,group=type, y=cor, fill=type)) +
  geom_bar(position = 'dodge', stat = 'identity', col='grey90') +
  theme_classic() +
  xlab('') +
  ylab('Correlation') +
  scale_fill_manual(values = c('black','grey20','grey50','grey70') %>% set_names(label_types)) +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=4,byrow=TRUE))

ggsave(plot=p,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/barplot_cor_vsx2_by_timepoint.pdf'),
       dpi = 250, unit='mm', width = 400.05/2, height = 215.22)

p2 <- ggplot(NULL) +
  geom_bar(data=df_plot_fish_4i_protein_rna %>%
  mutate(type=factor(type, levels = types, labels = label_types)),
           stat = 'identity',
           position = 'dodge',
           aes(y = reorder(feature, cor), x = cor, group = type, fill = type)) +
  xlab('Correlation') +
  ylab('') +
  ggtitle('Correlation FISH - predicted transcriptomes - 4i protein intensity') +
  facet_wrap(~age) +
  geom_vline(data=df_conf, aes(xintercept=conf),linetype="dotted") +
  geom_vline(data=df_conf, aes(xintercept=abs(conf)),linetype="dotted") +
  scale_fill_manual(values = c('black','grey20','grey50','grey70') %>% set_names(label_types)) +
  new_scale_fill() +
  geom_tile(data=df_transcript_abundance %>% filter(feature %in% unique(df_plot_fish_4i_protein_rna$feature)),
            aes(y=feature,x=0.8, width=0.025, fill=log_count), col='black') +
  scale_fill_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys") +
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6), limits=c(NA,0.85)) +
  theme_minimal()

ggsave(plot=p2,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/correlation_matched_nuclei_protein_overlap_by_timepoint.pdf'),
       dpi = 250, unit='mm', width = 400.05, height = 215.22)

df <- df_plot_fish_4i_rna %>% left_join(df_transcript_abundance, by=c('age','feature'))

ggplot(df, aes(x=log_count, y=cor, col=type)) +
  geom_point() +
  facet_wrap(~age) +
  coord_fixed(ratio = max(df$log_count)/max(df$cor)) +
  geom_hline(data=df_conf, aes(yintercept=abs(conf)),linetype="dotted") +
  theme_linedraw() +
  scale_color_manual(values = c('grey70','black') %>% set_names(unique(df_plot_fish_4i_rna$type))) +
  theme(aspect.ratio = 1) +
  geom_smooth(method = 'lm')






df_nuclei_matched <- map(results, function(x){
  map(x$matched_trajectory, magrittr::extract, 'df_matched_nuclei') %>% unlist(recursive = FALSE) %>%
  bind_rows()
}) %>% bind_rows(.id='age') %>%
  left_join(df_4i_windows %>%
              select(dpt_rank, age) %>%
              rename(age_4i = age), by = 'dpt_rank') %>%
  mutate(age_4i=factor(age_4i, levels=c('week6','week12','week18','week24','week39','adult')))

p3 <- ggplot(df_nuclei_matched %>% distinct(fish_window,dpt_rank, age),aes(x=dpt_rank, group=age, fill=age)) +
  geom_histogram(position='dodge') +
  theme_bw() +
  scale_fill_manual(values=timepoints_all)
p4 <- ggplot(df_nuclei_matched %>% distinct(fish_window,mean_cor, max_cor, age, age_4i),aes(x=age_4i, y=max_cor, fill=age_4i)) + geom_boxplot() + facet_wrap(~age) +
  scale_fill_manual(values=col_timepoint)

df_cor_trajectory <- map(results, function(x) {
  map(x$matched_trajectory, function(y) {
    map(y$all_max_cors, median) %>%
      enframe() %>%
      unnest(value) %>%
      rename(dpt = name, mean_cor = value) %>%
      mutate(dpt = as.numeric(dpt))
  }) %>% bind_rows(.id = 'window_fish')
}) %>%
  bind_rows(.id = 'age') %>%
  mutate(dpt_rank = as.numeric(dpt)) %>%
  left_join(df_4i_windows %>%
              select(dpt_rank, age) %>%
              rename(age_4i = age), by = 'dpt_rank')

p5 <- ggplot(df_cor_trajectory %>%
               group_by(age) %>%
               mutate(similarity_score=rescale(mean_cor)) %>%
               group_by(dpt_rank, age) %>%
               summarise(similarity_score=mean(similarity_score)), aes(x=dpt_rank, y=similarity_score, group=age, col=age, fill=age)) +
  geom_smooth() +
  scale_color_manual(values=timepoints_all)

p6 <- ggplot(df_cor_trajectory, aes(x=age_4i, y=mean_cor, fill=age_4i)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~age) +
  scale_fill_manual(values=col_timepoint)


df_cor_trajectory %>% group_by(age_4i, age) %>% summarise(median=median(mean_cor), mad=mad(mean_cor)) -> df_summary


ggplot(df_summary, aes(x=str_remove(age_4i,'week') %>% as.numeric(), y=median, group=age, col = age)) + geom_line() + geom_point()



p1 <- results$week13$matched_trajectory %>%
  enframe() %>%
  mutate(test=map(value,~.x$all_dists)) %>%
  mutate(test=map(test,unlist), test=map(test,function(x){enframe(x)%>%mutate(name=as.numeric(name))})) %>%
  select(name, test) %>%rename(window=name) %>% unnest(test) %>%
  rename(dpt_rank = name) %>%
  left_join(df_4i_windows %>% select(dpt_rank, age)) %>%
  mutate(age=factor(age, levels=c('week6','week12','week18','week24','week39','adult'))) %>%
  ggplot(aes(x=age, y=value, fill=age)) +
  geom_boxplot() +
  ylab('Distance') +
  xlab(NULL) +
  scale_fill_manual(values=col_timepoint) + ggtitle('Week 13') +
  scale_y_reverse() +
  theme_bw()

p2 <- results$week32$matched_trajectory %>%
  enframe() %>%
  mutate(test=map(value,~.x$all_dists)) %>%
  mutate(test=map(test,unlist), test=map(test,function(x){enframe(x)%>%mutate(name=as.numeric(name))})) %>%
  select(name, test) %>%rename(window=name) %>% unnest(test) %>%
  rename(dpt_rank = name) %>%
  left_join(df_4i_windows %>% select(dpt_rank, age)) %>%
  mutate(age=factor(age, levels=c('week6','week12','week18','week24','week39','adult'))) %>%
  ggplot(aes(x=age, y=value, fill=age)) +
  geom_boxplot(outlier.shape = NA) +
  ylab('Distance') +
  xlab(NULL) +
  scale_fill_manual(values=col_timepoint) +
  ggtitle('Week 32') +
  scale_y_continuous(limits = quantile(dfr$y, c(0.1, 0.9)))
  scale_y_reverse() +
  theme_bw()

p1 | p2



df_plot <- map(results, function (x){x$matched_trajectory %>%
  enframe() %>%
  mutate(test=map(value,~.x$all_dists)) %>%
  mutate(test=map(test,unlist), test=map(test,function(x){enframe(x)%>%mutate(name=as.numeric(name))})) %>%
  select(name, test) %>%rename(window=name) %>% unnest(test) %>%
  rename(dpt_rank = name) %>%
  left_join(df_4i_windows %>% select(dpt_rank, age)) %>%
  mutate(age=factor(age, levels=c('week6','week12','week18','week24','week39','adult')))}) %>% bind_rows(.id = 'age_fish')

p <- ggplot(df_plot,aes(x=age, y=value, fill=age)) +
  geom_boxplot(outlier.shape = NA) +
  ylab('Distance') +
  xlab(NULL) +
  scale_fill_manual(values=col_timepoint) +
  facet_wrap(~age_fish) +
  scale_y_reverse() +
  theme_bw()

ggsave(plot=p,
       filename = paste0('/links/groups/treutlein/DATA/imaging/charmel/plots/distances_laminar_window_fish_4i_by_timepoint.pdf'),
       dpi = 250, unit='mm', width = 400.05, height = 215.22)