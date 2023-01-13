library(tidyverse)

# prepare data structure for baysor bash-script
dir_results <- 'data/raw/fish/spot_data'
dir_output <- 'data/processed/fish/baysor'
dir.create(dir_output)

# move and rename spot data
files_spots <- list.files(dir_results, pattern = '.txt', full.names = TRUE)
names(files_spots) <- list.files(dir_results, pattern = '.txt')
df_spots <- map(files_spots,read_tsv, col_names=FALSE) %>%
  bind_rows(.id='ID') %>%
  select(-X5) %>%
  rename(x=X1, y=X2, z=X3, transcript=X4) %>%
  mutate(ID = str_remove(ID,'32955-slide928-2_') %>% str_remove('_results.txt'),
         z= z * (300/138)) %>%
  filter(!ID %in% c('D2-1','D2-2')) %>% rename(organoid = ID)

for (i in unique(df_spots$organoid)){
   dir <- paste0('data/processed/fish/baysor/', i)
   dir.create(dir)
   write_csv(df_spots %>%
               filter(organoid == i) %>%
               rename(gene=transcript), paste0(dir,'/','spots.csv'))
 }

# move segmentation label images
dir <- paste0('data/processed/fish/segmented')
files <- list.files(dir, full.names = TRUE)
names(files) <- list.files(dir)
for (file in seq_along(files)){
  sample <- str_remove(names(files[file]),'32955-slide928-2_') %>% str_remove('_DAPI.tiff')
  file.copy(files[file], paste0(dir_output,'/',sample,'/segmentation.tiff'))
}




