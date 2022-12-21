source('~/4i_organoid_pipeline/retina/rotator_analysis/utils.R')
dir_results <- '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission'
dir_rotator <- str_c(dir_results,'/rotator_results')



z_scoring_age_sample <- function(df){
    message("Z-scoring across stain, age and sample...")
    df <- df %>% mutate(age=ifelse(str_detect(organoid,'A|B'),'week32','week13'))
    df_age <- df %>% group_by(age, stain) %>% summarize(mean_age = mean(intensity),
                                             sd_age = sd(intensity))
    df_sample <- df %>% group_by(stain, organoid) %>% summarize(mean_sample = mean(intensity),
                                             sd_sample = sd(intensity))
    message("Merging datasets and calculating intensities...")
    df <- df  %>% left_join(df_age, by=c('stain','age')) %>%
      left_join(df_sample, by=c('stain','organoid')) %>%
      mutate(intensity_scaled = ((intensity - mean_sample)/sd_sample) * sd_age + mean_age)

    message('Scaling data...')
    df <- df %>% group_by(stain) %>% mutate(intensity = rescale(intensity_scaled, to = c(0,1))) %>%
      ungroup() %>%
      mutate(intensity=ifelse(is.nan(intensity),0,intensity)) %>%
      select(id,organoid,radial_position,stain,window,intensity)
    return(df)
}



read_and_process_rotator_results <- function(dir,samples=NULL,dir_output='/local2/USERS/charmel/rotator_dataset_fish_2/'){
    directories <- list.files(dir, full.names = 1)
    names(directories) <- list.files(dir)
    if(!is.null(samples)){
        directories <- directories[samples]
    }
    directories <- directories[!str_detect(directories,'df_contours')]
    file_paths <- map(directories, list.files, pattern = ".csv", full.names=1)
    files <- map(file_paths, function(x){x[!str_detect(x,'df_meta.csv')]})
    message('Reading in data...')
    df_intensity <- map(files,function(x){
       map(x, read_csv) %>% bind_rows()
    }) %>% bind_rows(.id='organoid') %>%
      unite('id',c('organoid','window'),remove=FALSE) %>%
      z_scoring_age_sample()
    message('Getting meta data and detecting cut areas...')
    df_meta <- map(directories, function(x){read_csv(str_c(x,'df_meta.csv',sep='/'))}) %>%
      bind_rows(.id = 'organoid') %>% filter(!organoid %in% c('D2-1','D2-2')) %>%
      unite('id',c('organoid','window'),remove=FALSE)
    message('Filtering intensity profiles')
    ids_selected <- df_meta %>% pull(id)

    df_intensity <- df_intensity %>% filter(id %in% ids_selected)

    shared_transcripts <- df_intensity %>% distinct(stain,organoid) %>% ungroup() %>%
      mutate(n_organoids=length(unique(organoid))) %>% group_by(stain) %>%
      mutate(n_stain=n()) %>% filter(n_stain==n_organoids) %>% distinct(stain) %>%
      pull(stain)

    for (i in seq_along(shared_transcripts)){
      message('Saving intensity profiles for: ', shared_transcripts[i])
      df_intensity %>% filter(stain == shared_transcripts[i]) %>% write.csv(paste0(dir_output,shared_transcripts[i],'.csv'))
    }

    message('Saving meta files...')
    df_meta %>% saveRDS(paste0(dir_output,'df_meta.rds'))
    message('Done.')
}

read_and_process_rotator_results(dir_rotator)