import sys
sys.path.extend(['/home/charmel/4i_organoid_pipeline/retina'])
import os
from pathlib import Path
import numpy as np
from laminator_analysis import laminator
import pandas as pd
from time import gmtime, strftime
from skimage import io
from tqdm import tqdm


dir_masks = '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission/masks_from_spots'
dir_output = '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission/laminator_results'
files = [file for file in os.listdir(dir_masks) if 'DAPI.tiff' in file]
dir_tables = '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission'
files_tables = [file for file in os.listdir(dir_tables) if not file.startswith('.') and 'results.txt' in file]

for file in files:
    print('Started laminator...')
    path_mask = Path(dir_masks, file)
    point = file.replace('32955-slide928-2_','').replace('_DAPI.tiff','')
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Processing sample: ' + str(point))
    if 'A' in point or 'B' in point:
        print('Skipped sample.')
        continue
    file_table = [file for file in files_tables if point in file][0]
    df = pd.read_csv(Path(dir_tables,file_table),sep='\t', names=['x','y','z','transcript'], index_col=False)

    # load mask
    mask = laminator.load_mask(path_mask)
    if np.sum(mask) == 0:
        print('Empty mask detected. Skipped sample.')
        continue
    dir_results_point = Path(dir_output, point)
    #if dir_results_point.exists():
    #    print('Already processed. Skipped sample.')
    #    continue
    dir_results_point.mkdir(parents=True, exist_ok=True)
    img_hoechst = io.imread(Path(dir_tables,file))
    # generate binary array stacks from spot dataframe
    stack = []
    transcripts = df['transcript'].unique().tolist()
    dict_transcripts = dict((i,j) for i,j in enumerate(transcripts))
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Constructing spot arrays...')
    stack = np.zeros((mask.shape[0], mask.shape[1],len(transcripts)),dtype='uint8')
    for i,transcript in dict_transcripts.items():
        df_tmp = df[df['transcript'] == transcript]
        for x,y in zip(df_tmp['x'].values,df_tmp['y'].values):
            stack[y,x,i] = 1

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Calculating contour and orienting windows...')
    # get contours
    contours = laminator.get_contour(mask, smoothing=False)

    # calculate distance transform of mask
    dist_transform = laminator.distance_transform(mask, smoothing=False)

    # calculate angles
    df_meta = laminator.calculate_angles(contours,
                                       distance_transform=dist_transform,
                                       window_size=1000,
                                       slice_width=50)

    # assess contours and angles
    laminator.assess_oriented_windows(df_meta, laminator.scale_image(img_hoechst),
                                    show_angles=True,
                                    save_plot=True,
                                    show_plot=False,
                                    file=str(Path(dir_results_point, 'contour_angles')))

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Orienting and retrieving intensity profiles...')
    # get intensity profiles from rotated slices
    df_intensity_profile = laminator.rotate_slices(df_meta, stack)
    df_intensity_profile = df_intensity_profile.replace({"stain": dict_transcripts})
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Saving results...')
    df_meta.to_csv(Path(dir_results_point, 'df_meta.csv'), index=False)
    for i,transcript in tqdm(dict_transcripts.items(), ascii=True):
        df_intensity_profile[df_intensity_profile['stain']==transcript].to_csv(Path(dir_results_point, 'df_intensity_profiles_{}.csv'.format(transcript)), index=False)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - laminator analysis for sample ' + str(point) + ' completed.')
    print()
