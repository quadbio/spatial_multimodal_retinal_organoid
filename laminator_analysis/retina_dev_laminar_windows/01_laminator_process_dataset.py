import os
from pathlib import Path
import numpy as np
from laminator_analysis import laminator
import pandas as pd
from time import gmtime, strftime
import re

dir_images = '/links/groups/treutlein/DATA/imaging/charmel/shiny_input'
points = os.listdir('/links/groups/treutlein/DATA/imaging/charmel/shiny_input')
dir_results = '/links/groups/treutlein/DATA/imaging/charmel/laminator_analysis_wo_scaling_2'
dir_masks = '/links/groups/treutlein/DATA/imaging/charmel/refined_masks'
masks = os.listdir(dir_masks)
masks = [mask for mask in masks if '_pw_' in mask]



for point in points:
    print('Started laminator...')
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Processing sample: ' + str(point))

    dir_results_point = Path(dir_results, str(point))
    Path(dir_results_point).mkdir(parents=True, exist_ok=True)

    # set directories and get paths for masks and images
    path_mask = Path(dir_masks, masks[[re.findall(r'\d+', path)[0] for path in masks].index(point)])
    dir_images_point = Path(dir_images, str(point))
    image_paths = os.listdir(dir_images_point)
    image_paths = [Path(dir_images_point, path) for path in image_paths if 'channel' in path]
    index_hoechst = [str(path) for path in image_paths].index(str(Path(dir_images_point, 'channel_hoechst.tif')))
    pd.DataFrame({'file':[str(path) for path in image_paths]}).reset_index().to_csv(Path(dir_results_point, 'stain_paths.csv'), index=False)

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Loading images...')

    # load mask
    mask = laminator.load_mask(path_mask)
    # load images
    stack = laminator.load_images(image_paths)

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Calculating contour and orienting windows...')
    # get contours
    contours = laminator.get_contour(mask)

    # calculate distance transform of mask
    dist_transform = laminator.distance_transform(mask)

    # calculate angles
    df_meta = laminator.calculate_angles(contours,
                                       distance_transform=dist_transform,
                                       window_size=1000,
                                       slice_width=50)

    # assess contours and angles
    laminator.assess_oriented_windows(df_meta, laminator.scale_image(stack[...,index_hoechst]),
                                    show_angles=True,
                                    save_plot=True,
                                    show_plot=False,
                                    file=str(Path(dir_results_point, 'contour_angles')))

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime())+' - Orienting and retrieving intensity profiles...')
    # get intensity profiles from rotated slices
    df_intensity_profile, stack_slices = laminator.rotate_slices(df_meta, stack)

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime())+' - Saving results...')
    df_meta.to_csv(Path(dir_results_point, 'df_meta.csv'), index=False)
    df_intensity_profile.to_csv(Path(dir_results_point, 'df_intensity_profiles.csv'), index=False)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime())+' - laminator analysis for sample ' + str(point) + ' completed.')
    print()

