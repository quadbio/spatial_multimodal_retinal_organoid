import os
from pathlib import Path
import numpy as np
from laminator_analysis import laminator
import pandas as pd
from time import gmtime, strftime
import re
from tqdm import tqdm
from skimage import io
from skimage.color import rgba2rgb, rgb2gray

dir_images = 'data/processed/4i/MTU_results'
dir_results = 'data/processed/4i/laminator/results_mtu'
dir_masks = 'data/processed/4i/masks'
masks = os.listdir(dir_masks)
points = os.listdir(dir_images)

for point in points:
    print('Started Laminator...')
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Processing sample: ' + str(point))

    dir_results_point = Path(dir_results, str(point))
    Path(dir_results_point).mkdir(parents=True, exist_ok=True)

    # set directories and get paths for masks and images

    path_mask = Path(dir_masks, masks[[re.findall(r'\d+', path)[0] for path in masks].index(point)])
    dir_images_point = Path(dir_images, str(point), 'single_MTU_imgs')
    hoechst = io.imread(Path('data/raw/4i/images',point,'channel_hoechst.tif'))
    image_paths = os.listdir(dir_images_point)

    pd.DataFrame({'file':[str(path) for path in image_paths]}).to_csv(Path(dir_results_point, 'stain_paths.csv'))

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Loading images...')
    # load mask
    mask = laminator.load_mask(path_mask)

    # load images
    stack = []
    for path in tqdm(image_paths):
        img = io.imread(path)
        img = rgb2gray(rgba2rgb(img, background=(0,0,0)))
        img = np.where(img > 0, 1, 0)
        stack.append(img)
    stack = np.dstack(stack)

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' - Calculating contour and orienting windows...')
    # get contours
    contours = laminator.get_contour(mask)

    # calculate distance transform of mask
    dist_transform = laminator.distance_transform(mask)

    # load angles and contour
    df_meta = pd.read_csv(Path('/nas/groups/treutlein/DATA/imaging/charmel/laminator_analysis_wo_scaling_2',
                               point,
                               'df_meta.csv'))

    # assess contours and angles
    laminator.assess_oriented_windows(df_meta, laminator.scale_image(hoechst),
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

