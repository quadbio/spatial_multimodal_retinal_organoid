import os
from pathlib import Path
from laminator_analysis import laminator
import pandas as pd
import numpy as np
from tqdm import tqdm

# set directories and get paths for masks and images
dir_masks = 'data/processed/fish/masks_from_spots'
files = [file for file in os.listdir(dir_masks) if 'DAPI.tiff' in file]

contours = []
for file in tqdm(files):
    path_mask = Path(dir_masks, file)
    point = file.replace('32955-slide928-2_', '').replace('_DAPI.tiff', '')
    # load mask
    mask = laminator.load_mask(path_mask)
    if (np.sum(mask) == 0):
        continue
    # get contours
    contours.append(laminator.get_contour(mask, smoothing=False).assign(organoid=point))

contours = pd.concat(contours)

dir_results = Path('data/processed/fish/laminator/analysis_results')
dir_results.mkdir(parents=True, exist_ok=True)
contours.to_csv(Path(dir_results, 'df_contours.csv'), index=False)
