import os
from pathlib import Path
import pandas as pd
from skimage import io
import random
import matplotlib as mpl
from multiprocessing import Pool
from tqdm import tqdm
from skimage.transform import rotate
from skimage.color import gray2rgba
import numpy as np
from skimage.color import rgba2rgb
from skimage.color import rgb2gray

def export_slices(df_meta, img, dir_results, rgba):
    positions = df_meta['position'].values
    slices = []
    for i in positions:
        df_tmp = df_meta[df_meta['position'] == i]
        id = df_tmp['id'].values[0]
        window_size = df_tmp['window_size'].values[0]
        slice_width = df_tmp['slice_width'].values[0]
        angle = df_tmp['angle'].values[0]
        x = df_tmp['x'].values[0]
        y = df_tmp['y'].values[0]
        x_offset = df_tmp['x_offset'].values[0]
        y_offset = df_tmp['y_offset'].values[0]
        img_rotated = rotate_slice(img, angle, x, y, x_offset, y_offset, window_size, slice_width, rgba=rgba)
        io.imsave(Path(dir_results, id + '.png'), img_rotated, check_contrast=False)
    return slices


def rotate_slice(img, angle, x, y, x_offset, y_offset, window_size, slice_width, rgba=True):
    stack = img
    # Crop image centered to point with size that keeps circle of window
    stack_cropped = stack[y - window_size + y_offset:y + window_size, x - window_size + x_offset:x + window_size, ...]
    if rgba:
        stack_padded = gray2rgba(np.zeros([2 * window_size, 2 * window_size], dtype='float64'))
    else:
        stack_padded = np.zeros([2 * window_size, 2 * window_size])

    if rgba:
        stack_padded[y_offset:stack_cropped.shape[0] + y_offset,
        x_offset:stack_cropped.shape[1] + x_offset] = stack_cropped
    else:
        stack_padded[y_offset:stack_cropped.shape[0] + y_offset,
        x_offset:stack_cropped.shape[1] + x_offset] = stack_cropped

    # Apply rotation to original cropped image
    stack_rotated = rotate(stack_padded, angle, preserve_range=True)[:,
                    window_size - slice_width:window_size + slice_width, ...]

    return stack_rotated


dir_images = 'data/processed/4i/MTU_results'
dir_results = 'data/processed/4i/laminar_windows'
df_meta = pd.read_csv('data/processed/4i/laminator/results_mtu/df_meta.csv')

points = os.listdir(dir_images)

def process_point_single_img(point):
    n_clusters = 32
    np.random.seed(seed=5)
    palette = mpl.colors.ListedColormap(np.random.rand(n_clusters, 3), N=n_clusters)
    palette.set_under("black")

    random.seed(5)
    random.shuffle(palette.colors)
    norm = mpl.colors.Normalize(vmin=1, vmax=n_clusters - 1)

    img_path = Path(dir_images, point, 'labelimg.tif')
    label_img = io.imread(img_path)
    for i in range(1, n_clusters):
        dir_results_mtu = dir_results + '_mtu_' + str(i)
        Path(dir_results_mtu).mkdir(exist_ok=True)
        label_mask = np.where(label_img == i, 1, 0)
        img = palette(norm(label_img))
        img = img * label_mask[..., None]
        df_selected = df_meta[df_meta['organoid'] == int(point)]
        export_slices(df_selected, img, dir_results_mtu, rgba=True)


# Start the process pool and do the computation
with Pool(processes=50) as pool:
    pool.map(process_point_single_img, points)

# Detect wedges with wholes
wedge_files = os.listdir(dir_results)

def quantify_background(img):
    img_gray = rgb2gray(rgba2rgb(img))
    img_binary = np.where(img_gray == 0, 1, 0)
    pixel_count_bg = np.sum(img_binary[int(img_binary.shape[0] / 2):, ])
    return pixel_count_bg

backgrounds = []
for file in tqdm(wedge_files):
    img = io.imread(Path(dir_results, file))
    backgrounds.append(quantify_background(img))

df = pd.DataFrame({'wedge': wedge_files, 'bg_count': backgrounds})
df.to_csv('data/processed/4i/laminator/analysis_results/df_bg_counts.csv', index=False)
