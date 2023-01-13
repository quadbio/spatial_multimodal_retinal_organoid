import os
from skimage.filters import gaussian, threshold_otsu
from skimage import img_as_uint
import matplotlib.pyplot as plt
from skimage import io
from skimage import measure
from scipy import ndimage
from pathlib import Path
from skimage.filters.rank import median
from skimage.morphology import disk
from tqdm import tqdm
from image_processing_pipeline.modules import scale_image
import numpy as np
import pandas as pd

def run_masking(dir_output='data/processed/fish/masks_from_spots',
                save=True,
                plot=True,
                save_plot=True,
                show_plot=False):

    dir_images = 'data/raw/fish/dapi_images'
    files = os.listdir(dir_images)
    files = [file for file in files if '._' not in file]

    dir_tables = 'data/raw/fish/spot_data'
    files_tables = [file for file in os.listdir(dir_tables) if not file.startswith('.') and 'results.txt' in file]

    for file in tqdm(files, ascii=True):
        img_init = io.imread(Path(dir_images, file))
        point = file.replace('32955-slide928-2_', '').replace('_DAPI.tiff', '')
        file_table = [file for file in files_tables if point in file][0]
        df = pd.read_csv(Path(dir_tables, file_table), sep='\t', names=['x', 'y', 'z', 'transcript'], index_col=False)

        # generate binary array stacks from spot dataframe
        transcripts = df['transcript'].unique().tolist()
        dict_transcripts = dict((i, j) for i, j in enumerate(transcripts))

        img_spots = np.zeros((img_init.shape[0], img_init.shape[1]))
        for i, transcript in dict_transcripts.items():
            df_tmp = df[df['transcript'] == transcript]
            for x, y in zip(df_tmp['x'].values, df_tmp['y'].values):
                img_spots[y, x] = 1

        struct = ndimage.generate_binary_structure(2, 2)
        img_spots = np.pad(img_spots, [50, 50], mode='constant')
        img_spots = ndimage.binary_dilation(img_spots, structure=struct, iterations=50).astype('uint8')
        img_spots = ndimage.binary_fill_holes(img_spots).astype('uint8')
        img_spots = ndimage.binary_erosion(img_spots, structure=struct, iterations=50).astype('uint8')

        # get mask of empty tiles
        mask_tiles = np.where(img_init > 0, 0, 1)

        mask_tiles = ndimage.binary_erosion(mask_tiles, structure=struct, iterations=10).astype('uint8')
        mask_tiles = ndimage.binary_dilation(mask_tiles, structure=struct, iterations=10).astype('uint8')
        img_init = scale_image(img_init)
        img = median(img_spots, disk(30))

        img = img[50:-50,50:-50].astype('int8') - mask_tiles
        img = np.where(img > 0, 1, 0)
        # Select biggest area
        img = measure.label(img)
        regions = measure.regionprops(img)
        regions.sort(key=lambda x: x.area, reverse=True)
        if len(regions) > 1:
            for rg in regions[1:]:
                img[rg.coords[:, 0], rg.coords[:, 1]] = 0

        # Smoothing selected area
        # Dilation
        img = ndimage.binary_dilation(img, structure=struct, iterations=50).astype('uint8')

        # Apply gaussian
        img = gaussian(img, sigma=50)

        # Otsu thresholding
        thr = threshold_otsu(img)
        img = (img > thr).astype(int)
        # Fill holes
        img = ndimage.binary_fill_holes(img).astype(int)
        img = img - mask_tiles
        img = np.where(img > 0, 1, 0)
        # Save mask
        if save:
            Path(dir_output).mkdir(parents=True, exist_ok=True)
            io.imsave(Path(dir_output, file), img_as_uint(img), check_contrast=False)

        # Plot mask
        if plot:
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 10))

            ax[0].imshow(img_init)
            ax[0].axis('off')
            ax[0].set_title('Init DAPI')

            ax[1].imshow(img)
            ax[1].axis('off')
            ax[1].set_title('Mask')

            ax[2].imshow(img_init * img)
            ax[2].axis('off')
            ax[2].set_title('Masked image')
            if save_plot:
                Path(dir_output, 'plots').mkdir(parents=True, exist_ok=True)
                plt.savefig(Path(dir_output, 'plots',file.replace('.tiff','') + '.png'))
            if show_plot:
                plt.show()
            plt.close('all')


if __name__ == '__main__':
    run_masking()

