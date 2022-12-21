def run_masking(dir_input='/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission',
                dir_output='/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission/masks_3',
                save=True,
                plot=True, save_plot=True, show_plot=False):
    # Import functions
    import os
    from skimage.filters import (gaussian, threshold_otsu)
    from skimage import img_as_uint
    import matplotlib.pyplot as plt
    from skimage import io
    from skimage import measure
    from scipy import ndimage
    from pathlib import Path
    from skimage.filters.rank import median
    from skimage.morphology import disk
    from tqdm import tqdm
    from image_processing.modules import scale_image
    import numpy as np

    files = os.listdir(dir_input)
    files = [file for file in files if '._' not in file]
    files = [file for file in files if 'raw.tiff' in file]

    for file in tqdm(files):
        img_init = io.imread(Path(dir_input, file))
        img_init = scale_image(img_init)
        # get mask of empty tiles
        mask_tiles = np.where(img_init > 0, 0, 1)
        struct = ndimage.generate_binary_structure(2, 2)
        mask_tiles = ndimage.binary_erosion(mask_tiles, structure=struct, iterations=10).astype('uint8')
        mask_tiles = ndimage.binary_dilation(mask_tiles, structure=struct, iterations=10).astype('uint8')

        # Apply gaussian
        img = gaussian(img_init, sigma=100)

        # Otsu thresholding
        thr = threshold_otsu(img)
        img = (img > thr).astype('uint8')

        # Fill holes
        img = ndimage.binary_fill_holes(img).astype('uint8')
        img = img + mask_tiles
        img = ndimage.binary_dilation(img, structure=struct, iterations=50).astype('uint8')
        img = ndimage.binary_fill_holes(img).astype('uint8')
        img = ndimage.binary_erosion(img, structure=struct, iterations=50).astype('uint8')
        img = img - mask_tiles
        # Median filter to remove outliers
        img = median(img, disk(30))
        # Select biggest area
        img = measure.label(img)
        regions = measure.regionprops(img)
        regions.sort(key=lambda x: x.area, reverse=True)
        if len(regions) > 1:
            for rg in regions[1:]:
                img[rg.coords[:, 0], rg.coords[:, 1]] = 0

        # Smoothing selected area
        # Dilation
        img = ndimage.binary_dilation(img, structure=struct, iterations=25).astype('uint8')

        # Apply gaussian
        img = gaussian(img, sigma=50)

        # Otsu thresholding
        thr = threshold_otsu(img)
        img = (img > thr).astype(int)
        # Fill holes
        img = ndimage.binary_fill_holes(img).astype(int)
        img = img - mask_tiles

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