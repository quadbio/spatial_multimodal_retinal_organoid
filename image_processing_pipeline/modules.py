import yaml
from skimage import img_as_uint
from scipy import ndimage
from skimage.filters.rank import median
from skimage.morphology import disk
from time import gmtime, strftime
import re
import cv2
import SimpleITK as sitk
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.filters import gaussian, threshold_otsu, threshold_multiotsu
from skimage import measure
from scipy import ndimage as ndi
from fnnls import fnnls
from tqdm import tqdm
from cellpose import models
import matplotlib.patches as mpatches
from scipy.stats import zscore
import seaborn as sns
import random
import os.path
import copy
import numpy as np
import phenograph
import matplotlib as mpl
import multiprocessing
from joblib import Parallel, delayed
from sklearn.metrics.pairwise import euclidean_distances
from skimage.measure import regionprops_table
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io

# load global variables
with open("params.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
globals().update(cfg)


# generic functions
def sorted_nicely(l):
    # sort list of strings alphanumerically as intuited by humans
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def get_metadata(dir_images):
    images = os.listdir(dir_images)
    images = [image for image in images if '.tif' in image]
    images.sort()
    regex = r'cycle(?P<cycle_id>\d+)_well(?P<well_id>\d{2})_channel(?P<channel>\d{1})_(?P<stain>[a-zA-Z0-9]+)'
    df = pd.DataFrame({'file': images})
    df = df.join(df['file'].astype(str).str.extractall(regex).groupby(level=0).last())
    df['well_id'] = df['well_id'].apply(lambda x: int(x))
    df['cycle_id'] = df['cycle_id'].apply(lambda x: int(x))
    df['channel'] = df['channel'].apply(lambda x: int(x))
    return df


def scale_image(image, percentile=1):
    image = np.interp(image, (np.percentile(image, percentile), np.percentile(image, 100 - percentile)), (0, +65535))
    return image


#################-----------------------------------------------------------------########################
# get a simple initial mask that helps elastix alignment
def simple_mask(point, dir_input, dir_output, cycle='cycle1', save=True, plot=False, show_plot=False, save_plot=False):
    """ Creating an initial simple mask

        :param point: (int) id of organoid
        :param dir_input: (str) path to input directory
        :param dir_output: (str) path to output directory
        :param cycle: (str) cycle to create mask from
        :param save: (bool) save created mask
        :param plot: (bool) create plot
        :param save_plot: (bool) save plot
        :param show_plot: (bool) return plot

        :return: depending on input parameters the function saves the created mask and/or plots and/or shows the overview plot

    """


    if not os.path.isfile(str(Path(dir_output, 'simple_initial_mask' + str(point) + '.tif'))):

        # get filename
        img_df = get_metadata(str(Path(data_path, 'raw_data')))
        filename = img_df.loc[(img_df['well_id'] == point) & (img_df['cycle_id'] == cycle) & (img_df['stain'] == 'hoechst')][
            'file'].values[0]

        # Load image

        img_init = io.imread(str(Path(dir_input, 'raw_data', filename)))

        img_init = np.apply_over_axes(scale_image, img_init, 2)

        # Average all channels
        # img_init = np.max(img_init, axis=2)

        # Apply gaussian
        img = gaussian(img_init, sigma=100)

        # Otsu thresholding
        thr = threshold_otsu(img)
        img = (img > thr).astype(int)

        # Fill holes
        img = ndimage.binary_fill_holes(img).astype('uint8')

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
        struct = ndimage.generate_binary_structure(2, 2)
        img = ndimage.morphology.binary_dilation(img, structure=struct, iterations=50).astype(int)

        # Apply gaussian
        img = gaussian(img, sigma=50)

        # Otsu thresholding
        thr = threshold_otsu(img)
        img = (img > thr).astype(int)

        # Fill holes
        img = ndimage.binary_fill_holes(img).astype(int)

        # Save mask
        if save:
            Path(dir_output).mkdir(parents=True, exist_ok=True)
            io.imsave(str(Path(dir_output, 'simple_initial_mask' + str(point) + '.tif')), img.astype('bool'),
                      check_contrast=False)

        # Plot mask
        if plot:
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 10))

            ax[0].imshow(img_init)
            ax[0].axis('off')
            ax[0].set_title('Max over all channels')

            ax[1].imshow(img)
            ax[1].axis('off')
            ax[1].set_title('Mask')

            ax[2].imshow(img_init * img)
            ax[2].axis('off')
            ax[2].set_title('Masked image')
            if save_plot:
                Path(dir_output, 'plots').mkdir(parents=True, exist_ok=True)
                plt.savefig(str(Path(dir_output, 'plots', str(point) + '.png')))
            if show_plot:
                plt.show()
    else:
        print('sample ' + str(point) + ' already has a simple mask, moving on')


#################-----------------------------------------------------------------########################
# run elastix alignment

def run_elastix(point, ref_cycles, dir_mask, dir_output, cycles=cycles,
                dir_input='/links/groups/treutlein/DATA/imaging/charmel/publication_example_data/raw_data/'):
    """
     Running Elastix and Transformix for dataset of one organoid and saves them to output directory
    :param point: (int) id of organoid
    :param ref_cycles: (int) cycles which are used as reference cycles
    :param dir_mask: (str) path to directory of initial simple masks
    :param dir_output: (str) path to output directory
    :param dir_input: (str) path to input directory
    :return: None
    """

    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    print('Started alignment for point:', point)
    print('Set initial reference to:', ref_cycles[0])

    # image info df
    img_df = get_metadata(str(Path(data_path, 'raw_data')))

    # Set reference cycle
    ref_imgs = img_df.loc[(img_df['well_id'] == point) & (img_df['cycle_id'] == ref_cycles[0])]['file'].to_list()

    # Load and safe initial reference cycle to output dir
    if not os.path.isfile(Path(dir_output, ref_imgs[3])):
        for n in ref_imgs:
            fixed_img = io.imread(str(Path(dir_input, n)))
            Path(dir_output).mkdir(parents=True, exist_ok=True)
            io.imsave(str(Path(dir_output, n)), fixed_img, check_contrast=False)

    # load fixed hoechst image
    r = re.compile(".*hoechst")
    filename = list(filter(r.match, ref_imgs))
    fixed_img = io.imread(str(Path(dir_input, filename[0])))

    # Load fixed mask
    fixed_mask = io.imread(str(Path(dir_mask, 'simple_initial_mask' + str(point) + '.tif'))).astype('int')
    initial_mask = fixed_mask.copy()
    print('Loaded fixed image and corresponding mask.')

    # Initialize elastix and transformix
    elastix_filter = sitk.ElastixImageFilter()
    elastix_filter.SetParameterMap(sitk.ReadParameterFile('param_maps/translation_pw.txt'))
    elastix_filter.AddParameterMap(sitk.ReadParameterFile('param_maps/affine_pw.txt'))
    elastix_filter.LogToConsoleOff()
    transformix_filter = sitk.TransformixImageFilter()
    print('Initialized Elastix.')

    # remove reference cycle from cycles (global params) to align
    align_cycles = copy.copy(cycles)
    align_cycles.remove(ref_cycles[0])

    # Iteratively apply alignment across all cycles
    for cycle in align_cycles:
        filenames = img_df[(img_df['cycle_id'] == cycle) & (img_df['well_id'] == point)]['file'].to_list()

        items = []
        for i in filenames:
            items.append(Path(data_path, 'aligned', i))

        if not (all(os.path.isfile(i) for i in items)):
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
            print('Aligning sample:', point, '... cycle ', cycle)

            # Load moving image
            print('Loading moving images and masks.')
            r = re.compile(".*hoechst")
            filename = list(filter(r.match, filenames))
            moving_img = io.imread(str(Path(dir_input, filename[0])))
            moving_mask = np.where(moving_img > 1, 1, 0)
            moving_mask = ndimage.binary_fill_holes(moving_mask).astype(int)

            # Elastix
            print('Running Elastix.')
            elastix_filter.SetFixedImage(sitk.GetImageFromArray(fixed_img))
            elastix_filter.SetFixedMask(sitk.Cast(sitk.GetImageFromArray(fixed_mask), sitk.sitkUInt8))
            elastix_filter.SetMovingImage(sitk.GetImageFromArray(moving_img))
            elastix_filter.SetMovingMask(sitk.Cast(sitk.GetImageFromArray(moving_mask), sitk.sitkUInt8))
            # elastix_filter.SetOutputDirectory(dir_output)
            elastix_filter.Execute()

            # Transformix
            print('Running Transformix.')
            channels = []
            for channel in range(4):
                # Update and execute Transformix
                moving_channel = io.imread(str(Path(dir_input, filenames[channel])))

                transformix_filter.SetTransformParameterMap(elastix_filter.GetTransformParameterMap())
                transformix_filter.SetMovingImage(sitk.GetImageFromArray(moving_channel))
                transformix_filter.SetLogToConsole(False)
                channel_aligned = transformix_filter.Execute()
                # Convert to numpy array
                channel_aligned = sitk.GetArrayFromImage(channel_aligned)
                # Cap values just in case
                channel_aligned[channel_aligned < 0] = 0
                channel_aligned[channel_aligned > 65535] = 65535
                # Convert to uint16
                channel_aligned = channel_aligned.astype(np.uint16)

                # Save aligned image
                io.imsave(str(Path(dir_output, filenames[channel])), channel_aligned, check_contrast=False)
                print('Saved aligned image.')

            # Set new reference cycle if necessary

        if cycle in ref_cycles:
            ref_cycle = cycle
            print('Setting reference to:', ref_cycle)

            r = re.compile(".*hoechst")
            filename = list(filter(r.match, filenames))
            fixed_img = io.imread(Path(dir_output, filename[0]))
            new_mask = fixed_img.copy()
            new_mask = np.where(new_mask > 1, 1, 0)
            new_mask = ndimage.binary_fill_holes(new_mask).astype(int)
            fixed_mask = initial_mask * new_mask
            print('Loaded new fixed image and mask.')
        # if (cycle in ref_cycles) & (os.path.isfile(Path(dir_output,filenames[2]))):
        #    r = re.compile(".*hoechst")
        #    filename = list(filter(r.match, filenames))
        #    print('sample ', point, ' cycle ',cycle, ' already done.loading reference image', filename[0])
        #    fixed_img = io.imread(Path(dir_output,filename[0]))


# denoise images, crop to size of simple mask to speed up the procedure
def crop_image(img='img', mask='mask'):
    mask = mask > 0
    return img[np.ix_(mask.any(1), mask.any(0))]


def denoise(image):
    sigma_est = np.mean(estimate_sigma(image, multichannel=False))

    patch_kw = dict(patch_size=5,  # 5x5 patches
                    patch_distance=6,  # 13x13 search area
                    multichannel=False)

    denoise2_fast = denoise_nl_means(image, h=0.8 * sigma_est, sigma=sigma_est,
                                     fast_mode=True, **patch_kw)

    return (denoise2_fast)


def denoise_well(well=0):
    mask = io.imread(data_path + 'masks/simple_initial_mask' + str(well) + '.tif').astype('int')

    if not os.path.isfile(data_path + 'masks/cropped_simple_mask' + str(well) + '.tif'):
        cropped_mask = crop_image(img=mask, mask=mask).astype('uint8')
        io.imsave(data_path + 'masks/cropped_simple_mask' + str(well) + '.tif', cropped_mask.astype('bool'),
                  check_contrast=False)

    else:
        cropped_mask = io.imread(data_path + 'masks/cropped_simple_mask' + str(well) + '.tif').astype('uint8')

    img_df = get_metadata(Path(data_path, 'raw_data'))
    filenames = img_df.loc[(img_df['well_id'] == well)]['file'].to_list()
    filenames = sorted_nicely(filenames)
    for i in filenames:
        if not os.path.isfile(Path(data_path, 'denoised', i)):
            print('denoising ' + i)
            img = io.imread(Path(data_path, 'aligned/', i))

            cropped = crop_image(img=img, mask=mask)
            tmp = denoise(cropped)

            del img
            os.makedirs(data_path + 'denoised/', exist_ok=True)
            io.imsave(str(Path(data_path, 'denoised', i)), tmp.astype('uint16'), check_contrast=False)


# create a tight mask in order to restrict the analysis to pixels inside the tissue that were imaged in every cycle
def refine_mask(well=0):
    if not os.path.isfile(data_path + '/masks/cropped_refined_mask' + str(well) + '.tif'):
        print('refining mask ' + str(well))
        old_mask = io.imread(data_path + '/masks/cropped_simple_mask' + str(well) + '.tif').astype('bool')
        imgs = []

        img_df = get_metadata(str(Path(data_path, 'denoised')))
        filenames = img_df.loc[(img_df['well_id'] == well) & (img_df['cycle_id'].isin([1, 2, 3, 4, 5, 7, 8, 9, 10, 11]))][
            'file'].to_list()
        filenames = sorted_nicely(filenames)
        for i in filenames:
            imgs.append(io.imread(Path(data_path, 'denoised/', i)))

        for i in np.arange(len(imgs)):
            tmp = scale_image(imgs[i]) / len(imgs)
            if i == 0:
                img_avg = tmp
            else:
                img_avg = tmp + img_avg

        threshold = threshold_multiotsu(old_mask * img_avg, classes=5)[0]

        img_thr_sharp = (img_avg >= threshold)
        img_mask = img_thr_sharp * old_mask

        img_mask = ndi.binary_fill_holes(img_mask).astype(int)

        struct = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (10, 10)).astype('bool')
        img_mask = ndi.morphology.binary_dilation(img_mask, structure=struct, iterations=15).astype('int')

        labels_mask = measure.label(img_mask)

        regions = measure.regionprops(labels_mask)
        regions.sort(key=lambda x: x.area, reverse=True)
        if len(regions) > 1:
            for rg in regions[1:]:
                labels_mask[rg.coords[:, 0], rg.coords[:, 1]] = 0
        labels_mask[labels_mask != 0] = 1
        mask2 = labels_mask

        threshold = threshold_multiotsu(mask2 * img_avg, classes=5)[0] - 100

        img_thr_sharp = (img_avg >= threshold)
        img_mask = img_thr_sharp * mask2

        img_mask = ndi.binary_fill_holes(img_mask).astype(int)

        struct = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)).astype('bool')
        umg_mask = ndi.morphology.binary_dilation(img_mask, structure=struct, iterations=2).astype('int')

        umg_mask = ndi.binary_fill_holes(umg_mask).astype(int)

        labels_mask = measure.label(umg_mask)
        regions = measure.regionprops(labels_mask)
        regions.sort(key=lambda x: x.area, reverse=True)
        if len(regions) > 1:
            for rg in regions[1:]:
                labels_mask[rg.coords[:, 0], rg.coords[:, 1]] = 0
        labels_mask[labels_mask != 0] = 1
        mask3 = labels_mask * old_mask

        # remove pixels from mask that were cut off during imaging in different cycles
        filenames = img_df.loc[(img_df['well_id'] == well) & (img_df['stain'] == 'hoechst')]['file'].to_list()
        for i in filenames:
            tmp = io.imread(Path(data_path, 'denoised/', i))
            tmp[tmp > 0] = 1
            mask3 = mask3 * tmp

        Path(data_path + '/masks/').mkdir(parents=True, exist_ok=True)
        io.imsave(data_path + '/masks/cropped_refined_mask' + str(well) + '.tif', mask3.astype('uint8'),
                  check_contrast=False)


# remove well and cycle specific backgrounds that are proportionally composed by combining the previous and
# consecutive elution cycle for every image per channel. Also move the reference dapi image to the same folder
# because this folder will have all the fully processed images used in the analysis
def bg_subtraction(well,
                   dir_input=Path(data_path, 'denoised'),
                   dir_masks=Path(data_path, 'masks'),
                   dir_output=Path(data_path, 'bg_subtracted'),
                   cycles_bg=[0, 6, 12],
                   cycle_hoechst=1):


    print('Subtracting backgrounds in sample :', well)
    # Get metadata, set paths and create directories
    df = get_metadata(dir_input)
    dir_output.mkdir(parents=True, exist_ok=True)

    # remove bg cycles from cycles
    cycles_fg = copy.copy(cycles)
    for i in cycles_bg:
        cycles_fg.remove(i)

    # Get file names of reference dapi and dfs for foreground and background images
    file_dapi = \
    df[(df['well_id'] == well) & (df['stain'] == 'hoechst') & (df['cycle_id'] == cycle_hoechst)]['file'].values[0]
    df_bg = df[(df['well_id'] == well) & (df['cycle_id'].isin(cycles_bg)) & (df['stain'] != 'hoechst')]
    df = df[(df['well_id'] == well) & (df['cycle_id'].isin(cycles_fg)) & (df['stain'] != 'hoechst')]
    df = df[df['channel'].isin([0, 1, 3])]
    files_images = df['file'].values

    # Load and save dapi image
    img = io.imread(Path(dir_input, file_dapi))
    io.imsave(str(Path(dir_output, file_dapi)), img, check_contrast=False)

    # Load mask
    mask = io.imread(Path(dir_masks, 'cropped_refined_mask' + str(well) + '.tif'))

    # check whether well is already fully processed
    items = []
    for i in files_images:
        items.append(Path(data_path, 'bg_subtracted', i))

    if not (all(os.path.isfile(i) for i in items)):
        for file_image in tqdm(files_images, ascii=True):
            if not os.path.isfile(Path(dir_output, file_image)):
                # read image
                img = io.imread(Path(dir_input, file_image))
                cycle = df[df['file'] == file_image]['cycle_id'].values[0]
                channel = df[df['file'] == file_image]['channel'].values[0]
                # figure out previous and subsequent elustion cycles
                cycles_bg_select = np.abs(cycle - np.array(cycles_bg)).argsort()
                cycle_closer = int(np.array(cycles_bg)[cycles_bg_select[0]])
                cycle_further = int(np.array(cycles_bg)[cycles_bg_select[1]])
                cycles_sorted = sorted([cycle_closer, cycle_further])
                leading_cycle_bg = cycles_sorted[0]
                # read leading elution image
                file_leading_cycle_bg = \
                df_bg[(df_bg['cycle_id'] == leading_cycle_bg) & (df_bg['channel'] == channel)]['file'].values[0]
                img_bg_leading = io.imread(Path(dir_input, file_leading_cycle_bg))

                # read lagging elution cycle
                lagging_cycle_bg = cycles_sorted[1]
                file_lagging_cycle_bg = \
                df_bg[(df_bg['cycle_id'] == lagging_cycle_bg) & (df_bg['channel'] == channel)]['file'].values[0]
                img_bg_lagging = io.imread(Path(dir_input, file_lagging_cycle_bg))
                # compute background while weighting leading and trailing background images appropriately
                cycle_range = lagging_cycle_bg - leading_cycle_bg
                factor_leading = (cycle_range - (cycle - leading_cycle_bg)) / cycle_range
                factor_lagging = (cycle_range - (lagging_cycle_bg - cycle)) / cycle_range
                img_bg = (img_bg_lagging * factor_lagging) + (img_bg_leading * factor_leading)
                # compute factor for background subtraction and subtract
                bg_factor = fnnls(img_bg[mask == 1][..., np.newaxis], img[mask == 1])
                img_bg = img_bg * bg_factor[0]
                img = img - img_bg
                img[img < 0] = 0
                img[img > 65535] = 65535
                img = img * mask
                img_subtracted = img.astype(np.uint16)

                io.imsave(str(Path(dir_output, file_image)), img_subtracted, check_contrast=False)
            else:
                print(file_image + 'already processed. skipping')
    else:
        print('well ' + str(well) + ' already fully background subtracted. skipping')


# segment nuclei. It is advised to perform this on the aligned and masked reference dapi image
def segment_nuclei(well,
                   ref_cycle,
                   dir_input=Path(data_path, 'denoised'),
                   dir_output=Path(data_path, 'segmented_nuclei')):


    if not os.path.isfile(data_path + 'segmented_nuclei/' + str(well) + '.tif'):
        df = get_metadata(dir_input)
        filename = \
        df[(df['cycle_id'] == ref_cycle) & (df['well_id'] == well) & (df['stain'] == 'hoechst')]['file'].values[0]
        print('segmenting nuclei; ' + filename)

        mask = io.imread(data_path + '/masks/cropped_refined_mask' + str(well) + '.tif')
        img = io.imread(Path(dir_input, filename)).astype('int') * mask

        thr = threshold_otsu(img)
        img[img < thr] = 0

        model = models.Cellpose(gpu=False, model_type="cyto")
        diam = 17
        flow_thr = 0.8
        cellprob_thr = 0
        channels = [0, 0]

        masks, flows, styles, diams = model.eval(img, diameter=diam, channels=channels, do_3D=False,
                                                 flow_threshold=flow_thr, cellprob_threshold=cellprob_thr)

        os.makedirs(dir_output, exist_ok=True)
        io.imsave(str(Path(dir_output, str(well) + '.tif')), masks)
    else:
        print('well ' + str(well) + ' nuclei already segmented')


# generate sample specific pixel matrices based on raw pixel intensities of the fully processed images.
# Optionally remove artifacts that reside in the collagen rich areas.
def generate_pixel_matrices(well=0,
                            dir_input=Path(data_path, 'bg_subtracted'),
                            dir_output=Path(data_path, 'pixel_matrices'),
                            cycles_bg=[0, 6, 12],
                            to_collagen_mask=to_collagen_mask,
                            mask_collagen=True,
                            overwrite_denoised_images_with_collagen_masked_images=True):


    mask = io.imread(Path(data_path, 'masks/cropped_refined_mask' + str(well) + '.tif')).astype('bool')

    stain_cycles = copy.copy(cycles)
    for i in cycles_bg:
        stain_cycles.remove(i)

    img_df = get_metadata(dir_input)
    filename_hoechst = img_df.loc[(img_df['well_id'] == well) & (img_df['stain'] == 'hoechst')]['file'].values[0]
    filenames = sorted_nicely(img_df[(img_df['well_id'] == well) & (img_df['stain'] != 'hoechst')]['file'].to_list())
    stains = []
    for i in filenames:
        stains.append(img_df[img_df['file'] == i]['stain'].values[0])
    stains.insert(0, 'hoechst')

    with open('metadata.yml', 'w') as outfile:
        yaml.dump(stains, outfile, default_flow_style=False)

    to_collagen_mask = img_df[(img_df['well_id'] == well) & (img_df['stain'].isin(to_collagen_mask))]['file'].to_list()

    img = io.imread(Path(dir_input, filename_hoechst)).astype('int16')
    full_matrix = img[mask][:, None]

    if not os.path.isfile(Path(dir_output, str(well) + 'raw_pixel_matrix.npz')):
        print('generating pixel matrix ' + str(well))
        images = []
        for filename in filenames:

            if (mask_collagen == True) & (filename in to_collagen_mask):
                colname = img_df.loc[(img_df['well_id'] == well) & (img_df['stain'] == 'COL2A1')]['file'].values[0]

                img = io.imread(str(dir_input / filename)).astype('uint16')
                collagen_mask = get_collagen_mask(well, colname)

                img = img * collagen_mask

                mat = img[mask][:, None]
                full_matrix = np.concatenate((full_matrix, mat), axis=1)

                if overwrite_denoised_images_with_collagen_masked_images == True:
                    print('saving collagen masked image ' + filename)
                    io.imsave(str(dir_input / filename), img)

            else:
                img = io.imread(Path(dir_input, filename)).astype('uint16')
                mat = img[mask][:, None]
                full_matrix = np.concatenate((full_matrix, mat), axis=1)

        pixel_matrix = full_matrix
        dir_output.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(str(Path(dir_output, str(well) + 'raw_pixel_matrix.npz')), pixel_matrix)
    else:
        print(str(well) + ' pixel matrix already exists')


def get_collagen_mask(well, colname):
    if not os.path.isfile(data_path + 'refined_masks/collagen_mask' + str(well) + '.tif'):

        print('creating collagen mask ' + colname + ' this might take a while')

        collagen = io.imread(Path(data_path, 'bg_subtracted', colname))
        tmp = copy.copy(collagen)
        thresh = threshold_multiotsu(collagen[::10, ::10], 4)[0]
        tmp[tmp < thresh] = 0
        tmp[tmp >= thresh] = 1
        tmp = ndi.binary_fill_holes(tmp).astype(int)
        tmp = tmp.astype('uint8')

        nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(tmp, connectivity=8)
        sizes = stats[1:, -1];
        nb_components = nb_components - 1
        min_size = 25
        img2 = np.zeros((tmp.shape))

        for i in range(0, nb_components):
            if sizes[i] >= min_size:
                img2[output == i + 1] = 1

        kernel = np.ones((5, 5), np.uint8)
        colmask = cv2.erode(img2, kernel, iterations=5).astype('int8')

        Path(data_path + '/refined_masks/').mkdir(parents=True, exist_ok=True)
        io.imsave(str(Path(data_path, 'refined_masks/collagen_mask' + str(well) + '.tif')), colmask, check_contrast=False)
    else:
        colmask = io.imread(data_path + 'refined_masks/collagen_mask' + str(well) + '.tif')

    return colmask


# normalize matrices to maintain inter condition differences and reduce intra condition differences.
# Also generate downsampled, normalized and scaled fusion matrix to generate the reference for MTU clustering and
# scale individual matrices in the same way to allow MTU assignment to sample pixels
def run_matrix_normalization(wells,
                             dir_input=Path(data_path, 'pixel_matrices'),
                             dir_output=Path(data_path, 'pixel_matrices'),
                             fusion_matrix=True):  # not currently used

    # check whether all samples already have a final matrix
    items = []
    for i in wells:
        items.append(Path(dir_output, str(i) + 'scaled_z_normalised_pixel_matrix.npz'))

    if not (all(os.path.isfile(i) for i in items)):
        import random
        from scipy import stats
        matrices = []
        for well in wells:
            matrices.append(np.load(Path(dir_input, str(well) + 'raw_pixel_matrix.npz'))['arr_0'])

        long_matrix = np.concatenate(matrices, axis=0)

        mean = np.mean(long_matrix, axis=0)
        std = np.std(long_matrix, axis=0)

        subs = []
        normalized_mats = []
        for i in np.arange(len(matrices)):
            normalised_values = []
            for n in np.arange(matrices[i].shape[1]):
                zscores = stats.zscore(matrices[i][:, n])
                normalised_values.append(zscores * std[n] + mean[n])

            normalised_mat = np.stack(normalised_values, axis=1)
            normalized_mats.append(normalised_mat)
            # np.savez_compressed(Path(data_path, 'pixel_matrices', str(wells[i]) +'z_normalised_pixel_matrix.npz'), normalised_mat.astype('float32'))

            pick = round(len(matrices[i]) / 250)
            random.seed(1)
            rand = np.array([random.randint(0, len(matrices[i])) for _ in range(pick)])
            pixel_matrix_sub = matrices[i][rand, :]
            subs.append(pixel_matrix_sub)
            # np.savez_compressed(Path(dir_output, str(well) + '_sub250_z_normalised_pixel_matrix.npz'), pixel_matrix_sub.astype('float32'))

        # generate normalized subsampled fusion matrix
        long_matrix = np.concatenate(subs, axis=0)

        bottom_percentiles, top_percentiles = get_upper_and_lower_percentiles(long_matrix)

        # scale and save fusion matrix
        scaled_long_matrix = scale_matrix(long_matrix, bottom_percentiles, top_percentiles)
        np.savez_compressed(str(Path(dir_output, 'scaled_fusion_sub250_z_normalised_pixel_matrix.npz')),
                            scaled_long_matrix.astype('float32'))

        # scale and save sample matrices
        for i in np.arange(len(matrices)):
            mat = scale_matrix(normalized_mats[i], bottom_percentiles, top_percentiles)
            np.savez_compressed(str(Path(dir_output, str(wells[i]) + 'scaled_z_normalised_pixel_matrix.npz')),
                                mat.astype('float32'))


def get_upper_and_lower_percentiles(matrix, percentiles=[1, 99]):
    bottom_percentiles = []
    top_percentiles = []
    for i in np.arange(matrix.shape[1]):
        bottom_percentiles.append(np.percentile(matrix[:, i], percentiles[0]))
        top_percentiles.append(np.percentile(matrix[:, i], percentiles[1]))
    return [bottom_percentiles, top_percentiles]


def scale_matrix(matrix, bottom_percentiles, top_percentiles):
    tmp = []
    for i in np.arange(matrix.shape[1]):
        tmp.append(np.interp(matrix[:, i], (bottom_percentiles[i], top_percentiles[i]), (0, +1)))
        scaled_mat = np.stack(tmp, axis=1)
    return scaled_mat


#################-----------------------------------------------------------------########################
# create MTUs and associated plots.
def MTU_assignment(well=43,
                   k=20,
                   dir_matrix=Path(data_path, 'pixel_matrices'),
                   dir_output=Path(data_path, 'MTU_results'),
                   dir_fSOM=Path(data_path, 'fSOM_output')):

    print('assigning MTUs in sample ' + str(well))

    class MidpointNormalize(mpl.colors.Normalize):
        def __init__(self, vmin, vmax, midpoint=0, clip=False):
            self.midpoint = midpoint
            mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

        def __call__(self, value, clip=None):
            normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
            normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
            normalized_mid = 0.5
            x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
            return np.ma.masked_array(np.interp(value, x, y))

    # load data
    os.makedirs(dir_output, exist_ok=True)
    som_codes = np.load(Path(dir_fSOM, 'som_codes.npz'))['arr_0']
    mask = io.imread(Path(data_path, 'masks', 'cropped_refined_mask' + str(well) + '.tif')).astype('bool')

    pixel_matrix = np.load(Path(dir_matrix, str(well) + 'scaled_z_normalised_pixel_matrix.npz'))['arr_0']

    # exclude unwanted stains from pixel matrix
    with open(Path("metadata.yml"), 'r') as ymlfile:
        colnames = yaml.load(ymlfile, Loader=yaml.FullLoader)
    keep = [colnames.index(x) for x in colnames if x != 'excluded']
    marker_names = [x for x in colnames if x != 'excluded']
    pixel_matrix = pixel_matrix[:, keep]

    # load or compute phenograph reference clusters to which image pixels will be assigned
    if not os.path.isfile(Path(dir_output, 'phgr_node_labels.npz')):
        phgr_node_labels = phenograph.cluster(som_codes, k=k)[0]
        np.savez_compressed(Path(dir_output, 'phgr_node_labels.npz'), phgr_node_labels)
    else:
        phgr_node_labels = np.load(Path(dir_output, 'phgr_node_labels.npz'))['arr_0']
    print('phenograph done')

    # make sure phenograph node labels are not negative and start at 1
    cluster_labels_unique = np.unique(phgr_node_labels)
    phgr_node_labels = phgr_node_labels + cluster_labels_unique.min() * -1 + 1
    cluster_labels_unique = np.unique(phgr_node_labels)
    n_clusters = len(cluster_labels_unique)
    # save n_MTUs to yaml file for later access
    with open('n_MTUs.yml', 'w') as outfile:
        yaml.dump(n_clusters, outfile, default_flow_style=False)

    # map pixels to phenograph clusters by euclidean distance using parallel computing (creating MTUs)
    split = np.array_split(pixel_matrix, 500)
    num_cores = multiprocessing.cpu_count()

    def my_function(i):
        tmp = euclidean_distances(i, som_codes)
        tmp2 = np.argmin(tmp, axis=1)
        return tmp2

    indices = Parallel(n_jobs=50)(delayed(my_function)(i) for i in split)

    min_indices = []
    for l in indices:
        min_indices.extend(l)

    pixel_cluster_labels = np.asarray([phgr_node_labels[min_indices[i]] for i in range(len(min_indices))])
    print('MTUs were assigned. ' + str(well))

    # recompose sample image by putting MTU labeled pixels back into the mask. Save labelimage to use for computation
    cluster_img = np.zeros((mask.shape))
    cluster_img[mask] = pixel_cluster_labels

    os.makedirs(Path(dir_output, str(well)), exist_ok=True)
    io.imsave(str(Path(dir_output, str(well), 'labelimg.tif')), cluster_img.astype('uint8'), check_contrast=False)

    # compute sample specific heatmap showing how stains relate to MTUs

    np.random.seed(seed=5)
    palette = mpl.colors.ListedColormap(np.random.rand(n_clusters, 3), N=n_clusters)
    palette.set_under("black")
    patches = [mpatches.Patch(color=palette.colors[i], label=i + 1) for i in range(n_clusters - 1)]

    random.seed(5)
    random.shuffle(palette.colors)

    cluster_avg_list = []

    for cluster in cluster_labels_unique:
        idx = np.where(pixel_cluster_labels == cluster)[
            0]  # tuple where the second element is empty, so only take the first
        mean_vals = np.mean(pixel_matrix[idx], axis=0)
        cluster_avg_list.append(mean_vals)

    cluster_avgs = np.stack(cluster_avg_list, axis=1)
    cluster_avgs = np.nan_to_num(cluster_avgs)
    cluster_avgs_zscored = zscore(cluster_avgs, axis=1)
    df_zscored = pd.DataFrame(cluster_avgs_zscored, columns=cluster_labels_unique, index=marker_names)

    cluster_min = cluster_avgs_zscored.min()
    cluster_max = cluster_avgs_zscored.max()
    norm = MidpointNormalize(vmin=cluster_min, vmax=cluster_max, midpoint=0)

    mpl.use("Agg", force=True)
    print('creating_heatmap')
    cm = sns.clustermap(df_zscored,
                        figsize=(45, 60),
                        yticklabels=1,
                        cbar_pos=(1.22, 0, 0.03, 0),
                        cbar_kws=dict(label="Marker intensity (z-score)", use_gridspec=False),
                        cmap=plt.cm.get_cmap('RdBu').reversed(),
                        norm=norm
                        )

    cm.ax_heatmap.set_xlabel("Meta Cluster number", fontsize=45, labelpad=30)
    cm.ax_heatmap.tick_params(labelsize=45, axis="y", labelrotation=0)
    cm.ax_heatmap.tick_params(labelsize=45, labelcolor="white", axis="x", pad=15)

    ticks = cm.ax_heatmap.get_xticklabels()
    for i in range(n_clusters):
        # Get cluster ID at current x tick
        cluster_id = int(ticks[i].get_text())
        colour = palette.colors[cluster_id - 1]  # i-1 because cluster labels start at 1 and index starts at 0
        # Define and apply colour settings (fc = face colour, ec = edge colour)
        bbox = dict(boxstyle="square", fc=colour, ec=colour, alpha=0.9)
        plt.setp(ticks[i], bbox=bbox)

    hm = cm.ax_heatmap.get_position()  # get heatmap position
    # Decrease height of column dendrogram
    col_dend = cm.ax_col_dendrogram.get_position()
    cm.ax_col_dendrogram.set_position([col_dend.x0, col_dend.y0, col_dend.width, col_dend.height * 0.15])
    # Decrease width of row dendrogram
    row_dend = cm.ax_row_dendrogram.get_position()
    cm.ax_row_dendrogram.set_position(
        [row_dend.x0 + row_dend.width * 0.5, row_dend.y0, row_dend.width * 0.5, row_dend.height])
    # Get color bar
    cbar = cm.ax_heatmap.collections[0].colorbar
    # Make sure color bar is flush with heatmap
    cbar_pos = cbar.ax.get_position()
    cbar.ax.set_position([cbar_pos.x0, hm.y0, cbar_pos.width, hm.height])
    # Set colorbar label and tick size
    cbar.set_label("Marker intensity (z-score)", fontsize=60)
    cbar.ax.tick_params(labelsize=45)

    cm.savefig(str(Path(dir_output, str(well), "MTU_heatmap.pdf")), bbox_inches="tight")
    plt.close()

    print('MTU heatmap created')

    # plot and save pixelperfect MTU tissue image in pdf
    # determine size necessary to plot every image pixle into a single pixle in the pdf
    dpi = 80
    height, width = cluster_img.shape
    figsize = width / float(dpi), height / float(dpi)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # Create colorpalette matching Heatmap annotation
    norm = mpl.colors.Normalize(vmin=1, vmax=n_clusters - 1)
    organoid_img = palette(norm(cluster_img))
    ax.imshow(organoid_img)

    # Add legend optionally
    # legend = plt.legend(handles=patches, loc="upper left", borderaxespad=0, fontsize=50, facecolor="0.25", edgecolor="white", fancybox=False)
    # for text in legend.get_texts():
    #    text.set_color("white")

    ax.set(xlim=[-0.5, width - 0.5], ylim=[height - 0.5, -0.5], aspect=1)
    fig.savefig(str(Path(dir_output, str(well), "MTU_image.pdf")), dpi=dpi, transparent=True)
    plt.close(fig)

    print('MTU image pdf done')

    # create pixelperfect single MTU plots.
    os.makedirs(Path(dir_output, str(well) + '/single_MTU_imgs'), exist_ok=True)

    for cluster in cluster_labels_unique:
        print(cluster)
        cluster_mask = np.ma.masked_where(cluster_img != cluster, cluster_img)
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(cluster_mask, cmap=palette,
                  norm=mpl.colors.Normalize(vmin=1, vmax=n_clusters - 1))
        ax.axis("off")
        ax.set(xlim=[-0.5, width - 0.5], ylim=[height - 0.5, -0.5], aspect=1)

        fig.savefig(str(Path(dir_output, str(well) + '/single_MTU_imgs', str(cluster) + '.png')), dpi=dpi, transparent=True)
        plt.close(fig)

    print(str(well) + ' MTUs done')


#################-----------------------------------------------------------------########################
# generate nuclear feature tables
def percentiles(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], q=(5, 50, 95))


def std(regionmask, intensity_image):
    return np.std(intensity_image[regionmask])


def pixelcount(regionmask, intensity_image):
    return np.sum(intensity_image[regionmask])


def skewness(regionmask, intensity_image):
    from scipy.stats import skew, kurtosis
    return skew(intensity_image[regionmask])


def kurtosisness(regionmask, intensity_image):
    from scipy.stats import skew, kurtosis
    return kurtosis(intensity_image[regionmask])


def clusterbynuc(regionmask, intensity_image):
    tmp = intensity_image[regionmask]
    frequency = np.unique(tmp, return_counts=True)
    n_clusters = 26  # automate
    index = np.asarray(np.where(~np.in1d(np.arange(n_clusters) + 1, frequency[0])))
    for n in index[0]:
        frequency = np.insert(frequency, n, [0], axis=1)

    return frequency[1][1:]


marker_independent_properties = ('area',
                                 'bbox_area',
                                 'centroid',
                                 'convex_area',
                                 'eccentricity',  # [0-1]
                                 'equivalent_diameter',
                                 'extent',
                                 'feret_diameter_max',
                                 # longest distance between to points in convex hull ? = major axis length?
                                 'label',
                                 'major_axis_length',
                                 'minor_axis_length',
                                 'perimeter',
                                 'solidity')


def run_nuclear_features_table(well=43,
                               dir_input=Path(data_path, 'bg_subtracted'),
                               dir_nuclei=Path(data_path, 'segmented_nuclei'),
                               dir_MTU_img=Path(data_path, 'MTU_results'),
                               dir_output=Path(data_path, 'feature_tables'),
                               ):


    if not os.path.isfile(Path(dir_output, str(well) + '_feature_table.csv')):
        print('extracting nuclei features ' + str(well))

        nuclei = io.imread(Path(dir_nuclei, str(well) + '.tif'))

        props = regionprops_table(nuclei,
                                  properties=marker_independent_properties)

        marker_independent_property_table = pd.DataFrame(props)
        marker_independent_property_table['well'] = well

        with open(Path("metadata.yml"), 'r') as ymlfile:
            stains = yaml.load(ymlfile, Loader=yaml.FullLoader)

        img_df = get_metadata(dir_input)
        tmp = {}
        for stain in stains:
            filename = img_df[(img_df['stain'] == stain) & (img_df['well_id'] == well)]['file'].values[0]
            img = io.imread(Path(dir_input, filename))

            props = regionprops_table(nuclei, intensity_image=img, properties=['mean_intensity'],
                                      extra_properties=(percentiles, std, pixelcount, skewness, kurtosisness))

            internal_tmp = pd.DataFrame(props)
            numes = internal_tmp.columns
            newnames = [s + '_' + stain for s in numes]
            internal_tmp.columns = newnames

            tmp[stain] = internal_tmp

        marker_dependent_property_table = pd.concat(tmp, axis=1)
        marker_dependent_property_table.columns = marker_dependent_property_table.columns.droplevel()

        with open(Path("n_MTUs.yml"), 'r') as ymlfile:
            n_clusters = yaml.load(ymlfile, Loader=yaml.FullLoader)

        def clusterbynuc(regionmask, intensity_image, n_clusters=n_clusters):
            tmp = intensity_image[regionmask.astype('bool')]
            frequency = np.unique(tmp, return_counts=True)
            index = np.asarray(np.where(~np.in1d(np.arange(n_clusters + 1), frequency[0])))
            for n in index[0]:
                frequency = np.insert(frequency, n, [0], axis=1)

            return frequency[1][1:]

        labelimg = io.imread(Path(dir_MTU_img, str(well), 'labelimg.tif'))
        props = regionprops_table(nuclei, intensity_image=labelimg, extra_properties=(clusterbynuc,))
        clustercount = pd.DataFrame(np.vstack(props['clusterbynuc']))
        MTUs = np.arange(n_clusters) + 1
        string = 'MTU_count_'
        colnames = [string + str(s) for s in MTUs]
        clustercount.columns = colnames

        mask = io.imread(Path(data_path, 'masks', 'cropped_refined_mask' + str(well) + '.tif')).astype('bool')
        dist = ndimage.distance_transform_edt(mask)

        props = regionprops_table(nuclei, intensity_image=dist,
                                  properties=['mean_intensity'])
        radial_distance = pd.DataFrame(props)
        radial_distance.columns = ['radial distance']
        radial_distance['organoid_size'] = sum(sum(mask))

        property_table = pd.concat(
            [marker_independent_property_table, marker_dependent_property_table, clustercount, radial_distance], axis=1)

        os.makedirs(Path(dir_output), exist_ok=True)
        property_table.to_csv(str(Path(dir_output, str(well) + '_feature_table.csv')), sep=',', line_terminator='\n',
                              encoding="ISO-8859-1")
    else:
        print('nuclei features already exist' + str(well))


######above this line are functions that were used in the example processing. below not yet###########

def check_refined_mask(point, dir_mask='/links/groups/treutlein/DATA/imaging/charmel/masks',
                       dir_mask_refined='/links/groups/treutlein/DATA/imaging/charmel/refined_masks'):
    from skimage import io
    from matplotlib import pyplot as plt
    from pathlib import Path
    mask = io.imread(str(Path(dir_mask, 'mask' + str(point) + '.tif')))

    mask_refined = io.imread(str(Path(dir_mask_refined, 'mask' + str(point) + '.tif')))

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 10))

    ax[0].imshow(mask)
    ax[0].axis('off')
    ax[0].set_title('Mask')

    ax[1].imshow(mask_refined)
    ax[1].axis('off')
    ax[1].set_title('Mask refined')

    Path(dir_mask_refined, 'plots').mkdir(parents=True, exist_ok=True)
    plt.savefig(str(Path(dir_mask_refined, 'plots', str(point) + '.png')))


def show_image(image, scale=True):
    import matplotlib.pyplot as plt
    if scale:
        image = scale_image(image)
    plt.figure(figsize=(20, 20))
    plt.imshow(image)
    plt.show()


def show_staining(point, stain, dir_input='/links/groups/treutlein/DATA/imaging/PW/4i/plate14/', show_nucleus=False,
                  nucleus_cycle=1):


    channels = {
        "red": 0,
        "green": 1,
        "far_red": 3,
        "blue": 2
    }
    timepoints = pd.read_csv(
        '/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt')
    stainings = pd.read_csv("/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/ABs_table.txt", sep='\t')
    df_order = pd.read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/ID_ABorder_well.txt')
    df_order = df_order[df_order['ID'] == point]
    order = df_order['AB-order'].values[0]

    if order == 2:
        stainings = stainings.drop('cycle', 1).rename(columns={"permuted cycle": "cycle"})

    channel = stainings[stainings['AB'] == stain]['color channel'].values[0]

    cycle = 'cycle' + stainings[stainings['AB'] == stain]['cycle'].values[0].astype(str)

    # list images
    images = os.listdir(Path(dir_input, cycle, 'stitched'))
    images.sort()
    images = [image for image in images if "._" not in image]  # exclude .files that get listed sometimes
    # load image
    img = io.imread(str(Path(dir_input, cycle, 'stitched', images[point])))

    plt.figure(figsize=(20, 20))
    if show_nucleus:
        cycle = cycle + str(nucleus_cycle)
        img = scale_image(img[..., channels["blue"]])
        plt.suptitle(
            ' '.join(["Point:", str(point), "Stain:", "Hoechst", "from", cycle, timepoints.iloc[point].values[0]]),
            size=20)
    else:
        img = scale_image(img[..., channels[channel]])
        plt.suptitle(' '.join(["Point:", str(point), "Stain:", stain, "from", cycle, timepoints.iloc[point].values[0]]),
                     size=20)

    plt.imshow(img)
    plt.show()


def make_collage(point, stains, dir_input='/links/groups/treutlein/DATA/imaging/PW/4i/plate14/', nucleus_cycle=1,
                 crop=False, x=None, y=None, size=1000):


    # set up channel dictionary
    channels = {
        "red": 0,
        "green": 1,
        "far_red": 3,
        "blue": 2
    }

    # read in dataframes
    timepoints = pd.read_csv(
        '/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/conditions_per_sample.txt')
    stainings = pd.read_csv("/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/ABs_table.txt", sep='\t')
    df_order = pd.read_csv('/links/groups/treutlein/DATA/imaging/PW/4i/plate14/accessory_data/ID_ABorder_well.txt')
    df_order = df_order[df_order['ID'] == point]
    order = df_order['AB-order'].values[0]

    if order == 2:
        stainings = stainings.drop('cycle', 1).rename(columns={"permuted cycle": "cycle"})

    # set up plot
    fig = plt.figure(figsize=(18, 12))
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    fig.suptitle(' '.join(["Point:", str(point), 'at', timepoints.iloc[point].values[0]]), size=20)
    # iterate over stainings, load image and add to plot
    for i, stain in enumerate(stains):
        if stain == 'Hoechst':
            # select channel
            channel = 'blue'
            # select cycle
            cycle = 'cycle' + str(nucleus_cycle)
        else:
            # select channel
            channel = stainings[stainings['AB'] == stain]['color channel'].values[0]
            # select cycle
            cycle = 'cycle' + stainings[stainings['AB'] == stain]['cycle'].values[0].astype(str)
        print(stain, cycle)
        # list images
        images = os.listdir(Path(dir_input, cycle, 'stitched'))
        images.sort()

        # load image and get channel
        img = io.imread(str(Path(dir_input, cycle, 'stitched', images[point])))
        img = scale_image(img[..., channels[channel]])
        if crop:
            img = img[y - size:y + size, x - size:x + size]
        # add image to plot
        plt.subplot(2, 3, i + 1)
        plt.imshow(img)
        plt.axis('off')
        plt.title(stain)

    # return plot
    plt.show()


def overlay_channels(img, img_2, img_3=None):
    img = img / 65535
    img_2 = img_2 / 65535
    if img_3 is None:
        img_3 = np.zeros(img.shape)
    else:
        img_3 = img_3 / 65535
    return np.dstack([img, img_2, img_3])


def check_alignment(point, cycles, dir_input='/links/groups/treutlein/DATA/imaging/charmel/aligned', show_plot=True,
                    save_plot=False, dir_output=None):


    imgs = [str(Path(dir_input, str(point), cycle + '.tif')) for cycle in cycles]
    img_1 = io.imread(imgs[0])
    img_2 = io.imread(imgs[1])
    if len(cycles) == 2:
        overlay = overlay_channels(scale_image(img_1[..., 2]), scale_image(img_2[..., 2]))
    else:
        img_3 = io.imread(imgs[2])
        overlay = overlay_channels(scale_image(img_1[..., 2]), scale_image(img_2[..., 2]), scale_image(img_3[..., 2]))

    plt.figure(figsize=(20, 20))
    plt.imshow(overlay)
    plt.suptitle('_'.join(cycles), size=20)
    if save_plot:
        if dir_output is None:
            print('Specify dir_output in order to save plot.')
        else:
            Path(dir_output, str(point)).mkdir(parents=True, exist_ok=True)
            plt.savefig(str(Path(dir_output, str(point), '_'.join(cycles) + '.png')))
    if show_plot:
        return plt.show()


def crop_image(img, mask):
    mask = mask > 0
    return img[np.ix_(mask.any(1), mask.any(0))]
