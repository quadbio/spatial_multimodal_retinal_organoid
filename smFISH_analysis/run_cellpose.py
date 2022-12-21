from cellpose import models
from cellpose.io import imread
import os
from pathlib import Path
from tqdm import tqdm
from skimage import io
model = models.Cellpose(model_type='cyto2')

dir_images = '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission/DAPI-lower exposure_submission'
dir_output = '/links/groups/treutlein/DATA/imaging/charmel/32955-slide928-2_submission/segmented'
Path(dir_output).mkdir(exist_ok=True, parents=True)
files = os.listdir(dir_images)
files = [file for file in files if '._' not in file]

for file in tqdm(files):
    img = imread(Path(dir_images,file))
    mask, flows, styles, diams = model.eval(img, diameter=50, channels=[[0,0]], do_3D=False,
                                             flow_threshold=-3, cellprob_threshold=0.8)
    io.imsave(Path(dir_output,file),mask,check_contrast=False)




