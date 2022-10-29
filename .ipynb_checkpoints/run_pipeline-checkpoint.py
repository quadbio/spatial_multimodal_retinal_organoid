
import yaml 
import importlib
import modules 
from time import gmtime, strftime
import os
from pathlib import Path
import subprocess

#load global variables and parameters
with open("params.yml", 'r') as ymlfile: 
   cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

globals().update(cfg)
r_path = str(Path('/usr/local/R-4.0.3/bin/Rscript')) # path to r interpreter
#################-----------------------------------------------------------------########################
for well in wells:
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    print('simple masking sample:' + well)
    modules.simple_mask(point=well, cycle=ref_cycles[0],
                dir_input=data_path,
                dir_output=data_path + 'masks/',
                save=True,
                plot=False,
                show_plot=False,
                save_plot=True)

    
    img_df = modules.get_metadata(Path(data_path,'raw_data'))
    filenames = img_df.loc[(img_df['well_id'] == well)]['file']
    
    items = []
    for i in filenames:
        items.append(Path(data_path,'aligned', i)) 
    
    print('aligning well '+ well)
    if not (all(os.path.isfile(i) for i in items)):
        modules.run_elastix(point = well,
                    ref_cycles = ref_cycles,
                    dir_mask = Path(data_path,'masks'),
                    dir_input = Path(data_path,'raw_data'),
                    dir_output = Path(data_path,'aligned'))
    else:
        print('well ' + str(well) + ' already aligned')

    print('denoising well '+ well)
    modules.denoise_well(well = well)
    
    
    print('refigning mask well '+ well)
    modules.refine_mask(well = well)
    
    
    img_df = modules.get_metadata(Path(data_path,'denoised'))
    filenames = img_df.loc[(img_df['well_id'] == well)]['file']
    
    items = []
    for i in filenames:
        items.append(Path(data_path,'bg_subtracted', i)) 
        
    print('subtracting backgrounds well '+ well)
    if not (all(os.path.isfile(i) for i in items)):
        modules.bg_subtraction(well = well,
                        dir_input=Path(data_path,'denoised'),
                        dir_masks=Path(data_path,'masks'),
                        dir_output=Path(data_path,'bg_subtracted'),
                        cycles_bg=[0, 6, 12],
                        cycle_hoechst = 1)
    
    print('segmenting nuclei well '+ well)
    modules.segment_nuclei(well = well, 
                           ref_cycle = 0, 
                           dir_input = Path(data_path, 'denoised'), 
                           dir_output = Path(data_path, 'segmented_nuclei'))
   
    print('generating pixel matrix well '+ well)
    modules.generate_pixel_matrices(well = well,
                         dir_input = Path(data_path,"bg_subtracted"),
                         dir_output = Path(data_path,'pixel_matrices/'),
                         cycles_bg = cycles_bg,
                         to_collagen_mask = to_collagen_mask,
                         mask_collagen = False,
                         overwrite_denoised_images_with_collagen_masked_images = False)

conditions = {'week39' : [43,44]}
for condition in conditions:
    modules.run_matrix_normalization(wells, 
                             dir_input = Path(data_path, 'pixel_matrices'),
                             dir_output = Path(data_path,'pixel_matrices'), 
                             fusion_matrix = True)


if not os.path.isfile(Path(data_path, 'fSOM_output/fSOM.rds')):
                      subprocess.run([r_path,'fSOM.R'])

for well in wells:
    modules.run_nuclear_features_table(well = well,
                             dir_input = Path(data_path,'bg_subtracted'),
                             dir_nuclei = Path(data_path,'segmented_nuclei'),
                             dir_MTU_img = Path(data_path, 'MTU_results'),
                             dir_output = Path(data_path,'feature_tables'),
                             )