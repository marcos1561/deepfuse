# Add the path were utils.py is
import sys
sys.path.append('./lib')

import os
import utils

u = utils.u

# Path to the current working directory
cwd = os.path.join(os.getcwd())
print(cwd)

''' Setup '''
# main directory that will containg everything
main_dir = os.path.join(cwd, "test_region") # Change this to save in other location.
os.makedirs(main_dir, exist_ok=True)

ims_dir = os.path.join(main_dir, "images") # where images are saved
ims_to_download_file = os.path.join(main_dir, "ims_to_download.csv") # file contaning the images to download
ims_id_file = os.path.join(main_dir, "ims_id.csv") # file with images already correctly downloaded.

# Images to download setup
ps = 0.7 * u.arcsec
# length = ps.value/(60*60)*(256*5) * u.deg # Region length
# height = ps.value/(60*60)*(256*4) * u.deg # Region hight
length, height = utils.find_angular_shape(3, 3, ps)
center_pos = (170.5717 * u.deg, 14.1776 * u.deg ) # (ra, dec) - Region center

# Cutouts setup
# NOTE: All elements from "kernel_shape" must be > 1.
kernel_shape = (2, 2) # num lines, num cols. Shape of the images to generate the mosaic for cutouts 
cutout_size = int(0.5 * utils.LEGACY_IM_SIZE) # pixels
overlap_percentage = 0.1
slice = 0


''' Creating the file that will containing the images to download '''
# Generates images center equatorial coordiantes
ims_eq_coords, ims_idx, region_shape = utils.generate_equatorial_coords(center_pos, length, height, ps)

# Check if the operation should continue
continue_operation = utils.should_continue(region_shape)
if continue_operation:
    # Generates the file containing the information to download the images
    utils.generate_sky_region_file(ims_eq_coords, ims_idx, ps, ims_to_download_file)

    ''' Downloading images '''
    download_error = utils.download_images(ims_to_download_file, ims_dir, ims_id_file, survey="ls-dr9")

    if not download_error:
        ''' Creating cutouts from the images downloaded'''
        utils.walk_through_region(kernel_shape, region_shape, cutout_size, overlap_percentage, slice, ims_dir, main_dir)
