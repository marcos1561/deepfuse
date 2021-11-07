import os
import functions, fetch
u = functions.u

# Path to the current working directory
cwd = os.path.join(os.getcwd())
print(cwd)

''' Setup '''
# main directory that will containg everything
main_dir = os.path.join(cwd, "sky_region/test_region") # Change this to save in other location.
os.makedirs(main_dir, exist_ok=True)

ims_dir = os.path.join(main_dir, "images") # where images are saved
ims_to_download_file = os.path.join(main_dir, "ims_to_download.csv") # file contaning the images to download
ims_id_file = os.path.join(main_dir, "ims_id.csv") # file with images already correctly downloaded.

# Images to download setup
ps = 0.7 * u.arcsec
length = ps.value/(60*60)*(256*4) * u.deg # Region length
height = ps.value/(60*60)*(256*5) * u.deg # Region hight
center_pos = (210.8002 * u.deg, 54.3478 * u.deg ) # (ra, dec) - Region center

# Cutouts setup
kernel_shape = (3, 3) # num lines, num cols. Shape of the images to generate the mosaic for cutouts 
cutout_size = int(1 * functions.LEGACY_IM_SIZE) # pixels
overlap_percentage = 0.1
slice = 0


''' Creating the file that will containg the images to download '''
# Generates images center equatorial coordiantes
ims_eq_coords, ims_idx, region_shape = functions.generate_equatorial_coords(center_pos, length, height, ps)

# Generates the file containg the information to download the images
functions.generate_sky_region_file(ims_eq_coords, ims_idx, ps, ims_to_download_file)


''' Downloading images '''
# Try to download the images until every image in the file was downloaded
download_error = True
while download_error:
    download_error = fetch.fetch(ims_to_download_file, ims_dir, ims_id_file, survey="ls-dr9")
    if download_error:
        print(download_error)
        print("\n\n-------- Not every image was downloaded --------\n")


''' Creating cutouts from the images downloaded'''
functions.walk_throught_region(kernel_shape, region_shape, cutout_size, overlap_percentage, slice, ims_dir, main_dir)
