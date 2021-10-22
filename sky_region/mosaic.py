import functions as func
import fetch
os = func.os

# Path to the directory of this file
cd = os.path.join(os.getcwd())
print(cd)

# User configuration #####################

# main directory that will containg everything
main_dir = os.path.join(cd, "sky_region/test_region")
os.makedirs(main_dir, exist_ok=True)

# Configurations for the sky region to download 
central_pos = (210.8002, 54.3478) # degre, degree
height = 0.7/(60*60)*(256*2) # degree
length = 0.7/(60*60)*(256*4) # degree
ps = 0.7 # arcsec

# Cutouts configurations 3 
cutout_size = 200
overlap_percentage = 0.1
###########################################


# Paths ##################################

file_path = os.path.join(main_dir, "ims_to_downlaod.csv") # file with the images to download
mosaic_ids_file = os.path.join(main_dir, "mosaic_ids.csv") # file informing wich images were downloaded correctly
mosaic_im_dir = os.path.join(main_dir, "mosaic_ims") # directory tha will contain the images downloaded
mosaic_seg_dir = os.path.join(main_dir, "mosaic_segments") # directory with mosaic segments
###########################################


# Download sky region #####################

# Generating file with the images to download
func.generate_sky_region_file(central_pos, length, height, ps, file_path)

# Downloading the images

# TODO: Make sure every image was downloaded
fetch.fetch(file_path, mosaic_im_dir, mosaic_ids_file, "ls-dr9")
###########################################


# Generating mosaic from the images downloaded ###

array, footprint, wcs_out, shape_out = func.generate_mosaic(mosaic_im_dir, 0)
###########################################


# Generate segmets from de mosaic created ###

# Creating segments
nrows, ncols = func.generae_mosaic_segments(array, wcs_out, shape_out, mosaic_seg_dir, cutout_size, overlap_percentage)
###########################################


# Visualing mosaic ###

# All mosaic
func.visualize_mosaic(array, footprint)

# Segments 
func.visualize_mosaic_segments(mosaic_seg_dir, nrows, ncols, array, sharexy=False)

func.plt.show()
###########################################
