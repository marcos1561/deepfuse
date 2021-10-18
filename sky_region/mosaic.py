import functions as func
import fetch
os = func.os

# Path to the directory of this file
cd = os.path.join(os.getcwd())
print(cd)

# Download sky region #####################
# Configurations for the sky region to download 
central_pos = (210.8002, 54.3478) # degre, degree
height = 0.7/(60*60)*(256*2) # degree
length = 0.7/(60*60)*(256*4) # degree
ps = 0.7 # arcsec

file_path = os.path.join(cd, "ims_to_downlaod.csv")
mosaic_im_dir = os.path.join(cd, "mosaic_ims")
mosaic_ids_file = os.path.join(cd, "mosaic_ids.csv")


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
# Configurations 
cutout_size = 200
overlap_percentage = 0.1
mosaic_seg_dir = os.path.join(cd, "mosaic_segments")

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

