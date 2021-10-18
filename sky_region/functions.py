from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from astropy.wcs import WCS
import numpy as np
from math import ceil
import os
from astropy.nddata import Cutout2D

LEGACY_IM_SIZE = 256 # Size of the images downloaded from legacy

def generate_sky_region_file(central_pos, length, height, ps, file_path):
    print("\nGenerating files with the images to download...")

    # Calculating the number of pixel in both axis (They must be a multiple os LEGACY_IM_SIZE). 
    # It's used the linear aproximation in the tangence point to find the number of pixels, in other words, it's used the pixel scale.
    naxis1 = ceil(length / (ps/(60*60)))
    naxis2 = ceil(height / (ps/(60*60)))
    naxis1 = LEGACY_IM_SIZE * (naxis1//LEGACY_IM_SIZE) 
    naxis2 = LEGACY_IM_SIZE * (naxis2//LEGACY_IM_SIZE)

    # Creating WCS
    wcs_area = WCS(naxis=2)
    wcs_area.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs_area.wcs.crval = central_pos
    wcs_area.wcs.crpix = [naxis1/2, naxis2/2]
    wcs_area.wcs.cd = [[-ps/(60*60), 0], [0, ps/(60*60)]]
    wcs_area.wcs.cunit = ["deg", "deg"]
    print(wcs_area, "#"*10)

    # Calculating the number os images in both axis.
    num_im_1 = naxis1 // LEGACY_IM_SIZE
    num_im_2 = naxis2 // LEGACY_IM_SIZE

    print(f"\nMosaic image shape (nrows, ncols) = {num_im_2, num_im_1}")
    print(f"Mosaic pixel shape (nrows, ncols) = {naxis2, naxis1}")

    # Finding the center pixel coordinate of the images (Each image is LEGACY_IM_SIZE X LEGACY_IM_SIZE pixels).
    ims_central_coord = [[], []]
    for i in range(num_im_1):
        for j in range(num_im_2):
            axis1_coord = (LEGACY_IM_SIZE * (1/2 + i)) 
            axis2_coord = (LEGACY_IM_SIZE * (1/2 + j))
            ims_central_coord[0].append(axis1_coord)
            ims_central_coord[1].append(axis2_coord)

    # Finding the wrold coordinates of the central pixel coordinates and saving a .csv file in the for to downlaod the images
    ims_world_coord = wcs_area.pixel_to_world(*ims_central_coord)
    ims_world_coord = [(wc.ra.value, wc.dec.value, ps, "None", 0) for wc in ims_world_coord]

    pd.DataFrame(ims_world_coord, columns=("ra","dec","ps","host","label")).to_csv(file_path)


def generate_mosaic(filename_dir, slice):
    print("\nGenerating mosaic...")
    
    # Load the images
    mosaic_im = []
    for im_name in os.listdir(filename_dir):
        im_file = fits.open(os.path.join(filename_dir, im_name))
        im_file[0].data = im_file[0].data[slice]
        mosaic_im.append(im_file)

    print("Finding optimal wcs..")
    wcs_out, shape_out = find_optimal_celestial_wcs(mosaic_im)

    print("Creating mosaic...")
    array, footprint = reproject_and_coadd(mosaic_im, wcs_out, shape_out=shape_out, reproject_function=reproject_interp)

    # Closing open files
    for im in mosaic_im:
        im.close()
    del mosaic_im

    return array, footprint, wcs_out, shape_out


def generae_mosaic_segments(array, wcs_out, shape_out, dest_dir, cutout_size, overlap_percentage):
    print("\n Generating mosaic segments...")

    orig_header = wcs_out.to_header()

    # Dict to store segmets positions for later use in plotting.
    mosaic_pos = {"x": [], "y": [], "name": []}

    os.makedirs(dest_dir,  exist_ok=True)
    cutout = (cutout_size, cutout_size) 
    overlap = cutout_size * (1-overlap_percentage) # if cutout_size = 500 px, setting overlap = 0.2 would leave 100 pixels overlapping 
    num_images_per_row = ceil(shape_out[1]/overlap)
    num_images_per_column = ceil(shape_out[0]/overlap)

    print("mosaic shape: ", shape_out)
    print("num_images_per_row: ", num_images_per_row)
    print("num_images_per_column: ", num_images_per_column)

    row_pos = ceil(cutout_size/2)

    for i in range(num_images_per_row): # changes the row
        col_pos = ceil(cutout_size/2)
        
        for j in range(num_images_per_column): # changes the column

            cen = (row_pos, col_pos)
            print("processing cutout at position ", str(cen))

            segment = Cutout2D(array, position=cen, size=cutout, wcs=wcs_out, copy=True)
            
            # checks if the segment has values other than 0, ignores the CCD parts with no good pixels
            if not np.all((segment.data == 0)): # returns True if all zeros, use "not" in front of it to look at the segments with other values only
                # Position for later use (plot all mosaic elements using matplotlib)
                mosaic_pos["name"].append(str(cen)+".fits")
                pos = (num_images_per_column-j-1, i) 
                mosaic_pos["x"].append(pos[0])
                mosaic_pos["y"].append(pos[1])
                
                header_new = segment.wcs.to_header()
            
                hdu_p = fits.PrimaryHDU(header=orig_header)
                hdu_i = fits.ImageHDU(segment.data, header=header_new)
                hdulist = fits.HDUList([hdu_p,hdu_i])

                output_file = os.path.join(dest_dir, str(cen)+".fits")
                hdulist.writeto(output_file, overwrite=True)

            col_pos += overlap
            
        row_pos += overlap

        # Save position for later use. (plot all mosaic elements using matplotlib)
        pd.DataFrame(mosaic_pos).to_csv(os.path.join(dest_dir, "mosaic_pos.csv"))

    return num_images_per_column, num_images_per_row


def visualize_mosaic(mosaic_array, footprint):
    std = mosaic_array.std()
    mean = mosaic_array.mean()
    
    mosaic_array = (mosaic_array-mean+1)**2

    fig = plt.figure()

    ax1 = plt.subplot(2, 1, 1)
    ax1.imshow(mosaic_array, cmap="gray", origin='lower', vmin=1-2*std, vmax=1+10*std)
    ax1.set_title('Mosaic')

    ax2 = plt.subplot(2, 1, 2)
    ax2.imshow(footprint, origin='lower')
    ax2.set_title('Footprint')


def visualize_mosaic_segments(mosaic_seg_path, nrows, ncols, mosaic_array, sharexy=False):
    # File containg the correct postion for plotting
    segmets_pos_df = pd.read_csv(os.path.join(mosaic_seg_path, "mosaic_pos.csv"))


    positions = zip(segmets_pos_df["x"].astype(int), segmets_pos_df["y"].astype(int))
    
    # Loading mosaic segments data
    mosaic_segmets_data = []
    for name in segmets_pos_df["name"]:
        im_data = fits.getdata(os.path.join(mosaic_seg_path, name))
        mosaic_segmets_data.append(im_data)

    # Midifyng the scale for better visualization
    std = mosaic_array.std()
    mean = mosaic_array.mean()
    mosaic_segmets_data = [(im_i - mean + 1)**2 for im_i in  mosaic_segmets_data]


    # Plotting segmets
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharexy, sharey=sharexy)

    for pos, data, name in zip(positions, mosaic_segmets_data, segmets_pos_df["name"]):
        ax[pos[0], pos[1]].imshow(data, origin="lower", cmap="gray", vmin=1-1*std, vmax=1+1*std)
        ax[pos[0], pos[1]].set_title(name)
