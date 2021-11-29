from genericpath import exists
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.utils.misc import isiterable

from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import floor, ceil
import os

LEGACY_IM_SIZE = 256 # Size of the images downloaded from legacy

def rotate_y(x, angle):
    '''Rotate points around the y axis.

    Parameters:
    -----------
    x: 1-d array, 2-d array
        The points to rotate. If it's 2-d array, the points coordinates must in the columns.
    
    angle: Astropy Quantaty
        Angle of rotation.

    Return:
    -------
    Return the points given in x rotated.
    '''
    angle = angle.to(u.rad).value
    matrix = np.array([[np.cos(angle), 0, -np.sin(angle)], [0, 1, 0], [np.sin(angle), 0, np.cos(angle)]])
    return np.dot(matrix, x)


def rotate_z(x, angle):
    '''Rotate points around the z axis.

    Parameters:
    -----------
    x: 1-d array, 2-d array
        The points to rotate. If it's 2-d array, the points coordinates must in the columns.
    
    angle: Astropy Quantaty
        Angle of rotation.

    Return:
    -------
    Return the points given in x rotated.
    '''
    angle = angle.to(u.rad).value
    matrix = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle),np.cos(angle), 0], [0, 0, 1]])  
    return np.dot(matrix, x)


def view_points(points):
    ''' Function for visualizing the center coordiantes of the images in a unit sphere 
        for debugging purposes.
    '''

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    if isiterable(points[0]):
        xyz_range = []
        xyz_middle = []
        for i in range(3):
            xyz_range.append(abs(max(points[i]) - min(points[i])))
            xyz_middle.append((max(points[i]) + min(points[i])) / 2)    
        radius = max(xyz_range) / 2

        prc_space = 0
        ax.set_xlim((xyz_middle[0]-radius)*(1-prc_space), (xyz_middle[0] + radius)*(1+prc_space))
        ax.set_ylim((xyz_middle[1]-radius)*(1-prc_space), (xyz_middle[1] + radius)*(1+prc_space))
        ax.set_zlim((xyz_middle[2]-radius)*(1-prc_space), (xyz_middle[2] + radius)*(1+prc_space))

    ax.plot([0, 1], [0, 0], [0,0], color="red")
    ax.plot([0, 0], [0, 1], [0,0], color="green")
    ax.plot([0, 0], [0, 0], [0,1], color="blue")

    ax.scatter(points[0], points[1], points[2])
    
    plt.show()


def generate_equatorial_coords(center_pos, length, height, ps_original):
    '''Generates the images center equatorial coordiantes in the specify region
    to download in legacy survey.
    First the points are created centered in the origin. The limits are 
    
    lenght/2 < ra < lenght/2, height/2 < dec < height/2

    Then, the points are rotated to the given center. 

    Parameters:
    -----------
        center_pos: (Astropy Quantaty, Astropy Quantaty)
            Center of the region in equatorial coordinates.

        lenght: Astropy Quantaty
            Angular lenght of the region.

        height: Astropy Quantaty
            Angular height of the region.

        ps: float
            pixel scale of the images to download.

    Return:
    ------
        eq_coords: np-array
            Array cointaining the images center equatorial coordinates.
        
        ims_idx: np-array 
            Array containg the images indexes.
        
        region_shape: (int, int)
            Tuple containing the number of images in each column and row, respectively.
    '''
    # Reduce pixel scale for preventing empty pixel between images
    ps = (LEGACY_IM_SIZE-2)/LEGACY_IM_SIZE * ps_original

    # Angular length (in rad) of the images to be downloaded
    img_length = ps.to(u.rad).value * LEGACY_IM_SIZE 
    
    lenght_rad = length.to(u.rad).value 
    heiht_rad = height.to(u.rad).value

    ### Grid of images center equatorial coordinates ###
    ra_range = np.arange(-lenght_rad/2 + img_length/2, lenght_rad/2, img_length)
    dec_range = np.arange(-heiht_rad/2 + img_length/2, heiht_rad/2, img_length)

    total_ims_per_col = dec_range.size
    total_ims_per_row = ra_range.size
    region_shape = (total_ims_per_col, total_ims_per_row)

    meshgrid = np.meshgrid(ra_range, dec_range)
    ra = meshgrid[0].flatten() # Array with all the RA positions.
    dec = meshgrid[1].flatten() # Array withh all the DEc positions.
    ######

    # Array containing images integer coordinates (origin is in the lower left corner)
    # The image in the ith column and jth row will have the coordinate (i, j)
    idx_meshgrid = np.meshgrid(np.arange(total_ims_per_row-1, -1, -1), np.arange(total_ims_per_col))
    img_idx = np.array([idx_meshgrid[0].flatten(), idx_meshgrid[1].flatten()])

    # Points in the unit sphere with cartesian coordinates correponding to the equatorial coordiantes.
    points = np.zeros(shape=(3, total_ims_per_col*total_ims_per_row))
    points_cartesian = SkyCoord(ra=ra * u.rad, dec=dec * u.rad).cartesian

    # Rotating points, such that the center of the region will go the center center_pos parameter.
    points[:] = np.array([points_cartesian.x, points_cartesian.y, points_cartesian.z])
    points = rotate_z(rotate_y(points, center_pos[1]), center_pos[0])
    # view_points(points)

    # Tranforming to equatorial coordinates 
    eq_coords = SkyCoord(points[0], points[1], points[2], representation_type="cartesian").galactic.transform_to("icrs")
    eq_coords = np.array([eq_coords.ra.value, eq_coords.dec.value])

    print(f"Num images to download: {eq_coords[0].size}")
    print(f"Num_cols: {total_ims_per_row} | Num_rows: {total_ims_per_col}")

    return eq_coords, img_idx, region_shape


def generate_sky_region_file(eq_coords, img_idx, ps, file_path):
    '''Create the file used to download the images.

    Parameters:
    -----------
        eq_coords: 1-d np array
            Array containing the images coordintes.
        
        img_idx: 1-d np-array
            Array containing the images indexes.
        
        ps: float
            Pixel scale of the images
        
        file_path: str
            Path where the file will be saved. 
    
    Return:
    -------
        None
    '''
    df_data = np.append(eq_coords, np.zeros((3, img_idx[0].size)), axis=0)
    index = [f"{img_idx[0][i]}_{img_idx[1][i]}" for i in range(img_idx[0].size)]
    columns = ("ra","dec","ps","host","label")

    ims_df = pd.DataFrame(np.transpose(df_data), columns=columns, index=index)
    ims_df["ps"] = ps.to(u.arcsec).value
    ims_df.to_csv(file_path)


def generate_mosaic(mosaic_im):
    '''Generates a mosaic of .fits images and return the mosaic data array.

    Parameters:
    -----------
        mosaic_im: 2-d list or 1-d list of ~astropy.io.fits.HDUList 
            The images that will form the mosaic. If a 2-d list is given, it is flattened.

    Return:
    -------
        array: 2d-array
            Array containing mosaic data.

        footprint: 2d-array
            Array containg the number of contribution of each pixel.

        wcs : ~astropy.wcs.WCS
            The optimal WCS determined from the input images. 
         
        shape : tuple
            The optimal shape required to cover all the output.
    '''
    print("\nGenerating mosaic...")
    
    # Flattening the list if necessary
    if type(mosaic_im[0]) == list: 
        mosaic_im = [item for sublist in mosaic_im for item in sublist]

    print("Finding optimal wcs..")
    wcs_out, shape_out = find_optimal_celestial_wcs(mosaic_im)

    print("Creating mosaic...")
    array, footprint = reproject_and_coadd(mosaic_im, wcs_out, shape_out=shape_out, reproject_function=reproject_interp)
    print("Mosaic created")

    return array, footprint, wcs_out, shape_out


def generate_mosaic_cutouts(array, wcs_out, shape_out, cuts, dest_dir, cutout_size, overlap_percentage, lower_left_grid_coord = (0, 0)):
    '''Generates cutouts with overlap in the mosaic provided by "array".

    Parameters:
        array: 2-d array
            Mosaic data tha will be used to generate the cutouts.
        
        wcs : ~astropy.wcs.WCS
            Mosaic WCS.

        shape_out: tuple
            Mosaic_shape (num_rows, num_cols)
        
        cuts: dic {"line": int, "col": int}
            First lines and columns to ignore when making the cutouts.

        dest_dir: str
            path to the directory to save the cutouts.

        cutout_size: int
            Side dimension of the cutout. For example: 1000 outputs a 1000x1000 cutouts.
        
        overlap_percentage: float
            Percentage of overlap between cutouts (0. = completely new cutout, no overlap; 1. = same cutout, total overlap).

        last_lower_left_grid_coord: tuple
            Last grid coordinte of the lower left cutout of this segment. This is usefull if we 
            have more than one segment.
        
        Return:
        -------
            num_images_per_column: int
                Number of cutouts per column.
            
            num_images_per_row: int
                Number of cutouts per row.
            
            cen: tuple
                Center position of upper right cutout.
    '''

    print("\nGenerating mosaic cutouts...")

    orig_header = wcs_out.to_header()

    # Dict to store segmets positions for later use in plotting.
    # mosaic_pos = {"row": [], "col": [], "name": []}

    cutout = (cutout_size, cutout_size) 
    overlap = cutout_size * (1-overlap_percentage) # if cutout_size = 500 px, setting overlap = 0.2 would leave 100 pixels overlapping 
    num_images_per_column = floor((shape_out[0] - cuts["line"])/overlap)
    num_images_per_row = floor((shape_out[1] - cuts["col"])/overlap)
    
    print("mosaic shape: ", shape_out)
    print("num_images_per_row: ", num_images_per_row)
    print("num_images_per_column: ", num_images_per_column)
    print(f"Mosaic cuts: {cuts}\n")

    # The coordinate origin is in the mosaic lower left corner 
    # and the segments center coordinates must be in the form (x, y)
    x_pos = cuts["col"] + ceil(cutout_size/2)
    for i in range(num_images_per_row): # changes the x coodr
        y_pos = cuts["line"] + ceil(cutout_size/2)
        for j in range(num_images_per_column): # changes the y coord
            cen = (x_pos, y_pos)
            print("processing cutout at position ", str(cen))

            segment = Cutout2D(array, position=cen, size=cutout, wcs=wcs_out, copy=True)
            
            # checks if the segment has values other than 0, ignores the CCD parts with no good pixels
            if not np.all((segment.data == 0)): # returns True if all zeros, use "not" in front of it to look at the segments with other values only
                seg_name = f"{i + lower_left_grid_coord[0]}_{j + lower_left_grid_coord[1]}.fits"

                # Position for later use (plot all mosaic elements using matplotlib)
                # mosaic_pos["name"].append(seg_name)
                # pos = (num_images_per_column -1 - j, i)
                # mosaic_pos["row"].append(pos[0])
                # mosaic_pos["col"].append(pos[1])
                
                header_new = segment.wcs.to_header()
            
                hdu_p = fits.PrimaryHDU(header=orig_header)
                hdu_i = fits.ImageHDU(segment.data, header=header_new)
                hdulist = fits.HDUList([hdu_p,hdu_i])

                output_file = os.path.join(dest_dir, seg_name)
                hdulist.writeto(output_file, overwrite=True)

            y_pos += overlap
            
        x_pos += overlap

        # Save position for later use. (plot all mosaic elements using matplotlib)
        # pd.DataFrame(mosaic_pos).to_csv(os.path.join(dest_dir, "mosaic_pos.csv"))
    
    if num_images_per_column > 0 and num_images_per_row > 0:
        print("Cutouts Created!")
    else:
        print("It was not possible to generate cutouts in this segment"
              "(probably cutout_size is bigger than some mosaic dimention)\n")

    return num_images_per_column, num_images_per_row, cen


def visualize_mosaic(mosaic_array, footprint):
    '''Plot mosaic and footprint for debbuging purposes
    '''

    # std = mosaic_array.std()
    # mean = mosaic_array.mean()
    
    # mosaic_array = (mosaic_array-mean+1)**2

    fig = plt.figure()

    ax1 = plt.subplot(2, 1, 1)
    # ax1.imshow(mosaic_array, cmap="gray", origin='lower', vmin=1-2*std, vmax=1+10*std)
    ax1.imshow(mosaic_array, cmap="gray", origin='lower')
    ax1.set_title('Mosaic')

    ax2 = plt.subplot(2, 1, 2)
    ax2.imshow(footprint, origin='lower')
    ax2.set_title('Footprint')


def visualize_mosaic_cutouts(mosaic_seg_path, nrows, ncols, mosaic_array, sharexy=False):
    '''Plot mosaic segments for debugging purposes
    '''
    if nrows == 0 or ncols == 0:
        return

    # File containg the correct postion for plotting
    segmets_pos_df = pd.read_csv(os.path.join(mosaic_seg_path, "mosaic_pos.csv"))


    positions = zip(segmets_pos_df["row"].astype(int), segmets_pos_df["col"].astype(int))
    
    # Loading mosaic segments data
    mosaic_segmets_data = []
    for name in segmets_pos_df["name"]:
        im_data = fits.getdata(os.path.join(mosaic_seg_path, name))
        mosaic_segmets_data.append(im_data)

    # Midifyng the scale for better visualization
    # std = mosaic_array.std()
    # mean = mosaic_array.mean()
    # mosaic_segmets_data = [(im_i - mean + 1)**2 for im_i in  mosaic_segmets_data]


    # Plotting segmets
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharexy, sharey=sharexy)
    
    if nrows == 1 and ncols == 1:
        ax = np.array([ax])

    if ax.ndim == 1:
        ax = ax.reshape(nrows, ncols)

    for pos, data, name in zip(positions, mosaic_segmets_data, segmets_pos_df["name"]):
        # ax[pos[0], pos[1]].imshow(data, origin="lower", cmap="gray", vmin=1-1*std, vmax=1+1*std)
        ax[pos[0], pos[1]].imshow(data, origin="lower", cmap="gray")
        # ax[pos[0], pos[1]].set_title(name)


def walk_throught_region(kernel_shape, region_shape, cutout_size, overlap, slice, ims_dir, main_dir):
    '''Walks between segments of the region formed by the images in "ims_dir" and
    creates cutouts with overlap in each segment. 
    
    Segments are formed by a group of images of the region. The shape of the segmets are
    specify by "kernel_shape". It will only form "cutout_size" X "cutout_size" cutouts. 
    
    Start taking the segmente in the lower left corner, next take the segment above and so on.
    When the top of the row is reached, a smaller segmegment is created if possible, 
    and when the row is finished, starts the same process in the next row. 

    The cutouts are saved in "main_dir", such that the coutouts of each segment are 
    in different directorys, named by the position of the lower left image of the segment.
    
    The name of the images must be in the form "i_j.fits", where i is the image column position
    and j the image row position, given that the images are organized in a grid and the origin 
    is in the lower left corner.

    Parameters:
    -----------
        kernel_shape: (int, int)
            The shape of the segmets in images, for example, (2, 3) will generate segments
            of 2 images in length and 3 images in hight.
        
        region_shape: (int, int)
            List containing the number of images per column and row in the region, respectively.
        
        cutout_size: int
            Size of the cutouts.
        
        overlap: float
            Percentage of overlap.
        
        slice: int
            The slice to use in the images (They can have three channels)
        
        ims_dir: str
            Path to the directory containg the images that form the region.

        main_dir: str
            Path to the directory where the cutout will be saved.

    Return:
    -------
        None
    '''
    print("\n####### Walking throught images #########")
    total_ims_per_col = region_shape[0]
    total_ims_per_row = region_shape[1]



    # This loop walks between the segmetns. col_i and row_i containg the grid coordinate (for 
    # example (0,3) is the image in the first column and the 4th row) of the lower
    # left image of the segment.
    #
    # line_cut and col_cut are the limits of the segment that are not necessary for the cutouts 
    # ( it was already used in previous segments). For example, line_cut = 10, col_cut = 50 will 
    # cut the first 10 lines and 50 columns. (origin is in lower left corner)
    #
    # Since there is overlap, we need to maintein some rows of the last segment and some
    # columns of the left neighbor segment. This information is stored in num_col_keep and num_row_keep
    
    # This is used to correctly name the cutouts
    lower_left_grid_coord = (0,0) 
    
    # Path to store the segments
    mosaic_seg_dir = os.path.join(main_dir, f"cutouts") 
    if os.path.exists(mosaic_seg_dir):
        for f in os.listdir(mosaic_seg_dir): # Removes existing cutouts if it exists
            os.remove(os.path.join(mosaic_seg_dir, f))
    else:
        os.makedirs(mosaic_seg_dir)

    col_i = 0
    num_col_keep = 0 
    col_cut = 0
    while col_i < total_ims_per_row - num_col_keep:
        if col_i != 0:
            new_col = True # This variable is used to correctly name the cutouts
        else:
            new_col = False

        mosaic_ims = [] # List containing the images of the current segment
        ims_names_debug = [] # List containing the images grid coordinates for debugging purposes
        
        row_i = 0
        num_row_keep = 0
        line_cut = 0
        while row_i < total_ims_per_col - num_row_keep:
            print(f"\n#### Mosaic lower left corner: {col_i}, {row_i} ###")
            print(f"\nNum_row_keeped: {num_row_keep} | Num_col_keeped: {num_col_keep}")
            
            # Get the names of the images that will form the current segment. The names are
            # stored in a grid (2d array), equal to the grid they form when combined to generate the segment.
            ims_names_list = []
            ims_names_debug = ims_names_list[-num_row_keep:]
            for i in range(row_i + num_row_keep, row_i + kernel_shape[0]):
                ims_names_row = []
                for j in range(col_i, col_i + kernel_shape[1]):
                    if i > (total_ims_per_col-1) or j > (total_ims_per_row-1): # Segments in the top and right limits can be smaller
                        continue
                    ims_names_row.append(f"{j}_{i}.fits")
                ims_names_list.append(ims_names_row)
            
            ims_names_debug = [*ims_names_debug, *ims_names_list]
            print(f"Ims: {ims_names_debug}")

            # Loads the images that will form the current segment. The images are stored
            # in the same structure of "ims_names_list"
            # Keeps rows from last segment that will be needed for this segment 
            mosaic_ims = mosaic_ims[-num_row_keep:]
            for ims_row in ims_names_list:
                mosaic_ims_row = []
                for im_name in ims_row:
                    with fits.open(os.path.join(ims_dir, im_name), memmap=False) as im_file:
                        im_file[0].data = im_file[0].data[slice]
                        mosaic_ims_row.append(fits.HDUList(im_file[0].copy()))
                mosaic_ims.append(mosaic_ims_row)
            

            # Based on the upper right cutout generated in the last mosaic and the last mosaic shape, calculates the rows and
            # columns tha will not be needed for this mosaic 
            if row_i > 0: # line cut
                line_cut = last_up_right_pos[1] + cutout_size*(0.5 - overlap) - last_shape_out[0] + num_row_keep * LEGACY_IM_SIZE

            # Since mosaics are created from the bottom up, giving them coordinates ( 0,0 is the lower left mosaic ), 
            # column cuts of mosaic that are in column n are calculated based on the last mosaic in column n-1.
            if col_i > 0 and row_i == 0: # col cut
                col_cut = last_up_right_pos[0] + cutout_size*(0.5 - overlap) - last_shape_out[1] + num_col_keep * LEGACY_IM_SIZE
            
            cuts = {"line":line_cut, "col":col_cut}

            # Combine segment's images to generate a mosaic
            array, footprint, wcs_out, shape_out = generate_mosaic(mosaic_ims)
            last_shape_out = shape_out
            
            # Generate cutouts
            # mosaic_seg_dir = os.path.join(main_dir, f"cutout_dir/{col_i}_{row_i}") # mosaic cutouts path
            
            # Recalculates the grid coordinate of the lower left cutout of this segment.
            # It only needs to recalculate if it's not the first time generating cutouts. This is
            # done checking if we are not in the first line and first column of the images. 
            if row_i != 0 or col_i != 0:
                if new_col:
                    x_coord = last_lower_left_grid_coord[0] + nrows
                    y_coord = 0
                else:
                    x_coord = last_lower_left_grid_coord[0]
                    y_coord = last_lower_left_grid_coord[1] + ncols
                lower_left_grid_coord = (x_coord, y_coord)

            ncols, nrows, up_right_pos = generate_mosaic_cutouts(array, wcs_out, shape_out, cuts, mosaic_seg_dir, cutout_size, overlap, lower_left_grid_coord)
            
            last_lower_left_grid_coord = lower_left_grid_coord

            # TODO: It's possible to ncols our nrows be 0 when cutout_size > LEGACY_IM_SIZE. Fix this pls.   
            if not(ncols == 0 or nrows == 0):
                last_up_right_pos = up_right_pos

            # TODO: Some pixels are missing from the edge of the images, a possible solution 
            # would be to make the equatorial coordinates of the center of the images a little bit closer.
            # visualize_mosaic(array, footprint)
            # visualize_mosaic_cutouts(mosaic_seg_dir, ncols, nrows, array)

            plt.show()

            # Number of rows to keep in the next segment 
            # upper_cutout_y is the y position (in this segment) from lower edge of the next cutout 
            # above the upper right cutout in this segment
            if not(ncols == 0 or nrows == 0):
                upper_cutout_y = up_right_pos[1] + cutout_size*(0.5 - overlap)
                num_row_keep = ceil((shape_out[0] - upper_cutout_y) / LEGACY_IM_SIZE)
                row_i += kernel_shape[0] - num_row_keep
            else:
                row_i += kernel_shape[0]

            # We no longer are in a new column
            new_col = False
        
        # Number of columns to keep in the next segment 
        # right_cutout_x is the x position (in this segment) from the left edge 
        # of the next cutout to the right of the top right cutout in this segment
        if not(ncols == 0 or nrows == 0):
            right_cutout_x = up_right_pos[0] + cutout_size*(0.5 - overlap)
            num_col_keep = ceil((shape_out[1] - right_cutout_x) / LEGACY_IM_SIZE)
            col_i += kernel_shape[1] - num_col_keep
        else:
            col_i += kernel_shape[1]
