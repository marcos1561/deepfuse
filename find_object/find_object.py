import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

def visualize_image(data):
    '''Plot mosaic segments for debugging purposes
    '''

    # Modifyng the scale for better visualization
    std = data.std()
    mean = data.mean()
    data = (data - mean + 1)**2

    plt.imshow(data, origin="lower", cmap="gray", vmin=1-1*std, vmax=1+1*std)
    plt.show()


LEGACY_IM_SIZE = 256

# Setup used to genereta the cutouts
ps = 0.7 * u.arcsec
length = ps.value/(60*60)*(256*5) * u.deg # Region length
height = ps.value/(60*60)*(256*4) * u.deg # Region hight
center_pos = (170.5717 * u.deg, 14.1776 * u.deg )
cutout_size = int(1 * LEGACY_IM_SIZE) # pixels

# Object equatorial coordinates
obj_ra = 170.6212 * u.deg
obj_dec = 14.1713 * u.deg

# Lower left equatorial coordinate of the region
lower_left_coord = (center_pos[0] + length/2, center_pos[1] - height/2) # (ra, dec)
cutout_angular_size = cutout_size * ps

# Based on the object coordinates, calculates cutout grid coordinates that probably have the object.
x_pos = (lower_left_coord[0] - obj_ra).to(u.rad) / cutout_angular_size.to(u.rad)
y_pos = (obj_dec - lower_left_coord[1]).to(u.rad) / cutout_angular_size.to(u.rad)
x_pos, y_pos = int(x_pos.value), int(y_pos.value)

'''
    Iterate over every neighbor cutout of the cutout find above and generate an object cutout if it's center
    coordinte is in the cutout.
'''
obj_cut_size = int(cutout_size * 0.3)
obj_skycoors = SkyCoord(ra=obj_ra, dec=obj_dec, frame='icrs')

for i in range(x_pos-1, x_pos+2):
    for j in range(y_pos-1, y_pos+2):
        print(f"\nImage: {i}_{j}.fits")
        im_path = f"./find_object/cutouts/{i}_{j}.fits"
        try:
            with fits.open(im_path) as im_file:
                im_wcs = WCS(im_file[1].header)
                im_data = im_file[1].data 
                # visualize_image(im_data)
                
                try:
                    obj_cut = Cutout2D(im_data, position=obj_skycoors, size=obj_cut_size, wcs=im_wcs, copy=True)
                    visualize_image(obj_cut.data)
                except ValueError as e:
                    print(e)

        except FileNotFoundError as e:
            print(e)