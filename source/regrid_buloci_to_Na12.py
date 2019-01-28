import scipy.interpolate
import numpy as np
import xarray as xr
import os

import matplotlib.pyplot as plt

from regrid import get_ease_coords

def read_mask():
    """
    Returns Northern hemisphere high resolution LOCI mask
    """
    
    diri = '/oldhome/apbarret/projects/ancillary/masks'
    fili = 'Nh_loci_land50_coast0km.1441x1441.bin'

    xdim = 1441
    ydim = 1441
    
    mask = np.fromfile( os.path.join(diri, fili), dtype='uint8').reshape(xdim,ydim)

    # Set mask to 0 for land [0] (and ice [101]) and 1 for ocean [255].  Cells off-globe
    # are set to NaN
    newmask = np.where( (mask == 0) | (mask == 101), 0, np.NaN )
    newmask = np.where( (mask == 255), 1, newmask )
    
    return newmask

def get_input_mask():
    """
    Returns points (npoints,2) and values (npoints) arrays for BU LOCI mask.
    """

    mask = read_mask()

    coords = get_ease_coords('Nh')

    values = mask[ np.isfinite(mask) ].flatten().astype('uint8')
    points = np.array( [coords[0][ np.isfinite(mask) ].flatten(),
                        coords[1][ np.isfinite(mask) ].flatten()] ).T
    
    return points, values
    
def main(verbose=True):

    if verbose: print ('Getting mask and coordinates')
    points, values = get_input_mask()

    if verbose: print ('Getting EASE grid definition')
    dstCoord = get_ease_coords('Na12')

    if verbose: print ('Regridding mask')
    mask = scipy.interpolate.griddata(points, values, dstCoord, method='nearest')

    np.save('BU_loci_Na12.npy', mask)

if __name__ == "__main__":
    main()
