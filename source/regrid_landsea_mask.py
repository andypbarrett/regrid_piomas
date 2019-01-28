import scipy.interpolate
import numpy as np
import xarray as xr
import os

from regrid import get_ease_coords

def read_mask():
    """
    Returns points (npoints,2) and values (npoints) arrays for land mask.
    Mask is subsetted to north of 40 N and for every other point to improve speed
    """
    
    diri = '/disks/arctic5_raid/abarrett/ETOPO1'
    fili = 'etopo1_land_ocean_mask.nc'

    ds = xr.open_dataset(os.path.join(diri,fili))
    tmp = ds['__xarray_dataarray_variable__'].values
    lat = ds['lat'].values
    lon = ds['lon'].values
    ds.close()
    
    xs = tmp[lat > 40., :]
    values = xs[::2,::2].flatten()

    lats = lat[lat > 40.][::2]
    lons = lon[::2]

    x, y = np.meshgrid(lons, lats)
    points = np.array([x.flatten(), y.flatten()]).T

    return (points, values)

def main(verbose=True):

    if verbose: print ('Getting mask and coordinates')
    points, values = read_mask()

    if verbose: print ('Getting EASE grid definition')
    dstCoord = get_ease_coords('Na12')

    if verbose: print ('Regridding mask')
    mask = scipy.interpolate.griddata(points, values, dstCoord, method='nearest')

    np.save('landsea_mask_Na12.npy', mask)

if __name__ == "__main__":
    main()
