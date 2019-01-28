import ESMF

import numpy as np
import os, glob
import re

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import cm, gridspec, rcParams
import matplotlib.pyplot as plt

from readers import read_grid, get_coords, get_mask
from regrid import get_ease_coords, ESMF_GenRegridWeights

def write_to_netcdf(coords, data, filo):

    from netCDF4 import Dataset, date2num
    import datetime as dt
    import numpy as np
    import re
    
    nt, nx, ny = data.shape

    # Create time
    yyyy = int( re.search('H(\d{4})\.', filo).groups()[0] )
    t = [dt.datetime(yyyy,mm,1) for mm in range(1,nt,1)]
        
    rootgrp = Dataset(filo, 'w')

    time = rootgrp.createDimension('time', None)
    x = rootgrp.createDimension('x', nx)
    y = rootgrp.createDimension('y', ny)

    times = rootgrp.createVariable( 'time', 'f8', ('time',) )
    lats = rootgrp.createVariable( 'latitude', 'f4', ('x','y',) )
    lons = rootgrp.createVariable( 'longitude', 'f4', ('x','y',) )
    sit = rootgrp.createVariable( 'sit', 'f4', ('time','x','y',), fill_value=1e20 )

    rootgrp.description = 'PIOMAS sea ice thickness'
    rootgrp.created = dt.datetime.now().strftime('%Y-%m-%d %H:%M')
    rootgrp.created_by = 'A.P.Barrett <apbarret@nsidc.org>'
    rootgrp.source = 'http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid'

    times.long_name = 'time'
    times.units = 'days since 1990-01-01 00:00:00'
    times.calendar = 'gregorian'
    lats.long_name = 'latitude'
    lats.units = 'degrees_north'
    lons.long_name = 'longitude'
    lons.units = 'degrees_east'
    sit.long_name = 'sea ice thickness'
    sit.units = 'm'

    times[:] = date2num(t, units=times.units, calendar=times.calendar)
    lons[:,:] = coords[1]
    lats[:,:] = coords[1]
    sit[:,:,:] = data

    rootgrp.close()
    
    return

def make_plot(srcCoord, srcData, dstCoord, dstData, month=1):
    
    it = month-1
    
    fig = plt.figure( figsize=(20,8) )
    gs = gridspec.GridSpec(2, 2, height_ratios=(30,1))
    ax1 = fig.add_subplot(gs[0], projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180,180,40,90], ccrs.PlateCarree())
    
    pcm = ax1.pcolormesh(srcCoord[0], srcCoord[1], srcData, cmap = 'gist_ncar',
                         transform=ccrs.PlateCarree(), vmin = 0, vmax = 6)
    ax1.add_feature(cfeature.COASTLINE)
    ax1.set_title('PIOMAS Mask - Original')

    ax2 = fig.add_subplot(gs[1], projection=ccrs.NorthPolarStereo())
    ax2.set_extent([-180,180,40,90], ccrs.PlateCarree())
    ax2.pcolormesh(dstCoord[0], dstCoord[1], dstData, cmap = 'gist_ncar',
                   transform=ccrs.PlateCarree(), vmin = 0, vmax = 6)
    ax2.add_feature(cfeature.COASTLINE)
    ax2.set_title('Regridded')

    plt.show()
    

def regrid_piomas(fileList, outgrid='Na12', verbose=False, doplot=False, nowrite=False):
    # Set up regridding
    # -----------------

    # Get PIOMAS coordinates and mask
    if verbose: print ('Getting PIOMAS grid definition')
    srcCoord = get_coords()
    srcMask = get_mask()

    dstMaskFile = 'combined_mask_{:s}.npy'.format(outgrid)
    if os.path.isfile(dstMaskFile):
        dstMask = np.load(dstMaskFile)
#    else:
#        print ('No mask exists for grid {:s}: Destination grid will not be masked\n'+
#               '      Use regrid_piomas_mask to create a mask')

    # Get EASE grid coordinates
    if verbose: print ('Getting EASE grid definition')
    dstCoord = get_ease_coords(outgrid)

    if verbose: print ('Calculating regridding weights')
    regrid, srcField, dstField = ESMF_GenRegridWeights(srcCoord, dstCoord,
                                                       srcMask=srcMask, method='nearest_stod')

    # Regrid fields
    if verbose: print ('Regridding...')
    srcField.data[...] = srcMask.T
    dstField = regrid(srcField, dstField)
    dstMask = dstField.data.T
    dstMask = np.where(dstCoord[1] < 50, 0, dstMask)
    
    # Write to netCDF
    if not nowrite:
        filo = 'piomas_mask_{:s}.npy'.format(outgrid)
        if verbose: print ('Writing regridded data to '+filo)
        np.save(filo, dstMask)
        #write_to_netcdf(dstCoord, dstData, filo)

    # Plot results
    if doplot:
        make_plot(srcCoord, srcMask, dstCoord, dstMask)
    
if __name__ == "__main__":
    import argparse
    
    # Add arguments
    parser = argparse.ArgumentParser(description='Regrid PIOMAS monthly thickness files')
    parser.add_argument('yearList', metavar='Y', type=str, nargs='+',
                        help='A list of years (Y0 Y1 Y2...) or year range (Y0..Yn) to regrid')
    parser.add_argument('--outgrid', type=str, default='Na12',
                        help='Name of output grid - must be same as EASE or EASE2 grid')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--doplot', action='store_true', help='Generate a comparison plot')
    parser.add_argument('--nowrite', action='store_true', help='Do not write to file')
    args = parser.parse_args()

    # generate year list
    m = re.search('(\d{4})\.\.(\d{4})', args.yearList[0])
    if m:
        yearList = range( int(m.groups()[0]), int(m.groups()[1])+1 )
    else:
        yearList = [int(y) for y in args.yearList]

    regrid_piomas(yearList, outgrid=args.outgrid, verbose=args.verbose, doplot=args.doplot, nowrite=args.nowrite)


    
