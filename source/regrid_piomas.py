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

varAttrs = {
            'heff': {'long_name': 'sea ice thickness',
                     'units': 'm'},
            'area': {'long_name': 'sea ice concentration',
                     'units': 'none'},
           }
                  
def write_to_netcdf(coords, data, filo, fill_value=1e20, variable='heff', attributes=None):

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
    var = rootgrp.createVariable( variable, 'f4', ('time','x','y',), fill_value=fill_value )

    rootgrp.description = 'PIOMAS sea ice thickness'
    rootgrp.created = dt.datetime.now().strftime('%Y-%m-%d %H:%M')
    rootgrp.created_by = 'A.P.Barrett <apbarret@nsidc.org>'
    rootgrp.source = 'http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid'

    times.long_name = 'time'
    times.units = 'days since 1900-01-01 00:00:00'
    times.calendar = 'gregorian'
    lats.long_name = 'latitude'
    lats.units = 'degrees_north'
    lons.long_name = 'longitude'
    lons.units = 'degrees_east'
    if attributes:
        var.long_name = attributes['long_name']
        var.units = attributes['units']

    times[:] = date2num(t, units=times.units, calendar=times.calendar)
    lons[:,:] = coords[1]
    lats[:,:] = coords[1]
    var[:,:,:] = np.where(np.isnan(data), fill_value, data)

    rootgrp.close()
    
    return

def make_plot(srcCoord, srcData, dstCoord, dstData, month=1):
    
    it = month-1
    
    fig = plt.figure( figsize=(20,8) )
    gs = gridspec.GridSpec(2, 2, height_ratios=(30,1))
    ax1 = fig.add_subplot(gs[0], projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180,180,40,90], ccrs.PlateCarree())
    
    pcm = ax1.pcolormesh(srcCoord[0], srcCoord[1], srcData[it,:,:], cmap = 'gist_ncar',
                         transform=ccrs.PlateCarree(), vmin = 0, vmax = 6)
    ax1.add_feature(cfeature.COASTLINE)
    ax1.set_title('PIOMAS Ice Thickness - Original')

    ax2 = fig.add_subplot(gs[1], projection=ccrs.NorthPolarStereo())
    ax2.set_extent([-180,180,40,90], ccrs.PlateCarree())
    ax2.pcolormesh(dstCoord[0], dstCoord[1], dstData[it,:,:], cmap = 'gist_ncar',
                   transform=ccrs.PlateCarree(), vmin = 0, vmax = 6)
    ax2.add_feature(cfeature.COASTLINE)
    ax2.set_title('Regridded')

    plt.show()
    

def regrid_piomas(fileList, variable='heff', outgrid='Na12', verbose=False, doplot=False, nowrite=False):
    # Set up regridding
    # -----------------

    # Get PIOMAS coordinates and mask
    if verbose: print ('Getting PIOMAS grid definition')
    srcCoord = get_coords()
    srcMask = get_mask()

    # Get output mask if one exists
    dstMaskFile = 'piomas_mask_{:s}.npy'.format(outgrid)
    if os.path.isfile(dstMaskFile):
        dstMask = np.load(dstMaskFile)
    else:
        print ('No mask exists for grid {:s}: Destination grid will not be masked\n'+
               '      Use regrid_piomas_mask to create a mask')

    # Get EASE grid coordinates
    if verbose: print ('Getting EASE grid definition')
    dstCoord = get_ease_coords(outgrid)

    if verbose: print ('Calculating regridding weights')
    regrid, srcField, dstField = ESMF_GenRegridWeights(srcCoord, dstCoord,
                                                       srcMask=srcMask)

    # Work though fileList
    for year in yearList:
        
        # Get PIOMAS ice thickness fields
        diri = '/disks/arctic5_raid/abarrett/PIOMAS/v2.1/{:s}'.format(variable)
        fili = os.path.join(diri, '{:s}.H{:4d}'.format(variable,year))
        if not os.path.isfile(fili):
            fili = fili+'.gz'
            if not os.path.isfile(fili):
                print ( 'regrid_piomas: {:s} or {:s} dp not exist.  Skipping regrid'.format(fili, fili.replace('.gz','')) )
                continue
        if verbose: print ('Opening {:s}...'.format(fili))
        srcData = read_grid(fili)

        # Regrid fields
        dstData = np.full([srcData.shape[0], dstCoord[0].shape[0], dstCoord[0].shape[1]], np.nan)
        if verbose: print ('Regridding...')
        for it in np.arange(0,srcData.shape[0]):
            srcField.data[...] = np.where(srcMask.T == 1, srcData[it,:,:].T, np.nan)
            dstField = regrid(srcField, dstField)
            if 'dstMask' in dir():
                dstData[it,:,:] = np.where(dstMask == 0, np.nan, dstField.data.T)
            else:
                dstData[it,:,:] = dstField.data.T
            

        # Write to netCDF
        if not nowrite:
            filo = fili.replace('{:s}.'.format(variable), '{:s}.{:s}.'.format(variable, outgrid)).replace('.gz','')+'.nc4'
            if verbose: print ('Writing regridded data to '+filo)
            write_to_netcdf(dstCoord, dstData, filo, variable=variable, attributes=varAttrs[variable])

        # Plot results
        if doplot:
            make_plot(srcCoord, srcData, dstCoord, dstData)
    
if __name__ == "__main__":
    import argparse
    
    # Add arguments
    parser = argparse.ArgumentParser(description='Regrid PIOMAS monthly thickness files')
    parser.add_argument('yearList', metavar='Y', type=str, nargs='+',
                        help='A list of years (Y0 Y1 Y2...) or year range (Y0..Yn) to regrid')
    parser.add_argument('--variable', type=str, default='heff',
                        help='PIOMAS variable name to regrid')
    parser.add_argument('--outgrid', type=str, default='Na25',
                        help='Name of output grid - must match EASE or EASE2 grid name')
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

    regrid_piomas(yearList, variable=args.variable, outgrid=args.outgrid,
                  verbose=args.verbose, doplot=args.doplot, nowrite=args.nowrite)


    
