#---------------------------------------------------------------------
# Tools to regrid PIOMAS data to EASE grids
#
# 2018-03-14 A.P.Barrett <apbarret@nsidc.org>
#---------------------------------------------------------------------

import numpy as np
import os
import ESMF

#import geopy.distance

gdir = '/oldhome/apbarret/projects/ancillary/maps'

meta = {
        'Na12': {'xdim': 722,
                 'ydim': 722,
                 'latfile': os.path.join(gdir,'ycenter.Na12.722x722x1.float'),
                 'lonfile': os.path.join(gdir,'xcenter.Na12.722x722x1.float')},
        'Nh':   {'xdim': 1441,
                 'ydim': 1441,
                 'latfile': os.path.join(gdir,'ycenter.Nh.1441x1441x1.float'),
                 'lonfile': os.path.join(gdir,'xcenter.Nh.1441x1441x1.float')},
        'AWI_25km': {'xdim': 129,
                      'ydim': 104,
                      'latfile': os.path.join(gdir,'ycenter.AWI_25km.104x129x1.float'),
                      'lonfile': os.path.join(gdir,'xcenter.AWI_25km.104x129x1.float')},
        'Nh50km': {'xdim': 360,
                   'ydim': 360,
                   'latfile': os.path.join(gdir,'ycenter.nh50km_nest.360x360x1.float'),
                   'lonfile': os.path.join(gdir,'xcenter.nh50km_nest.360x360x1.float')},
       }

def get_ease_coords(gridname):
    """
    Gets lon, lat coordinates for an EASE grid

    Argument
    --------
    gridname - name of EASE grid, e.g. Na12
    
    Returns
    -------
    Tuple containing lon, lat
    """

    latfile = meta[gridname]['latfile']
    lonfile = meta[gridname]['lonfile']
    xdim = meta[gridname]['xdim']
    ydim = meta[gridname]['ydim'] 

    lat = np.fromfile(latfile, dtype='f4').reshape(xdim,ydim)
    lon = np.fromfile(lonfile, dtype='f4').reshape(xdim,ydim)
 
    return (lon, lat)

def find_nearest(ingrid, outgrid):
    """
    Finds the index of nearest ingrid cells to outgrid cells

    Argument
    --------
    ingrid - tuple containing (lon, lat) for source grid
    outgrid - tuple containing (lon, lat) for destimation grid

    Returns
    -------
    1D numpy array containing indices
    """
 
    ilon = ingrid[0]
    ilat = ingrid[1]

    olon = outgrid[0]
    olat = outgrid[1]

    idx = []

    for ox, oy in zip(olon.flat[0:10], olat.flat[0:10]):
        if (ox < -9998.) | (oy < -9998,): continue
        d = np.array([geopy.distance.great_circle((ox,oy),(ix,iy)).km \
                      for ix, iy in zip(ilon.flat, ilat.flat)])
        print (ox, oy, ix, iy, d.min(), d.argmin())
        idx.append(d.argmin())

    return np.array(idx)

def make_regrid_index(easegridnm, SAVEFILE=True):
    """
    Creates a file containing indices based on nearest neighbour method.  Indices are
    written to a numpy file.
 
    Arguments
    ---------
    easegridnm - name of ease grid, e.g. Na12

    """

    from readers import get_coords

    ingrid = get_coords()
    outgrid = get_ease_coords(easegridnm)

    indices = find_nearest(ingrid, outgrid)

    if SAVEFILE:
        filo = 'piomas2{:}_indices.npy'
        np.save(filo, indices)

    return indices

def griddef(gridname):
    """
    Returns a tuple of grid parameters for an EASE grid

    gridname - standard name of grid
    """
    gdef = {
            'Nh': {'c': 12.5, 'nx': 1441, 'ny': 1441, 'r0': 720, 's0': 720}
            }

    try:
        result = gdef[gridname]
        return result
    except:
        print ('griddef: unknown grid name')
#        return -1

def ll2northEase(ingrid, gridname, radius=6371.228):
    """
    Calculates EASE grid column and row coordinates for a given grid definition

    Arguments
    ---------
    ingrid - (lon,lat) tuple for input grid
    gridname - standard ease grid name

    Returns
    -------
    EASE grid row column coordinates
    """

    prm = griddef(gridname)

    rlon = np.radians(ingrid[0])
    rlat = np.radians(ingrid[1])
    
    r = 2*radius/prm['c']*np.sin(rlon)*np.sin((np.pi/4.)-(rlat/2.)) + prm['r0']
    s = 2*radius/prm['c']*np.cos(rlon)*np.sin((np.pi/4.)-(rlat/2.)) + prm['s0']

    return (r, s)

def ESMF_GenRegridWeights(srcCoord, dstCoord, srcMask=None, dstMask=None,
                          src_mask_values=None, dst_mask_values=None,
                          method='bilinear'):
    """
    Generates an ESMF regrid object containing regrid weights

    Arguments
    ---------
    srcGrid - tuple containing source grid lons and lats
    dstGrid - tuple containing destination grid lons and lats
    srcMask - ndarray containing a mask for the source grid
    dstMask - ndarray containing a mask for the source grid
    src_mask_values - ndarray containing values to mask out in srcMask
    dst_mask_values - ndarray containing values to mask out in dstMask
    method - Interpolation method ('bilinear', 'nearest', 'conserve', 'patch'
             See <http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/RegridMethod.html#ESMF.api.constants.RegridMethod>

    N.B. masking does not appear to work currently

    Returns
    -------
    An ESMF regrid object
    """

    methods = {'bilinear': ESMF.RegridMethod.BILINEAR,
              'nearest_stod': ESMF.RegridMethod.NEAREST_STOD,
              'nearest_dtos': ESMF.RegridMethod.NEAREST_DTOS,
              'patch': ESMF.RegridMethod.PATCH,
              'conserve': ESMF.RegridMethod.CONSERVE}
              
    # Define grids
    # - grids are transposed to make ESMF efficient.  ESMF is Fortran based, so transposing grids
    #   makes them Fortran contiguous
    srcGrid_shape = srcCoord[0].T.shape
    dstGrid_shape = dstCoord[0].T.shape
    
    sourcegrid = ESMF.Grid(np.array(srcGrid_shape), staggerloc=ESMF.StaggerLoc.CENTER, coord_sys=ESMF.CoordSys.SPH_DEG)
    destgrid = ESMF.Grid(np.array(dstGrid_shape), staggerloc=ESMF.StaggerLoc.CENTER, coord_sys=ESMF.CoordSys.SPH_DEG)

    # Add mask to destination grid
    if isinstance(srcMask, np.ndarray):
        smask = sourcegrid.add_item(ESMF.GridItem.MASK, staggerloc=ESMF.StaggerLoc.CENTER)
        smask[...] = srcMask.T

    if isinstance(dstMask, np.ndarray):
        print ('HERE')
        dmask = destgrid.add_item(ESMF.GridItem.MASK, staggerloc=ESMF.StaggerLoc.CENTER)
        dmask[...] = dstMask.T
        
    # Assign grid coordinates
    source_lon = sourcegrid.get_coords(0)
    source_lat = sourcegrid.get_coords(1)
    source_lon[...] = srcCoord[0].T
    source_lat[...] = srcCoord[1].T

    dest_lon = destgrid.get_coords(0)
    dest_lat = destgrid.get_coords(1)
    dest_lon[...] = dstCoord[0].T
    dest_lat[...] = dstCoord[1].T

    # Define fields
    sourcefield = ESMF.Field(sourcegrid, name='PIOMAS Ice Thickness')
    destfield = ESMF.Field(destgrid, name='Regridded Ice Thickness')

    regrid = ESMF.Regrid(sourcefield, destfield, 
                         regrid_method=methods[method],
                         unmapped_action=ESMF.UnmappedAction.IGNORE,
                         src_mask_values=src_mask_values,
                         dst_mask_values=dst_mask_values )

    return regrid, sourcefield, destfield

def main():

    indices = make_regrid_index('Na12', SAVEFILE=True)

    return

if __name__ == "__main__":
    main()

