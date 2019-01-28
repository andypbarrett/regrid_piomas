#----------------------------------------------------------------------
# A module containing readers for PIOMAS data
#
# 2018-03-13 A.P.Barrett <apbarret@nsidc.org>
#----------------------------------------------------------------------
import numpy as np
import os
import re
import gzip

def read_grid(fili):
    """
    Reads PIOMAS data files.  time dimension is inferred from
    size of array.

    Argument
    --------
    fili - file path

    Returns
    -------
    numpy array with dimensions (t,x,y)
    """
    
    nrow = 360
    ncol = 120
    
    if re.search('gz$', fili):
        f = gzip.GzipFile(fili)
        s = f.read()
        f.close()
        data = np.frombuffer(s, dtype='f4')
    else:
        data = np.fromfile(fili, dtype='f4')

    len = data.shape[0]
    ntim = float(len) / (float(nrow)*float(ncol))

    return data.reshape(int(ntim),ncol,nrow)

def get_coords(fili=None):
    """
    Gets a coordinate array for the PIOMAS grid

    Keyword
    -------
    fili - filepath. Defaults to 
       /disks/arctic5_raid/abarrett/PIOMAS/utilities/grid.dat
    
    Returns
    -------
    lon, lat arrays
    """

    ncol = 120
    nrow = 360

    if not fili:
        diri = '/disks/arctic5_raid/abarrett/PIOMAS/utilities'
        fili = os.path.join(diri,'grid.dat')
        
    tmp = np.loadtxt(fili, dtype='float')
    lon = tmp.reshape(-1)[:ncol*nrow].reshape(ncol,nrow)
    lat = tmp.reshape(-1)[ncol*nrow:].reshape(ncol,nrow)

    return lon, lat

def get_mask(fili=None):
    """
    Gets the PIOMAS mask 0: Land, 1: Ocean.

    The file contains values Land <= 0, Ocean > 0.  I set these to
    0,1 for my purpose.

    Keyword
    -------
    fili - alternative filepath for mask. Default is:
         /disks/arctic5_raid/abarrett/PIOMAS/utilities/io.dat_360_120.output

    Returns
    -------
    Numpy array containing mask values
    """

    chunk_size = 2

    if not fili:
        diri = '/disks/arctic5_raid/abarrett/PIOMAS/utilities'
        fili = os.path.join(diri,'io.dat_360_120.output')

    mask = []
    f = open(fili)
    for s in f.readlines():
        mask.append([np.float(s[i:i+chunk_size]) \
                     for i in range(0,len(s.strip()),chunk_size)])
    return np.where( np.array(mask) > 0, 1, 0)


    
                    
