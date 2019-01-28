import xarray as xr
import datetime as dt
import pandas as pd

import os
import glob

code = {
        'Bering': 3,
        'Beaufort': 13,
        'Chukchi': 12,
        'East_Siberian': 11,
        'Laptev': 10,
        'Kara': 9,
        'Barents': 8, 
        'Central_Arctic': 15,
        'Greenland': 7,
        'CAA': 14,
        'Labrador': 6,
        'Hudson Bay': 4,
        'Okhutsk': 2,
        }

def convert_time(time):
    """Deals with non-CF compliant times"""
    base_time = dt.datetime(1900,1,1,0,0,0)
    return [base_time + dt.timedelta(days=d) for d in time]

def get_data():
    """
    Gets PIOMAS thickness data in Nh50km grid
    """
    dirpath = '/disks/arctic5_raid/abarrett/PIOMAS/v2.1/heff'
    filelist = glob.glob(os.path.join(dirpath,'heff.Nh50km.H????.nc4'))
    ds = xr.open_mfdataset(filelist, concat_dim='time', decode_times=False, data_vars='minimal')
    ds['time'] = pd.date_range('1979-01-01', periods=ds.time.shape[0], freq='MS')
    
    return ds

def get_mask():
    """
    Gets Nh50km Arctic region mask
    """
    filepath = '/oldhome/apbarret/projects/ancillary/masks/'+\
               'Arctic_region_mask_Meier_AnnGlaciol2007_Nh50km.tif'
    return xr.open_rasterio(filepath).squeeze()

def main():

    ds = get_data()

    mask = get_mask()

    df = pd.DataFrame({key: ds['heff'].where((mask == value) & (ds['heff'] > 0.)).mean(dim=['x','y']).to_series() for  key, value in code.items()})

    df.to_csv('arctic_region_mean_piomas_heff_nozero.csv')
    
if __name__ == "__main__":
    main()

    
