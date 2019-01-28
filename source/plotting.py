from cartopy import crs, feature
import matplotlib.pyplot as plt
import os

from readers import get_mask, get_coords, read_grid

def main():

    # Test data set
    ddir = '/disks/arctic5_raid/abarrett/PIOMAS/v2.1/heff'
    dfil = 'heff.H2017'
    data = read_grid(os.path.join(ddir,dfil))
    print (data.shape)
    
    lon, lat = get_coords()
    
    # Setup plotting
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='k',
                                        facecolor=feature.COLORS['land'])
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=crs.NorthPolarStereo())

    ax.set_extent([-180,180,65,90], crs.PlateCarree())
    ax.add_feature(land_50m)

    plt.contourf(lon, lat, data[0,:,:], transform=crs.PlateCarree())
    
    plt.show()

if __name__ == "__main__":
    main()


    
