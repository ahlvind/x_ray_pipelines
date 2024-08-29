import subprocess
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
from astropy import units as u

# threshold for how far appart the found source and SN coordinate should be in order to count as an image detection
threshold_sep = 2/3600 # deg



def func(obsid, name):
    # FITS file with counts from wavdetect of found
    fits_file = f'/Users/juliaahlvind/Documents/projekt_1/Chandra/data/{obsid}/repro/new_sources_{name}.fits'
    
    # region file to check if it is indeed one of the detected
    region_file = f'/Users/juliaahlvind/Documents/projekt_1/Chandra/data/{obsid}/repro/src_{name}.reg'

    # Load the FITS file
    hdulist = fits.open(fits_file)
    data = hdulist[0].data

    """# Display the image in DS9
    ds9_cmd = f'. ds9.sh {fits_file} -regions {region_file}'
    subprocess.run(ds9_cmd, shell=True, check=True)

    # Wait for the user to interact with DS9
    input("Press Enter after inspecting the image in DS9...")"""

    # Load the WCS information from the FITS header
    wcs = WCS(hdulist[0].header)
    
    # Read the region file and extract the values from the specified region
    with open(region_file, 'r') as f:
        region_str = f.read().strip()

    # exctract the ra and dec from the region
    RA = region_str.split('(')[1].split(',')[0]
    DEC = region_str.split('(')[1].split(',')[1]

    # Parse the region string to get the region coordinates
    region_coords = SkyCoord(RA,DEC, frame='fk5', unit=(u.hourangle, u.deg))

    # Convert the region coordinates to pixel coordinates
    x, y = wcs.world_to_pixel(region_coords)

    # Extract values from the specified region
    region_values = data[int(y), int(x)]
    
    """print(f'Values inside the specified region: {region_values}')
    
    # close the ds9 window
    subprocess.run('pkill ds9', shell=True, check=True)"""
    
    # get the live time
    expt = hdulist[0].header.get('LIVETIME')
    return region_values, expt

def func2(obsid,name):
    """ This function extracts the source detections significance in terms of sigma"""

    # FITS file with counts from wavdetect of found
    fits_file = f'/Users/juliaahlvind/Documents/projekt_1/Chandra/data/{obsid}/repro/outfile_{name}.fits'

    # region file to check if it is indeed one of the detected
    region_file = f'/Users/juliaahlvind/Documents/projekt_1/Chandra/data/{obsid}/repro/src_{name}.reg'

    # Load the FITS file
    hdulist = fits.open(fits_file)
    data = hdulist[1].data

    # Load the WCS information from the FITS header
    wcs = WCS(hdulist[0].header)

    # Read the region file and extract the values from the specified region
    with open(region_file, 'r') as f:
        region_str = f.read().strip()

    # exctract the ra and dec from the region
    RA = region_str.split('(')[1].split(',')[0]
    DEC = region_str.split('(')[1].split(',')[1]

    # Parse the region string to get the region coordinates
    region_coords = SkyCoord(RA,DEC, frame='fk5', unit=(u.hourangle, u.deg))

    # list all detected sources coordinates
    sources_coords = SkyCoord(data['RA'],data['DEC'], frame='fk5', unit=(u.deg, u.deg))
    
    # if fits file is empty
    if len(data['RA'])==0:
        sig=0
        return sig,0
    
    else:

        # find if there is a match with coordinates of detected sources and src_reg, by a separation limit
        print('fits file coord:', sources_coords)
        chandra_index_best, separation_best, _ = region_coords.match_to_catalog_sky(sources_coords)
        
        if separation_best.deg[0]<threshold_sep:
            # detection
            sig = data['SRC_SIGNIFICANCE'][chandra_index_best]
            print(' ')
            print('Source is detected with:',round(sig,1),'σ')
            return sig, separation_best.deg[0]
        else:
            sig = 0
            print(' ')
            print('The source is not detected :(')
            return sig, 0



# run for main_list_EjGalax.txt or main_list_Galax.txt
with open('main_list_EjGalax.txt', 'r') as dat:
    lines = dat.readlines()
    
f = open('/Users/juliaahlvind/Documents/projekt_1/Chandra/22feb/wavdetect_output_GALAX2.csv', '+w')
f.write('"Name"\t"obsid"\t"type"\t"epoch"\t"dist[Mpc]"\t"expt(ks)"\t"efft(ks)"\t"dateSN"\t"dateXray"\t"ra"\t"dec"\t"nH"\t"count_rate"\t"sigma"\t"sep(arcsec)"\n')
for l in lines[1:]:
    name = l.split()[0]
    obsid = l.split()[1]
    print(name,obsid)
    region_values, expt= func(obsid, name)
    sig, sep = func2(obsid, name)
    print(f'{name} {obsid}: ', round(region_values/expt,6),'counts/s')
    f.write('"'+l.replace('\n','\t').replace('\t','"\t"'))
    f.write(format(region_values/expt,'.4e')+'"\t')
    f.write(str(round(sig,1))+'\t')
    f.write(str(round(sep*3600,1))+'\n')

f.close()

