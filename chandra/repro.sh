#! /bin/bash
# This script creates eventfiles in the energy range 0.3 to 12keV from raw Chandra data. It also calles a python function, lightcurve_chandra.py, which plots the light curve and prints the effective time that is left after deflaring the observation. The final product, event_file_0312_gti.fits, is what is used to create spectra.

input="/Users/juliaahlvind/Documents/projekt_1/Chandra/main_list_download.txt"

cd ../../Chandra/data
while read Name obsid Type epoch dist expt dateSN dateXray ra dec
do
    # find filename, depending on the obsid there are differnet nr of 0 before the id
    FILE="$obsid"
    if [ -d "$FILE" ];
    then
        number="$obsid"
        if [ ${#number} -eq 3 ]
        then
            full_obsid="00$obsid"
            
        elif [ ${#number} -eq 4 ]
        then
            full_obsid="0$obsid"
            
        else
            full_obsid=$obsid
        fi
        
        echo "folder $obsid exists"
        echo $full_obsid
        
        chandra_repro $obsid"/" outdir=$obsid"/repro" clobber=yes
        
        cd $obsid
        cd "repro"
        
        dmcopy infile="acisf${full_obsid}_repro_evt2.fits[energy=300:12000]" outfile=event_file_0312.fits  clobber=yes

        dmextract infile="event_file_0312.fits[bin time=::200]" outfile='lightcurve_0312.fits' opt=ltc1  clobber=yes
        
        dmcopy "acisf${full_obsid}_repro_evt2.fits[bin X=1,Y=1][energy=300:12000]" image_0312_bin1.fits clobber=yes
        
        cd ../..
        pwd
        
        # plot light curves
        python lightcurve_chandra.py $obsid
        
        cd $obsid
        cd "repro"
        
        # deflare the observation/remove flaring time
        deflare lightcurve_0312.fits lightcurve_0312.gti method=clean save=deflare.png
        
        dmcopy "event_file_0312.fits[@lightcurve_0312.gti]" event_file_0312_gti.fits clobber=yes
        
        punlearn ardlib
        
        cd ../..
    fi
done < "$input"

cd /Users/juliaahlvind/Documents/projekt_1/pipelines/chandra
