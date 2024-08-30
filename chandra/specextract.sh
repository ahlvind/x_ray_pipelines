#! /bin/bash
# This script creates SN grouped spectra, source & background spectra, rmf and arf files from eventfiles for respective observation


input="/Users/juliaahlvind/Documents/projekt_1/Chandra/main_list_EjGALAX.txt"

cd ../../Chandra/data

{
    read -r # skips the first line in txt file, i.e column headers
    while read name obsid Type epoch dist expt efft dateSN dateXray ra dec nH

    do
    
        cd $obsid
        cd "repro"
        
        src="src_${name}.reg"
        bkg="bkg_${name}.reg"
        

        # creates a source spectrum, rmf, arf, bkg and grouped spectrum
        specextract "event_file_0312_gti.fits[sky=region(${src})]" bkgfile="event_file_0312_gti.fits[sky=region(${bkg})]" grouptype=NUM_CTS binspec=1 outroot=/Users/juliaahlvind/Documents/projekt_1/Chandra/data/$obsid/repro/spectra_${name}.fits clobber=yes weight=no weight_rmf=no correctpsf=yes
        
        mv spectra_${name}.fits_grp.pi src_${name}_acis_grp.pi
        
        
        cd ../..
    
    done
    
} < "$input"

cd /Users/juliaahlvind/Documents/projekt_1/pipelines/chandra
