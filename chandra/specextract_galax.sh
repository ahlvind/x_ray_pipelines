#! /bin/bash
# This script creates grouped spectra, source & background spectra, rmf and arf files from eventfiles for galaxies as source, for respective observation

input="/Users/juliaahlvind/Documents/projekt_1/Chandra/main_list_GALAX.txt"


cd ../../Chandra/data

{
    read -r
    while read name obsid typ dist epoch expt nH RA DEC counts nH eh
    do
        cd $obsid
        cd "repro"
        
        src="src_${name}.reg"
        galax="galax_${name}.reg"
        bkg="bkg_${name}.reg"

        # creates a source spectrum for the galaxy, rmf, arf, bkg and grouped spectrum
        specextract "event_file_0312_gti.fits[sky=region(${galax})]" bkgfile="event_file_0312_gti.fits[sky=region(${bkg})]" grouptype=NUM_CTS binspec=1 outroot=/Users/juliaahlvind/Documents/projekt_1/Chandra/data/$obsid/repro/spectra_galax_${name}.fits clobber=yes weight=no weight_rmf=no correctpsf=no
        
        mv spectra_galax_${name}.fits_grp.pi galax_${name}_acis_grp.pi
        
        
        cd ../..
        
    done
}< "$input"

cd /Users/juliaahlvind/Documents/projekt_1/pipelines/chandra
