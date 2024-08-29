#! /bin/bash
# This script runs ciao wavdetect, that finds sources in the Chandra images


# Start the timer
start_time=$(date +%s)

input="/Users/juliaahlvind/Documents/projekt_1/Chandra/main_list_Galax.txt"

cd ../../Chandra/data
{
    read -r
    while read name obsid epoch Type dist expt efft dateSN dateXray ra dec nH
    do
        
        cd $obsid
        cd "repro"
        
        # create psf map
        mkpsfmap infile="image03_12_${name}.fits" outfile="mypsfmap_${name}.fits" energy=1.49 ecf=0.90 clobber=yes

        
        fits_file="image03_12_${name}.fits"
        livetime_full=$(fitsheader "$fits_file" | grep "LIVETIME" | awk -F '=' '{print $2}' | tr -d ' ')
        livetime=$(echo $livetime_full | sed 's/[^0-9.E+-]//g')
        echo $livetime

        # run wavdetect for given interval using scalesizes 1,2 3 and a signal threshold 1e-5
        set clobber=yes
        wavdetect infile="image03_12_${name}.fits" outfile="outfile_${name}.fits" imagefile="image_wavdetect_${name}.fits" defnbkgfile="background03_12_${name}.fits" scellfile="new_sources_${name}.fits" scales="1.0 2.0 3.0" exptime=$livetime psffile="mypsfmap_${name}.fits" regfile="src_wavdetect_${name}.reg" sigthresh=1e-5 clobber=yes
            
        fi
        
        cd ../..
        
    done
}< "$input"

# Stop the timer
end_time=$(date +%s)

# Calculate and display the elapsed time
elapsed_time=$(((end_time - start_time)/60))
echo "Script execution time: $elapsed_time minutes"

cd /Users/juliaahlvind/Documents/projekt_1/pipelines/chandra
