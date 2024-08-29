#! /bin/bash
# This script is using srcflux to extract count rate for the SN position. The output are txt files with countrates and upper limits. These txt files are created and stored in the respective obsid folder. The count rates and txt files are created for 5 different energy ranges: 0.5-2, 2-7, 7-10, 0.5-10 and 2-10keV. This script uses different backgrounds since the SN with prominent galaxy background need other annulus bkg regions.

#input="/Users/juliaahlvind/Documents/projekt_1/Chandra/main_list_GALAX.txt"

punlearn srcflux

cd ../../Chandra/data
{
    read -r
    while read name obsid Type epoch dist expt efft dateSN dateXray RA DEC nH
    do
        
        cd $obsid
        cd "repro"
        

        RA_s=$RA
        DEC_s=$DEC
        
        if [ -f spectra_${name}.fits.arf ]
        then
            arf=spectra_${name}.fits.corr.arf
        else
            echo "${name} ${obsid} har ej arf"
        fi

        curent_dir="/Users/juliaahlvind/Documents/projekt_1/Chandra/data/${obsid}/repro"
        prefix='./'
        
        fovfile_path="$(find . -name '*repro_fov1.fits*')"
        fov_file=${fovfile_path/#$prefix}
        asolfile_path="$(find . -name '*asol1.fits*')"
        asol_file=${asolfile_path/#$prefix}
        mskfile_path="$(find . -name '*_msk1.fits*')"
        msk_file=${mskfile_path/#$prefix}
        bpixfile_path="$(find . -name '*repro_bpix1.fits*')"
        bpix_file=${bpixfile_path/#$prefix}
        
        energ_lo=(0.5 2.0 7.0 0.5 2.0)
        energ_hi=(2.0 7.0 10.0 10.0 10.0)
        # Get the length of the arrays (assuming both are the same length)
        length=${#energ_lo[@]}
        
        for ((i=0; i<length; i++))
        do
        
            dmtcalc $arf arf_weights expression="mid_energy=(${energ_lo[$i]}+${energ_hi[$i]})/2.0;weights=(mid_energy*specresp)" clob+
            dmstat "arf_weights[mid_energy=${energ_lo[$i]}:${energ_hi[$i]}][cols weights,specresp]" verbose=0
            res=$(pget dmstat out_sum)
            a=${res%,*}
            b=${res#*,}
            fin=$(python -c "print($a/$b)")
            echo "SÖK"
            echo $fin
            srcflux event_file_0312_gti.fits "${RA_s},${DEC_s}" outroot=rhooph${energ_lo[$i]}_${energ_hi[$i]}_${name} srcreg="region(src_${name}.reg)" bkgreg="region(galax_an_${name}.reg)" conf=0.997 bands=${energ_lo[$i]}:${energ_hi[$i]}:${fin} binsize=1 rmffile="spectra_${name}.fits.rmf" arffile="$arf" fovfile=$fov_file asolfile=$asol_file mskfile=$msk_file bpixfile=$bpix_file model=xsphabs.abs1*xspowerlaw.pow1 paramvals=abs1.nH=${nH} clobber=yes  ;
            
            mv rhooph${energ_lo[$i]}_${energ_hi[$i]}_${name}_summary.txt image_limit_${energ_lo[$i]}_${energ_hi[$i]}_${name}.txt
        
        done
        rm *rhooph*
        
        
        cd ../..
    done
}< "$input"


