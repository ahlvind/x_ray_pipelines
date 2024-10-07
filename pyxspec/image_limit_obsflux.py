# This script creates a csv file with the estimated observed SN luminosities for the non detected (or all) SNE using image limit results.

import xspec
from xspec import *
import numpy as np
from astropy.io import ascii,fits
import glob
import os, re, sys
from scipy.integrate import quad
import matplotlib.pyplot as plt
import pandas as pd

# Parameters to adjust befor running
# ––––––––––––––––––––––––––––––––––
folder = '22feb' # folder where output csv is saved
galax = '' # '' or '_GALAX' : with or without galaxy background
telescope = 'XMM_ny' # 'XMM_ny' or 'Chandra'

# Constants
Mpc_cm = 3.08567758*1e24  # mult with distance
kev_erg = 1/624150964.7  # mult with integral output keV


# Set up paths
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should giv a string like '/Users/juliaahlvind'
base = os.path.join(home, 'Documents', 'projekt_1')
wdir = os.path.join(base,telescope,'data')
output_folder = os.path.join(base,telescope,folder)
start_folder = os.getcwd()

# import script of functions
sys.path.insert(0,f'{base}/pipelines')
from pyxspec_functions import *

fileobsflux = read_csv(f'{output_folder}/resultat/detected_obsflux{galax}.csv')

# setting up cstatistics and limited iterations
Fit.nIterations = 100
Fit.statMethod = "cstat"

E_int = [[0.5,2.0],[2.0,7.0],[7.0,10.0]]
E_int_large = [[0.5,10.0],[2.0,10.0]]
E_int_ignore = ["**-0.5 2.0-**", "**-2.0 7.0-**", "**-7.0 10.0-**"]
ener_interval = ["0.5 2.0", "2.0 7.0", "7.0 10.0"]
ener_interval_large = ["0.5 10.0", "2.0 10.0"]

# Incorrectly classified SNe or IIn
err_class = ['1996cr', '1994W', '1995N', '1997ab', '1997eg', '2008fq', '1986J', '1978K','2014C', '2004dk', '2019yvr', '1954J','1997dc', '1998bw', '1993J'] # 1997dc är en AGN, 1954J false SN?


def master_func(name, obsid, i, up_countrate05_10_i, up_countrate2_10_i, folder, galax):
    # This function writes an output, a line in a csv file for each observation with results of "fitting" a galactic absorbed power law to a fake spectra. In reality the model is simply scaled such that predicted count rate matches the image limits. Since these are the observed SN luminosities/fluxes, the image limits are based on the full energy range of interest, that is 0.5-10.0 and 2.0-10.0 keV.
    
    if name not in err_class:
        
        if telescope=='XMM_ny':
            repro_pn = 'pn'
            src = f'{name}_fake.pha'
            arf = f'{name}_pn.arf'
            rmf = f'{name}_pn.rmf'
            bkg = f'{name}_bkg_pn.pha'

        else:
            repro_pn = 'repro'
            bkg = f'spectra_{name}.fits_bkg.pi'
            src = f'{name}_fake.pi'
            rmf = f'spectra_{name}.fits.rmf'
            arf = f'spectra_{name}.fits.corr.arf'
            
        # go into each repro file of the corresponding obsid
        list_files_in_repro = os.listdir(f'{wdir}/{obsid}/{repro_pn}')
        spec_dir = f'{wdir}/{obsid}/{repro_pn}'

        # find the spectral file
        current = os.getcwd()  # current folder we are in
        os.chdir(spec_dir) # same as cd but for "python environment"

        # if spectral file does exist
        if (bkg in list_files_in_repro):

            # write out on the list
            f.write(name+'\t')
            if telescope == 'XMM_ny':
                f.write(f's{obsid}\t')
            else:
                f.write(f'{obsid}\t')
            f.write(str(dist[i])+'\t')
            f.write(str(epoch[i])+'\t')
            f.write(str(type[i])+'\t')
                    
            hdul = fits.open(bkg)
            expt1 = hdul[1].header['EXPOSURE']
            f.write(str(efft[i])+'\t')
            f.write(str(expt1/1e3)+'\t')
            f.write(str(obsid_file['dateSN'][i])+'\t')
            f.write(str(obsid_file['dateXray'][i])+'\t')
            f.write(str(nH[i])+'\t')
            
            for m in range(2):  # Loop over each energy interval for cstat

                # apply model
                mo = Model("tbabs*po")
                
                # define model parameters
                # –––––––––––––––––––––––
                # tbabs
                mo(1).values = nH[i]
                mo(1).frozen = True
                
                #  powerlaw
                mo(2).values = 2
                mo(2).frozen = True
                mo(3).values = 1
                        
                if src in os.listdir(spec_dir): # first delete old fakit spectra
                    os.remove(src)
                
                hdul = fits.open(bkg)
                expt = hdul[1].header['EXPOSURE']
                fs1 = FakeitSettings(response=rmf, arf=arf, exposure=expt, fileName=src) #, background=name+"_bkg_pn.pha"
                AllData.fakeit(1,fs1)
                
                s = Spectrum(src)
                
                AllData.ignore("**-0.5 10.0-**")
            
                
                # extract values from fit
                # –––––––––––––––––––––––
                # photon index
                gamma = mo(2).values[0]

                # normalinsation of powerlaw
                norm = mo(3).values[0]

                # if we want upper limit or adjusted to limits
                if m==0:
                    image_limit = up_countrate05_10_i
                elif m==1:
                    image_limit = up_countrate2_10_i

            
                # change the norm of PL to reach countrate given by the image limit
                predicted_rate = float(s.rate[-1])
                if image_limit>-10:
                    mo(3).values = mo(3).values[0]*image_limit/predicted_rate
                else:
                    mo(3).values = 0
                mo(3).frozen = True
                norm = mo(3).values[0]
                new_predicted_rate = s.rate[-1]

                # write out the results
                f.write(str(image_limit)+'\t')
                f.write(format(new_predicted_rate,'.4e')+'\t')
                f.write(format(norm,'.4e')+'\t')
                
                # obsflux flux
                AllModels.calcFlux(ener_interval_large[m])
                flux = s.flux[0]
                f.write(format(flux,'.4e')+'\t')
                # calc obs. luminosity
                L, L_min, L_max = lum_funk(flux,flux,flux,dist[i])
                if m==1:
                    f.write(format(L,'.4e')+'\n')
                else:
                    f.write(format(L,'.4e')+'\t')
                    
            mo = None
            xspec.AllData.clear()
            s = None
            xspec.AllModels.clear()

    return

## ===================================
obsid_file = read_csv(f'{output_folder}/image_limit_output05_2{galax}.csv')
dist = obsid_file['dist(Mpc)']
type = obsid_file['type']
epoch = obsid_file['epoch']
efft = obsid_file['efft(ks)']
nH = obsid_file['nH']

file05_10 = read_csv(f'{output_folder}/image_limit_output05_10{galax}.csv')
file2_10 = read_csv(f'{output_folder}/image_limit_output2_10{galax}.csv')
up_countrate05_10 = file05_10['count_rate_up_lim']
countrate05_10 = file05_10['count_rate']
up_countrate2_10 = file2_10['count_rate_up_lim']
countrate2_10 = file2_10['count_rate']


f = open(f'{output_folder}/resultat/IMAGE_LIMITS_obsflux{galax}.csv','+w')
f.write('"Name"\t"obsid"\t"dist(Mpc)"\t"epoch"\t"type"\t"efft(ks)"\t"expt_bkg_file(ks)"\t"dateSN"\t"dateXray"\t"mean_nH(1e22atoms_cm-2)"\t"image_limit_counts(05_2)"\t"image_limit_counts(05_10)"\t"model_predicted_counts(05_10)"\t"norm_PL(05_10)"\t"obsFlux_05_10"\t"obslum_05_10"\t"image_limit_counts(2_10)"\t"model_predicted_counts(2_10)"\t"norm_PL(2_10)"\t"obsFlux_2_10"\t"obslum_2_10)"\n')

# Start of the actuall program. Loop over each obsid
# ––––––––––––––––––––––––––––––––––––––––––––––––––
for i in range(len(obsid_file)):
    
        # remove potential space in name
        name = obsid_file['Name'][i].replace(' ', '')
        
        if telescope=='XMM_ny':
            # previous csv file ar emissing 0 in front of XMM obsid, correcting for this here
            if obsid_file['obsid'][i]==51610101:
                obsid = f'00{obsid_s[i]}'
            else:
                obsid = f'0{obsid_s[i]}'
            
        else:
            obsid = obsid_file['obsid'][i]
            
        up_countrate05_10_i = up_countrate05_10[i]
        up_countrate2_10_i = up_countrate2_10[i]
        
        master_func(name, obsid, i, up_countrate05_10_i, up_countrate2_10_i, folder, galax)
            
f.close()


os.chdir(start_folder)
