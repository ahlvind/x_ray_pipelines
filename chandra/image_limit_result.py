
import numpy as np
from astropy.io import ascii
import os
from scipy.integrate import quad
import matplotlib.pyplot as plt


# Set up and files to read
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should giv a string like '/Users/juliaahlvind'
wdir = f'{home}/Documents/projekt_1/Chandra/data'  # path tp mother directory that includes all obsid folders
folder='22feb'
Mpc_cm = 3.08567758*1e24  # mult with distance
galax = '' #'_GALAX'

if galax=='_GALAX':
    file = ascii.read(f'{home}/Documents/projekt_1/Chandra/main_list_Galax.txt', delimiter='\t')
else:
    file = ascii.read(f'{home}/Documents/projekt_1/Chandra/main_list_EjGalax.txt', delimiter='\t')

def days_to_yrs(x):
    return np.divide(x, 365.25)
def yrs_to_days(x):
    return np.multiply(x, 365.25)




def main_func(file, E_int,galax):
    name = file['Name']
    obsid = file['obsid']
    dist = file['dist(Mpc)']
    type = file['type']
    epoch = file['epoch']
    expt = file['efft(ks)']
    nH = file['nH']

    if E_int=='05_2':
        l = 0.5 + (2-0.5)/2
    elif E_int=='2_7':
        l = 2.0 + (7-2)/2
    elif E_int=='05_10':
        l = 0.5 + (10-0.5)/2
    elif E_int=='2_10':
        l = 2.0 + (10-2)/2
    else:
        l = 7 + (10-7)/2

    f = open(f'{home}/Documents/projekt_1/Chandra/{folder}/image_limit_output{E_int}{galax}.csv','+w')
    f.write('Name \t obsid \t dist(Mpc) \t epoch \t type \t efft(ks) \t dateSN \t dateXray \t nH \t modFlux(3σ) \t lower_modFlux(3σ) \t upper_modFlux(3σ) \t Lum(3σ) \t lower_Lum(3σ) \t upper_Lum(3σ) \t countrate \t lower_countrate \t upper_countrate \t upper_countrate_keV \n ')

    for i in range(len(obsid)):

        f.write(name[i]+'\t')
        f.write(str(obsid[i])+'\t')
        f.write(str(dist[i])+'\t')
        f.write(str(epoch[i])+'\t')
        f.write(str(type[i])+'\t')
        f.write(str(expt[i])+'\t')
        f.write(str(file['dateSN'][i])+'\t')
        f.write(str(file['dateXray'][i])+'\t')
        f.write(str(nH[i])+'\t')
               
        
        # go into each repro file of the corresponding obsid
        list_files_in_repro = os.listdir(wdir+'/'+str(obsid[i])+"/repro")
        spec_dir = f'{wdir}/'+str(obsid[i])+'/repro/'
        txt_fil_05 = spec_dir+'image_limit_'+E_int+'_'+str(name[i])+'.txt'
        if 'image_limit_'+E_int+'_'+str(name[i])+'.txt' not in os.listdir(spec_dir):
            print('Finns ej ', name[i], obsid[i])
        else:
            with open(txt_fil_05, 'r') as dat:
                lines = dat.readlines()
            s_rate = lines[6].split()
            s_rate_upp = s_rate[-1].split(',')[1].replace(')','')
            s_rate_low = s_rate[-1].split(',')[0].replace('(','')
            s_rate_best = s_rate[7]
            if 'NAN' in s_rate_upp:
                s_rate_upp = 0
            
            if 'NAN' in s_rate_best:
                s_rate_best = s_rate_upp
            
            if 'NAN' in s_rate_low:
                s_rate_low = 0
            s_rate_upp = float(s_rate_upp)
            s_rate_low = float(s_rate_low)
            s_rate_best = float(s_rate_best)
            
            
            s_modflux = lines[8].split()
            s_obs_mod_flux = lines[9].split()
            midF = s_modflux[1]
            lowerF = s_modflux[-1].split(',')[0].replace('(','')
            upperF = s_modflux[-1].split(',')[1].replace(')','')
            
            if 'NAN' in upperF:
                upperF = 0
                
            if 'NAN' in lowerF:
                lowerF = 0
                
            if 'NAN' in midF:
                midF = upperF
            
            midF = float(midF)
            lowerF = float(lowerF)
            upperF = float(upperF)
            
            midL = midF*4*np.pi*(float(dist[i])*Mpc_cm)**2
            lowerL = lowerF*4*np.pi*(float(dist[i])*Mpc_cm)**2
            upperL = upperF*4*np.pi*(float(dist[i])*Mpc_cm)**2


            f.write(str(midF)+'\t')
            f.write(str(lowerF)+'\t')
            f.write(str(upperF)+'\t')
            f.write(format(midL,'.4e')+'\t')
            f.write(format(lowerL,'.4e')+'\t')
            f.write(format(upperL,'.4e')+'\t')
            f.write(str(s_rate_best)+'\t')
            f.write(str(s_rate_low)+'\t')
            f.write(str(s_rate_upp)+'\t')
            s_rate_upp_kev = s_rate_upp*( l+1/3-0.02/l )
            f.write(str(round(s_rate_upp_kev,2))+'\n')

    f.close()
    return


main_func(file, '05_2', galax)
main_func(file, '2_7', galax)
main_func(file, '05_10', galax)
main_func(file, '2_10', galax)
main_func(file, '7_10', galax)

