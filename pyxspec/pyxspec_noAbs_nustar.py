import sys
import os
import numpy.ma as ma
from nustar_plot_tool import nustar_3_window_plot

# Set up and files to read
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should giv a string like '/Users/juliaahlvind'
wdir = f'{home}/Documents/projekt_1/NuSTAR/data'  # path tp mother directory that includes all obsid folders
base = f'{home}/Documents/projekt_1/NuSTAR'  # path tp mother directory that includes all obsid folders
projekt_1 = f'{home}/Documents/projekt_1'
folder = '22feb'

sys.path.insert(0,f'{projekt_1}/pipelines')
from pyxspec_functions import *

sys.path.insert(0,f'{projekt_1}/pipelines/plotting')
from nustar_plot_tool import nustar_3_window_plot
                
                
# setting up cstatistics and limited iterations
Fit.nIterations = 100
Fit.statMethod = "cstat"
Fit.query = "yes"

# Chandra file with SN names and obsid
obsid_file = read_csv(f'{projekt_1}/{folder}/ObsFlux_test.csv')
obsid = obsid_file['obsid']
name_file = obsid_file['Name']
dist = obsid_file['dist(Mpc)']
type = obsid_file['type']
epoch = obsid_file['epoch']
date = obsid_file['dateXray']

nH_file = ascii.read(f"{projekt_1}/nH_60Mpc.csv", delimiter=';')
nH_all = nH_file['NH_tot_(Mean)_[atoms_cm-2]']
nH_names = nH_file['Name']
string_nH=[]
for k in nH_names:
    string_nH.append(k)
    
Mpc_cm = 3.08567758*1e24  # mult with distance
kev_erg = 1/624150964.7  # mult with integral output keV


E_int = [[2.0,10.0],[7.0,10.0],[2.0,78.0],[10.0,30.0], [10.0,78.0]]
    
f = open(f'{projekt_1}/{folder}/noAbs_test2.csv','+w')  # creates a csv document with all outputs
f.write('"Name"\t"obsid"\t"dist(Mpc)"\t"type"\t"epoch"\t"exp_timeA(ks)"\t"exp_timeB(ks)"\t"mean_nH(1e22atoms_cm-2)"\t"RA"\t"DEC"\t"dateXray"\t"countsA"\t"countsB"\t"PhoIndex"\t"confInt(1.6σ)"\t"norm"\t"confInt(1.6σ)"\t"kT"\t"confInt(1.6σ)"\t"norm"\t"confInt(1.6)"\t"cstat"\t"chi2"\t"dof"\t"PhoInd_2_10"\t"confInt(1.6σ)_2_10"\t"kT_2_10"\t"confInt(1.6σ)_2_10"\t"kT_norm_2_10"\t"confInt_2_10(1.6σ)"\t"absflux_2_10"\t"flux_confInt_2_10(1.6σ)"\t"Lum_2_10"\t"lum_confInt_2_10(1.6σ)"\t"flux_confInt_2_10(3σ)"\t"lum_confInt_2_10(3σ)"\t"PhoInd_3_78"\t"confInt(1.6σ)_3_78"\t"kT_3_78"\t"confInt(1.6σ)_3_78"\t"kT_norm_3_78"\t"confInt_3_78(1.6σ)"\t"absflux_3_78"\t"flux_confInt_3_78(1.6σ)"\t"Lum_3_78"\t"lum_confInt_3_78(1.6σ)"\t"flux_confInt_3_78(3σ)"\t"lum_confInt_3_78(3σ)"\t"PhoInd_2_78"\t"confInt(1.6σ)_2_78"\t"kT_2_78"\t"confInt(1.6σ)_2_78"\t"kT_norm_2_78"\t"confInt_2_78(1.6σ)"\t"absflux_2_78"\t"flux_confInt_2_78(1.6σ)"\t"Lum_2_78"\t"lum_confInt_2_78(1.6σ)"\t"flux_confInt_2_78(3σ)"\t"lum_confInt_2_78(3σ)"\t"PhoInd_10_30"\t"confInt(1.6σ)_10_30"\t"kT_10_30"\t"confInt(1.6σ)_10_30"\t"kT_norm_10_30"\t"confInt_10_30(1.6σ)"\t"absflux_10_30"\t"flux_confInt_10_30(1.6σ)"\t"Lum_10_30"\t"lum_confInt_10_30(1.6σ)"\t"flux_confInt_10_30(3σ)"\t"lum_confInt_10_30(3σ)"\t"PhoInd_10_78"\t"confInt(1.6σ)_10_78"\t"kT_10_78"\t"confInt(1.6σ)_10_78"\t"kT_norm_10_78"\t"confInt_10_78(1.6σ)"\t"absflux_10_78"\t"flux_confInt_10_78(1.6σ)"\t"Lum_10_78"\t"lum_confInt_10_78(1.6σ)"\t"flux_confInt_10_78(3σ)"\t"lum_confInt_10_78(3σ)"\n')



# Start of the actuall program. Loop over each obsid
# ––––––––––––––––––––––––––––––––––––––––––––––––––
for i in range(len(obsid)): 

    # remove SN in name
    name = name_file[i].replace(' ','')
    
    #if str(obsid[i]) in os.listdir(wdir):
        #if 'out' in os.listdir(f'{wdir}/{obsid[i]}/'):
            #if str(name) in os.listdir(f'{wdir}/{obsid[i]}/out/'):
    # go into each repro file of the corresponding obsid
    list_files_in_repro = os.listdir(f'{wdir}/{obsid[i]}/out/{name}')
    spec_dir = f'{wdir}/{obsid[i]}/out/{name}'

    # match nH value with SN
    index_nH = np.where(np.array(string_nH)=='SN'+name)[0]
    nH = float(np.array(nH_all)[index_nH][0])/1e22

    RA = np.array(nH_file['SN_RA'])[index_nH]
    DEC = np.array(nH_file['SN_DEC'])[index_nH]
    
    # find the spectral file
    grp_FPMA = [file for file in list_files_in_repro if file.endswith(f'nu{obsid[i]}A01_sr_grp.pha')][0]
    grp_FPMB = [file for file in list_files_in_repro if file.endswith(f'nu{obsid[i]}B01_sr_grp.pha')][0]
    
    spec_file_A = f'{spec_dir}/{grp_FPMA}'
    spec_file_B = f'{spec_dir}/{grp_FPMB}'
    
    current = os.getcwd()  # current folder we are in
    os.chdir(spec_dir) # same as cd but for "python environment"
    
    # if spectral file does exist
    if (grp_FPMA in list_files_in_repro):
        
        # write out on the list
        f.write(name+'\t')
        f.write(str(obsid[i])+'\t')
        f.write(str(dist[i])+'\t')
        f.write(str(type[i])+'\t')
        f.write(str(epoch[i])+'\t')

        
        # load spectrum
        xspec.AllData(' 1:1 '+spec_file_A + ' 2:2 '+ spec_file_B)
        
        s1 = AllData(1)
        s2 = AllData(2)
        
        print(s1.rate)
        countsA = float(s1.exposure)*float(s1.rate[0])
        countsB = float(s2.exposure)*float(s2.rate[0])
        
        f.write(str(round(s1.exposure/1e3,2))+'\t')
        f.write(str(round(s2.exposure/1e3,2))+'\t')
        f.write(str(nH)+'\t')
        f.write(str(RA[0])+'\t')
        f.write(str(DEC[0])+'\t')
        f.write(str(date[i].split('T')[0])+'\t')
        
        xspec.AllData(' 1:1 '+spec_file_A + ' 2:2 '+ spec_file_B)
        xspec.AllData.ignore("**-3.0 78.0-**")
        xspec.AllData.ignore("bad")
        
        f.write(str(round(countsA,2)) +'\t')
        f.write(str(round(countsB,2)) +'\t')
        
        if ma.is_masked(obsid_file['obsflux_2_10'].data[i]) or obsid_file['obsflux_2_10'].data[i]==' ':
            if i==len(obsid)-1:
                f.write(" \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t")
            else:
                f.write(" \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n")

        else:

            mo = None
            if float(obsid_file['flux_confInt_2_10(3σ)'][i].split(',')[0].replace('(',''))==0.0 or float(obsid_file['flux_confInt_2_10(3σ)'][i].split(',')[0].replace('(',''))==1.0:
                detected='no'
                print('Not detected, aka no mekal')
                
                # apply model
                mo = Model("const*tbabs*(po)")
            else:
                detected='yes'
                mo = Model("const*tbabs*(po+mekal)")
                
                # interaction
                mo(5).values = [0.4, 0.01, 0.1, 0.1, 15.0, 15.0]
                mo(6).values = 1
                mo(7).values = 1
                mo(8).values = 0
                mo(9).values = 1
                mo(10).values = [1e-6, 1e-9, 1e-20, 1e-20, 1e-1,1e-1]
                
            # define model parameters
            # –––––––––––––––––––––––
          
            # const
            mo.setPars({1:"1.0 -.1"})
                                    
            # tbabs
            mo(2).values = nH
            mo(2).frozen = True
            
            # powerlaw
            mo(3).values = [2.0, 0.01, 0.5, 0.5, 3.0, 3.0]  # value, fit delta, min, bot, top, max
            mo(4).values = [1e-6, 1e-9, 1e-20, 1e-20, 1e2,1e2]
            
            
            # second model
            m2 = AllModels(2)
            m2(1).values = [1, 0.01, 0.9, 0.9, 1.1, 1.1 ]

                
            # here we start the fit of the model for 90% confidence
            # –––––––––––––––––––––––––––––––––––––––––––––––––––––
            try:
                # fit the model to the spectrum
                Fit.perform()
                Fit.error("2.706 3")
                
                gamma_err = mo(3).error
                gamma_low , gamma_up = gamma_err[0], gamma_err[1]
                
                # photon index
                gamma = mo(3).values[0]
                
                # constrains so that photon index is not skyhigh or deepsea low
                if mo(3).values[0]<0.51 or mo(3).values[0]>2.99:
                    mo(3).values = 2
                    mo(3).frozen = True
                    
                    gamma = 2
                    gamma_low , gamma_up = 0,0
                    
                    Fit.perform()
                    

                if detected=='yes':

                    Fit.error("2.706 5")
                    kT = mo(5).values[0]
                    kT_err = mo(5).error
                    kT_low, kT_up = kT_err[0], kT_err[1]
                    kT_2= mo(6).values[0]
                    kT_3 = mo(7).values[0]
                    kT_4 = mo(8).values[0]
                    kT_5 = mo(9).values[0]
                    Fit.error("2.706 10")
                    kT_norm = mo(10).values[0]
                    kT_norm_err = mo(10).error
                    kT_norm_low, kT_norm_up = kT_norm_err[0], kT_norm_err[1]

                
                # extract values from fit
                # –––––––––––––––––––––––
                const1 = mo(1).values[0]
                const2 = m2(1).values[0]
                

                # normalinsation of powerlaw
                Fit.error("2.706 4")  #cstat 2.706=1σ double limit on the second parameter
                norm = mo(4).values[0]
                norm_err = mo(4).error
                gamma_norm_low, gamma_norm_up = norm_err[0], norm_err[1]

                
                # cstatistics info
                cstat = xspec.Fit.statistic  # cstat
                chi2 = xspec.Fit.testStatistic  # chi2 stat.
                dof = xspec.Fit.dof  # degrees of freedom = number of data points-number of free parameters

                # write out the results
                f.write(f'{round(gamma,4)}\t')
                f.write(f'({round(gamma_low,3)},{round(gamma_up,3)})\t')
                f.write(f"{format(norm,'.4e')}\t")
                f.write(f"({format(gamma_norm_low,'.4e')},{format(gamma_norm_up,'.4e')})\t")
                if detected=='yes':
                    f.write(f'{round(kT,4)}\t')
                    f.write(f'({round(kT_low,3)},{round(kT_up,3)})\t')
                    f.write(f"{format(kT_norm,'.4e')}\t")
                    f.write(f"({format(kT_norm_low,'.4e')},{format(kT_norm_up,'.4e')})\t")
            
                else:
                    f.write('0\t')
                    f.write('(0.0,0.0)\t')
                    f.write('0\t')
                    f.write('(0.0,0.0)\t')
                
                f.write(f'{round(cstat,2)}\t')
                f.write(f'{round(chi2,2)}\t')
                f.write(f'{dof}\t')
                try1 = True
            except:
                try1 = False
                if i==len(obsid)-1:
                    f.write(" \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t")
                else:
                    f.write(" \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n")
                
            if try1==True:
                # save the models for each line of sight within the repro folder of each obsid
                model_name = f'model_{name}_{obsid[i]}_noAbs.xcm'
                if (os.path.isfile(model_name)):
                    os.remove(model_name)
                Xset.save(f'model_{name}_{obsid[i]}_noAbs.xcm', info='m')
                
                save_data = f'./specfit_delchi_{name}_{obsid[i]}_noAbs.qdp'
                if (os.path.isfile(save_data)):
                    os.remove(save_data)

                xspec.Plot.addCommand(f'wd {save_data}')
                Plot.setRebin(minSig=4, maxBins=4)
                Plot.setGroup("1,2")
                #Plot.add=False
                xspec.Plot.device = '/xs'
                xspec.Plot.xAxis = "keV"
                xspec.Plot('ldata','model','del')
                xspec.Plot.commands = ()

                
                # plot and save figure
                nustar_3_window_plot(name,obsid[i])
                
                for j in range(5):
                    mo = None
                    if detected=='yes':
                        mo = Model("const*tbabs*(cflux*po+mekal)")
                        mo(8).values = [kT, 0.01, 0.1, 0.1, 15.0, 15.0]
                        mo(8).frozen = True
                        mo(9).values = kT_2
                        mo(10).values = kT_3
                        mo(11).values = kT_4
                        mo(12).values = kT_5
                        mo(13).values = [kT_norm, 1e-9, 1e-20, 1e-20, 1e-1,1e-1]
                        mo(13).frozen = True
                        
                    else:
                        mo = Model("const*tbabs*(cflux*po)")
                    
                    mo(1).values = const1
                    mo(1).frozen = True
                    
                    m2 = AllModels(2)
                    m2(1).values = const2
                    m2(1).frozen = True
                        
                    # tbabs
                    mo(2).values = nH
                    mo(2).frozen = True

                    # cflux
                    mo(3).values = E_int[j][0]
                    mo(3).frozen = True
                    mo(4).values = E_int[j][1]
                    mo(4).frozen = True
                    mo(5).values = -13
                    
                    mo(6).values =  gamma
                    mo(6).frozen = True
                    mo(7).values = norm
                    mo(7).frozen = True
                
                    
                    try:
                        Fit.perform()
                        if mo(6).values[0]!=2.0:
                            mo(6).frozen = False
                        if detected=='yes':
                            mo(8).frozen = False
                            mo(13).frozen = False
                    
                        Fit.perform()
                    
                        Fit.error("2.706 5")
                        cflux = mo(5).values[0]
                        cflux_err = mo(5).error
                        test=1
                        cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                        
                        # if PL index is out of range:
                        if mo(6).values[0]<=0.51 or mo(6).values[0]>2.99:
                            mo(6).values = 2
                            mo(6).frozen = True
                            
                            Fit.perform()
                            
                            cflux = mo(5).values[0]
                            cflux_err = mo(5).error
                            cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                
                        # if flux is an even number, do steppar instead
                        if (cflux % 1==0) or cflux<-20:
                            cflux, cflux_under, cflux_over = call_steppar(cstat=2.706,paramNr=5)
                        flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                        
                    except:
                        try:
                            cflux, cflux_under, cflux_over = call_steppar(cstat=2.706,paramNr=5)
                            flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                        except:

                            flux, flux_under, flux_over = 0,0,0
                                    
                    gamma2 = mo(6).values[0]
                    gamma_err2 = mo(6).error
                    gamma_low2, gamma_up2 = gamma_err2[0], gamma_err2[1]
                    f.write(f'{round(gamma2,4)}\t')
                    f.write(f'({round(gamma_low2,3)},{round(gamma_up2,3)})\t')

                    
                    if detected=='yes':
                        Fit.error("2.706 8") # 1σ error
                        kT2 = mo(8).values[0]
                        kT_err2 = mo(8).error
                        kT2_low2, kT2_up2 = kT_err2[0], kT_err2[1]
                        
                        Fit.error("2.706 13") # 1σ error
                        kT_norm2 = mo(13).values[0]
                        kT_norm_err2 = mo(13).error
                        kT_norm_low2, kT_norm_up2 = kT_norm_err2[0], kT_norm_err2[1]
                        
                        f.write(f'{round(kT2,4)}\t')
                        f.write(f'({round(kT2_low2,3)},{round(kT2_up2,3)})\t')
                        f.write(f"{format(kT_norm2,'.4e')}\t")
                        f.write(f"({format(kT_norm_low2,'.4e')},{format(kT_norm_up2,'.4e')})\t")
                
                    else:
                        f.write('0\t')
                        f.write('(0.0,0.0)\t')
                        f.write('0\t')
                        f.write('(0.0,0.0)\t')
        
                    
                    
                    L, L_min, L_max = lum_funk(flux,flux_under,flux_over,dist[i])
                                                                
                    f.write(f"{format(flux,'.4e')}\t")
                    f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                    f.write(f"{format(L,'.4e')}\t")
                    f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")
                    
                    try:
                        Fit.error("9.0 5")
                        cflux = mo(5).values[0]
                        cflux_err = mo(5).error
                        cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                        test2=1
                        
                        if cflux_under == 0:
                            Fit.error("7.74 5")
                            cflux = mo(5).values[0]
                            cflux_err = mo(5).error
                            cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                            
                        # if flux is an even number, do steppar instead
                        if (cflux % 1==0) or cflux<-20:
                            cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=5)
                            
                            if cflux_under==1:
                                cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=5)

                        flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                    except:
                        try:
                            cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=5)

                            if cflux_under==1:
                                cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=5)
                            flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                        except:
                            
                            flux, flux_under, flux_over = 0,0,0
            
                    L, L_min, L_max = lum_funk(flux,flux_under,flux_over,dist[i])
                            
                        
                    # write out the results
                    f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                    if j==4:
                        f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\n")
                    else:
                        f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")
                    
                       
    mo = None
    xspec.AllData.clear()
    s = None
                    
os.chdir(base)


f.close()









