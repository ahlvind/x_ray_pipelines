import sys
import os

# Set up and files to read
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should giv a string like '/Users/juliaahlvind'
wdir = f'{home}/Documents/projekt_1/NuSTAR/data'  # path tp mother directory that includes all obsid folders
base = f'{home}/Documents/projekt_1/NuSTAR'  # path tp mother directory that includes all obsid folders
projekt_1 = f'{home}/Documents/projekt_1'
folder = '22feb'

sys.path.insert(0,f'{projekt_1}/pipelines')
from pyxspec_functions import *

# setting up cstatistics and limited iterations
Fit.nIterations = 150
Fit.statMethod = "cstat"
Fit.query = "yes" # run fitting until this limit:Fit.nIterations is reached


# Chandra file with SN names and obsid
obsid_file = ascii.read(f'{projekt_1}/nustar_to_download.txt')
#obsid_file = ascii.read(f'./test.txt')
obsid = obsid_file['obsid']
name_file = obsid_file['Name']
dist = obsid_file['dist[Mpc]']
type = obsid_file['Type']
epoch = obsid_file['epoch']
expt = obsid_file['exp_time(ks)']
date = obsid_file['nustar_date']
# nH file to cross check with

nH_file = ascii.read(f"{projekt_1}/nH_60Mpc.csv", delimiter=';')
nH_all = nH_file['NH_tot_(Mean)_[atoms_cm-2]']
nH_names = nH_file['Name']
string_nH=[]
for k in nH_names:
    string_nH.append(k)
    
Mpc_cm = 3.08567758*1e24  # mult with distance
kev_erg = 1/624150964.7  # mult with integral output keV


E_int = [[2.0,10.0],[7.0,10.0],[2.0,78.0],[10.0,30.0], [10.0,78.0]]

    
f = open(f'{projekt_1}/{folder}/ObsFlux_test2.csv','+w')  # creates a csv document with all outputs
f.write('"Name"\t"obsid"\t"dist(Mpc)"\t"type"\t"epoch"\t"exp_timeA(ks)"\t"exp_timeB(ks)"\t"mean_nH(1e22atoms_cm-2)"\t"RA"\t"DEC"\t"dateXray"\t"countsA"\t"countsB"\t"frac_counts_src_to_tot_A(%)"\t"frac_counts_src_to_tot_B(%)"\t"gamma"\t"confInt(1.6σ)"\t"norm(1.6σ)"\t"confInt(1.6σ)"\t"constant"\t"cstat"\t"chi2"\t"dof"\t"obsflux_2_10"\t"flux_confInt_2_10(1.6σ)"\t"Lum_2_10"\t"lum_confInt_2_10(1.6σ)"\t"flux_confInt_2_10(3σ)"\t"lum_confInt_2_10(3σ)"\t"obsflux_3_78"\t"flux_confInt_3_78(1.6σ)"\t"Lum_3_78"\t"lum_confInt_3_78(1.6σ)"\t"flux_confInt_3_78(3σ)"\t"lum_confInt_3_78(3σ)"\t"obsflux_2_78"\t"flux_confInt_2_78(1.6σ)"\t"Lum_2_78"\t"lum_confInt_2_78(1.6σ)"\t"flux_confInt_2_78(3σ)"\t"lum_confInt_2_78(3σ)"\t"obsflux_10_30"\t"flux_confInt_10_30(1.6σ)"\t"Lum_10_30"\t"lum_confInt_10_30(1.6σ)"\t"flux_confInt_10_30(3σ)"\t"lum_confInt_10_30(3σ)"\t"obsflux_10_78"\t"flux_confInt_10_78(1.6σ)"\t"Lum_10_78"\t"lum_confInt_10_78(1.6σ)"\t"flux_confInt_10_78(3σ)"\t"lum_confInt_10_78(3σ)"\n')

# Start of the actuall program. Loop over each obsid
# ––––––––––––––––––––––––––––––––––––––––––––––––––
for i in range(13):#len(obsid)):
    # remove SN in name
    name = name_file[i].replace(' ','')
    
    
    if str(obsid[i]) in os.listdir(wdir):
        if 'out' in os.listdir(f'{wdir}/{obsid[i]}/'):
            if str(name) in os.listdir(f'{wdir}/{obsid[i]}/out/'):
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
                    #s1 = Spectrum(spec_file_A)
                    #s2 = Spectrum(spec_file_B)
                    
                    xspec.AllData(' 1:1 '+spec_file_A + ' 2:2 '+ spec_file_B)
                    
                    s1 = AllData(1)
                    s2 = AllData(2)

                    percentageA = round(float(s1.exposure)*float(s1.rate[0])/(float(s1.exposure)*float(s1.rate[2]))*100,2)
                    percentageB = round(float(s2.exposure)*float(s2.rate[0])/(float(s2.exposure)*float(s2.rate[2]))*100,2)
                                        
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

                    #s1 = AllData(1)
                    
                    f.write(str(round(countsA,2)) +'\t')
                    f.write(str(round(countsB,2)) +'\t')
                    f.write(str(percentageA) +'\t')
                    f.write(str(percentageB) +'\t')
                    
                    print(f'Counts: A:{countsA} and B:{countsB}')
                    if countsA<1 and countsB<1:
                        if i==len(obsid)-1:
                            f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t ") #48
                        else:
                            f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n ") #48
                    else:
                    
                        
                        mo = None
                        # apply model
                        mo = Model("const*tbabs*(po)")
                        XspecSettings.chatter = 0
                        
                        # define model parameters
                        # –––––––––––––––––––––––
                        
                        xspec.AllData.clear()
                        s = None
                        xspec.AllData(' 1:1 '+spec_file_A + ' 2:2 '+ spec_file_B)
                        xspec.AllData.ignore("**-3.0 78.0-**")
                        xspec.AllData.ignore("bad")
                        
                        # const
                        mo.setPars({1:"1.0 -.1"})
                                                
                        # tbabs
                        mo(2).values = nH
                        mo(2).frozen = True
                        
                        # powerlaw
                        mo(3).values = [2.0, 0.01, 0.5, 0.5, 3.0, 3.0]
                        mo(4).values = [1e-6, 1e-9, 1e-20, 1e-20, 1e2,1e2]

                        # second model
                        m2 = AllModels(2)
                        m2(1).values = [1, 0.01, 0.9, 0.9, 1.1, 1.1]
                            
                        
                        # here we start the fit of the model for 90% confidence
                        # –––––––––––––––––––––––––––––––––––––––––––––––––––––
                        try:
                            # fit the model to the spectrum
                            Fit.perform()
                            
                            # extract values from fit
                            # –––––––––––––––––––––––
                            const1 = mo(1).values[0]
                            const2 = m2(1).values[0]
                                
                            # extract values from fit
                            # –––––––––––––––––––––––
                            Fit.error("2.706 3") # 1σ error

                            gamma = mo(3).values[0]
                            gamma_err = mo(3).error
                            gamma_low, gamma_up = gamma_err[0], gamma_err[1]
                            
                            # constrains so that photon index is not skyhigh or deepsea low
                            if mo(3).values[0]<0.51 or mo(3).values[0]>2.99 or gamma_up==0:
                                mo(3).values = 2
                                mo(3).values[0] = 2
                                mo(3).frozen = True
                                                    
                                Fit.perform()
                                gamma = mo(3).values[0]
                                gamma_low, gamma_up = 0,0
                            
                            
                            Fit.error("2.706 4") # 1σ error
                            gamma_norm = mo(4).values[0]
                            gamma_norm_err = mo(4).error
                            gamma_norm_low, gamma_norm_up = gamma_norm_err[0], gamma_norm_err[1]
                    
                            # cstatistics info
                            cstat = xspec.Fit.statistic  # cstat
                            chi2 = xspec.Fit.testStatistic  # chi2 stat.
                            dof = xspec.Fit.dof  # degrees of freedom = number of data points-number of free parameters

                            # write out the results
                            f.write(f'{round(gamma,4)}\t')
                            f.write(f'({round(gamma_low,3)},{round(gamma_up,3)})\t')
                            f.write(f"{format(gamma_norm,'.4e')}\t")
                            f.write(f"({format(gamma_norm_low,'.4e')},{format(gamma_norm_up,'.4e')})\t")
                            f.write(f'{round(const2,2)}\t')
                            f.write(f'{round(cstat,2)}\t')
                            f.write(f'{round(chi2,2)}\t')
                            f.write(f'{dof}\t')
                            try1 = True
                        except:
                            try1 = False
                            if i==len(obsid)-1:
                                f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t ") #48
                            else:
                                f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n ") #48
                        
                        if try1==True:
                            # save the models for each line of sight within the repro folder of each obsid
                            model_name = f'model_{name}_{obsid[i]}_obsFlux.xcm'
                            if (os.path.isfile(model_name)):
                                os.remove(model_name)
                            Xset.save(f'model_{name}_{obsid[i]}_obsflux.xcm', info='m')
                            
                            #if sing_double=='double':
                            #    Plot.setGroup("1-2")
                            #    xspec.Plot.add= True
                            #Plot.add=False
                            xspec.Plot.device = '/xs'
                            print('testtest 2.7')
                            xspec.Plot.xAxis = "keV"
                            xspec.Plot('ldata','del')
                            # test to rebin just so it will be easier to see in the fit
                            #xspec.Plot.setRebin(minSig=20.0,maxBins=20)
                            xspec.Plot.commands = ()
                            
                            
                            print('testtest 3')
                            save_data = f'specfit_delchi_{name}_{obsid[i]}_obsFlux.qdp'
                            if (os.path.isfile(save_data)):
                                os.remove(save_data)
  
                            """
                            print('testtest 4')
                            
                            #if sing_double=='double':
                            #    Plot.setGroup("1-2")
                            #    xspec.Plot.add = True
                            xspec.Plot.device = '/null'
                            xspec.Plot.addCommand(f'wd {save_data}')
                            xspec.Plot("ld","del")
                            # test to rebin just so it will be easier to see in the fit
                            #xspec.Plot.setRebin(minSig=20.0,maxBins=20)
                            names = ['e','de','rate','rate_err','total']  # energi(x-axel), step in energy, y-axel data, y-axel error data, totala model
                            ncomp = len(mo.componentNames)
                            for a in range(ncomp):
                                names.append(f'model{a}')
                            df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsFlux.qdp',skiprows=3,names=names, delimiter=' ')
                            index = np.where(df.e=='NO')[0]
                            df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsFlux.qdp', skiprows=3, names=names, delimiter=' ',nrows=index[0])
  
                            print('index', index)
                            fig = plt.figure(1)
                            ax = fig.add_subplot(211)
                            # Plot using Matplotlib:
                            # plott first spectra
                            ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err,fmt='.',label='data', markersize=1, color='black')
                            ax.step(df.e, df.total, color='black',label='Total model',linewidth=2)
                                                        
                            # plot second spectra
                            df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsFlux.qdp',skiprows=index[0]+4, names=names, delimiter=' ',nrows=index[1]-index[0]-1)
                            ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err,fmt='.',label='data2', color='red')
                            ax.step(df.e, df.total, color='red',label='Total model2',linewidth=2)
                            ax.set_ylabel(r'count$~$s$^{-1}~$keV$^{-1}$')
                            ax.set_xscale("linear")
                            ax.set_yscale("log")
                            ax.legend()
                
                            
                            df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsFlux.qdp', skiprows=index[2]+4, names=names, delimiter=' ',nrows=index[3]-index[2]-1)
                            ax = fig.add_subplot(2,1,2)
                            ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err, fmt='.', color='black')
                            
                            df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsFlux.qdp', skiprows=index[3]+4, names=names, delimiter=' ',nrows=index[4]-index[3]-1)
                            ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err, fmt='.', color='red')
                            ax.set_ylabel('(data-model)/error')
                            ax.set_xlabel('Energy (keV)')
                            ax.axhline(y=0, linestyle='--', linewidth=1, color='black')
                            
                            #plt.show()
                            
                            fig.savefig(f'/Users/juliaahlvind/Documents/projekt_1/NuSTAR/fits/obsflux/fit_{name}_{obsid[i]}_obsFlux')
                            fig.clf()
                            
                            Plot.commands = ()
                            """
                       
                            for j in range(5):
                                mo = None
                                mo = Model("cflux*(const*tbabs*po)")

                                mo(4).values = 1
                                mo(4).frozen = True
                                      
                                # cflux
                                mo(1).values = E_int[j][0]
                                mo(2).values = E_int[j][1]
                                mo(3).values = -13
                                
                                # tbabs
                                mo(5).values = nH
                                mo(5).frozen = True

                                # powerlaw
                                mo(6).values = [2.0, 0.01, 0.5, 0.5, 3.0, 3.0]
                                mo(6).frozen = True
                                mo(7).values = [1e-6, 1e-9, 1e-20, 1e-20, 1e2,1e2]
                                mo(7).frozen = True
                                
                                m2 = AllModels(2)
                                m2(4).values = const2#[1, 0.01, 0.9, 0.9, 1.1, 1.1 ]
                                m2(4).frozen = True
                                
                                
                                try:
                                    print('first try')
                                    Fit.perform()
                                    mo(6).frozen = False
                                    Fit.perform()
                                    
                                    if mo(6).values[0]<0.51 or mo(6).values[0]>2.99 or gamma_up==0:
                                        mo(6).values = 2
                                        mo(6).values[0] = 2
                                        mo(6).frozen = True
                                                            
                                        Fit.perform()
                                
                                    if mo(3).values[0]<-20: # try fitting again
                                        mo(3).values = -15
                                        Fit.perform()
                                    
                                    Fit.error("2.706 3")
                                    cflux = mo(3).values[0]
                                    cflux_err = mo(3).error
                                    cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                                    
                                    # if flux is an even number, do steppar instead
                                    if (cflux % 1==0) or cflux<-20:
                                        try: # try doing steppar if fitting didn't work
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=2.706,paramNr=3,step_start=-18, step_end=-8)
                                        except:
                                            cflux, cflux_under, cflux_over = 0,0,0

                                except:
                                    
                                    mo(6).frozen = False
                                    try:
                                        print('second try')
                                        cflux, cflux_under, cflux_over = call_steppar(cstat=2.706,paramNr=3,step_start=-18, step_end=-8)
                                    except:
                                        cflux, cflux_under, cflux_over = 0,0,0
                                        crach_par = 0

                                
                                flux, flux_under, flux_over = 10**cflux, 10**cflux_under, 10**cflux_over
                                
                                # calculate luminosity
                                L, L_min, L_max = lum_funk(flux,flux_under,flux_over,dist[i])
                                
                                f.write(f"{format(flux,'.4e')}\t")
                                f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                                f.write(f"{format(L,'.4e')}\t")
                                f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")
                                
                                try:
                                    Fit.error("9.0 3")
                                    cflux = mo(3).values[0]
                                    cflux_err = mo(3).error
                                    cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                                        
                                                    
                                    if cflux_under == 0 or cflux_under == 1:
                                        Fit.error("7.74 3")
                                        cflux = mo(3).values[0]
                                        cflux_err = mo(3).error
                                        cflux_under, cflux_over = 0, cflux_err[1]
                                    
                                    # if flux is an even number, do steppar instead
                                    if (cflux % 1==0) or cflux<-20:
                                        cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=3,step_start=-18, step_end=-8)
                                        
                                        if cflux_under==0 or cflux_under == 1:
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=3,step_start=-18, step_end=-8)
                                            cflux_under = 0

                                    flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                                
                                # if fit not possible to steppar directly
                                except:
                                    try:
                                        cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=3,step_start=-18, step_end=-8)

                                        if cflux_under==0 or cflux_under==1:
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=3,step_start=-18, step_end=-8)
                                            cflux_under = 0
                                        flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                                    except:
                                        flux, flux_under, flux_over = 0,0,0
                                            
                                # calculate luminosity
                                L, L_min, L_max = lum_funk(10**cflux,10**cflux_under,10**cflux_over,dist[i])
                                    
                                # write out the results
                                if flux==1 or flux==0:
                                    if j==4:
                                        f.write('\t \n')
                                    else:
                                        f.write('\t \t')
                                else:
                                    f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                                    if j==4:
                                        f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\n")
                                    else:
                                        f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")

                                
                                #if j==4:
                                #    f.write('\t\t\t\n')
                                #else:
                                #    f.write('\t\t\t\t')
    mo = None
    xspec.AllData.clear()
    s = None

os.chdir(base)


f.close()









