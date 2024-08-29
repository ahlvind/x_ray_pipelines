import sys
import os
import numpy.ma as ma

# Set up and files to read
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should giv a string like '/Users/juliaahlvind'
wdir = f'{home}/Documents/projekt_1/NuSTAR/data'  # path tp mother directory that includes all obsid folders
base = f'{home}/Documents/projekt_1/NuSTAR'  # path tp mother directory that includes all obsid folders
projekt_1 = f'{home}/Documents/projekt_1/'
folder = '22feb'

sys.path.insert(0,f'{projekt_1}/pipelines')
from pyxspec_functions import *

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
# nH file to cross check with

nH_file = ascii.read(f"{projekt_1}/nH_60Mpc.csv", delimiter=';')
nH_all = nH_file['NH_tot_(Mean)_[atoms_cm-2]']
nH_names = nH_file['Name']
string_nH=[]
for k in nH_names:
    string_nH.append(k)
    
Mpc_cm = 3.08567758*1e24  # mult with distance
kev_erg = 1/624150964.7  # mult with integral output keV

line_of_sight_string = ["10th", "median", "90th"]

E_int = [[2.0,10.0],[7.0,10.0],[2.0,78.0],[10.0,30.0], [10.0,78.0]]

f = open(f'{projekt_1}/{folder}/tbvarabs_test.csv','+w')  # creates a csv document with all outputs
f.write('"Name"\t"obsid"\t"dist[Mpc]"\t"type"\t"epoch"\t"exp_timeA(ks)"\t"exp_timeB(ks)"\t"mean_nH(1e22atoms_cm-2)"\t"RA"\t"DEC"\t"dateXray"\t"countsA"\t"countsB"\t"line_of_sight"\t"PhoIndex"\t"confInt(1.6σ)"\t"norm"\t"confInt(1.6σ)"\t"kT"\t"confInt(1.6σ)"\t"norm"\t"confInt(1.6)"\t"cstat"\t"chi2"\t"dof"\t"PhoInd_2_10"\t"confInt(1.6σ)_2_10"\t"PL_norm_2_10"\t"kT_2_10"\t"confInt(1.6σ)_2_10"\t"kT_norm_2_10"\t"confInt_2_10(1.6σ)"\t"absflux_2_10"\t"flux_confInt_2_10(1.6σ)"\t"Lum_2_10"\t"lum_confInt_2_10(1.6σ)"\t"flux_confInt_2_10(3σ)"\t"lum_confInt_2_10(3σ)"\t"PhoInd_3_78"\t"confInt(1.6σ)_3_78"\t"PL_norm_3_78"\t"kT_3_78"\t"confInt(1.6σ)_3_78"\t"kT_norm_3_78"\t"confInt_3_78(1.6σ)"\t"absflux_3_78"\t"flux_confInt_3_78(1.6σ)"\t"Lum_3_78"\t"lum_confInt_3_78(1.6σ)"\t"flux_confInt_3_78(3σ)"\t"lum_confInt_3_78(3σ)"\t"PhoInd_2_78"\t"confInt(1.6σ)_2_78"\t"PL_norm_2_78"\t"kT_2_78"\t"confInt(1.6σ)_2_78"\t"kT_norm_2_78"\t"confInt_2_78(1.6σ)"\t"absflux_2_78"\t"flux_confInt_2_78(1.6σ)"\t"Lum_2_78"\t"lum_confInt_2_78(1.6σ)"\t"flux_confInt_2_78(3σ)"\t"lum_confInt_2_78(3σ)"\t"PhoInd_10_30"\t"confInt(1.6σ)_10_30"\t"PL_norm_10_30"\t"kT_10_30"\t"confInt(1.6σ)_10_30"\t"kT_norm_10_30"\t"confInt_10_30(1.6σ)"\t"absflux_10_30"\t"flux_confInt_10_30(1.6σ)"\t"Lum_10_30"\t"lum_confInt_10_30(1.6σ)"\t"flux_confInt_10_30(3σ)"\t"lum_confInt_10_30(3σ)"\t"PhoInd_10_78"\t"confInt(1.6σ)_10_78"\t"PL_norm_10_78"\t"kT_10_78"\t"confInt(1.6σ)_10_78"\t"kT_norm_10_78"\t"confInt_10_78(1.6σ)"\t"absflux_10_78"\t"flux_confInt_10_78(1.6σ)"\t"Lum_10_78"\t"lum_confInt_10_78(1.6σ)"\t"flux_confInt_10_78(3σ)"\t"lum_confInt_10_78(3σ)"\n') #90

# Start of the actuall program. Loop over each obsid
# ––––––––––––––––––––––––––––––––––––––––––––––––––
for i in range(len(obsid)):

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
                    
                    xspec.AllData(f' 1:1 {spec_file_A} 2:2 {spec_file_B}')
                    xspec.AllData.ignore("**-3.0 78.0-**")
                    xspec.AllData.ignore("bad")
                    tbvarabs_all = tbvarabs_funk(name,obsid[i],type[i],epoch[i])

                    
                    f.write(str(round(countsA,2)) +'\t')
                    f.write(str(round(countsB,2)) +'\t')
                    for j in range(3):
                        if j==0:
                            f.write(line_of_sight_string[j]+'\t')
                        else:
                            f.write(f'\t \t \t \t \t \t \t \t \t \t \t \t \t {line_of_sight_string[j]} \t')
                    
                        if ma.is_masked(obsid_file['obsflux_2_10'].data[i]) or obsid_file['obsflux_2_10'].data[i]==' ':
                            if i==len(obsid)-1 and j==2:#76
                                f.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t")
                            else:
                                f.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n")
                        else:

                            mo = None
                            if float(obsid_file['flux_confInt_2_10(3σ)'][i].split(',')[0].replace('(',''))==0.0 or float(obsid_file['flux_confInt_2_10(3σ)'][i].split(',')[0].replace('(',''))==1.0:
                                detected='no'
                                print('Not detected, aka no mekal')
                                
                                # apply model
                                mo = Model("const*tbabs*(tbvarabs*po)")
                            else:
                                detected = 'yes'
                                mo = Model("const*tbabs*(tbvarabs*po+mekal)")
                                
                                # interaction
                                mo(47).values = [0.4, 0.01, 0.1, 0.1, 15.0, 15.0]
                                mo(48).values = 1
                                mo(49).values = 1
                                mo(50).values = 0
                                mo(51).values = 1
                                mo(52).values = 1e-4
                                
                            # define model parameters
                            # –––––––––––––––––––––––
                        
                            # const
                            mo.setPars({1:"1.0 -.1"})
                                                    
                            # tbabs
                            mo(2).values = nH
                            mo(2).frozen = True
                            
                            # tbvarabs
                            tbvarabs = tbvarabs_all[j]
                            k=0
                            for k in range(len(tbvarabs)):
                                s1, s2, s3, s4, s5, s6 = tbvarabs[k].split()  # each row of tbvarabs aka each element
                                tbvarabs_line = [float(s1), float(s2), float(s3), float(s4), float(s5), float(s6)]
                                mo(k+3).values = tbvarabs_line  # set the model parameter
                                mo(k+3).frozen = True # freezes them

                            # powerlaw
                            mo(45).values = [2.0, 0.01, 0.5, 0.5, 3.0, 3.0]  # value, fit delta, min, bot, top, max
                            mo(46).values = [1e-6, 1e-9, 1e-20, 1e-20, 1e2,1e2]
                            
                            
                            # second model
                            m2 = AllModels(2)
                            m2(1).values = [1, 0.01, 0.9, 0.9, 1.1, 1.1 ]
                                
                                
                            # here we start the fit of the model for 90% confidence
                            # –––––––––––––––––––––––––––––––––––––––––––––––––––––
                            try:
                                # fit the model to the spectrum
                                Fit.perform()
                                Fit.error("2.706 45")
                                
                                # constrains so that photon index is not skyhigh or deepsea low
                                if mo(45).values[0]<0.51 or mo(45).values[0]>2.99:
                                    mo(45).values = 2
                                    mo(45).values[0] = 2
                                    mo(45).frozen = True
                                    
                                    gamma = 2
                                    err_up_gamma, err_low_gamma = 0,0
                                    
                                    Fit.perform()

                                else:
                                    gamma = mo(45).values[0]
                                    Fit.error("2.706 45")
                                    gamma_err = mo(45).error
                                    err_up_gamma, err_low_gamma = gamma_err[0], gamma_err[1]

                                if detected=='yes':
                                    Fit.error("2.706 47")

                                    kT = mo(47).values[0]
                                    kT_err = mo(47).error
                                    kT_err_under, kT_err_over = kT_err[0], kT_err[1]
                                    kT_2= mo(48).values[0]
                                    kT_3 = mo(49).values[0]
                                    kT_4 = mo(50).values[0]
                                    kT_5 = mo(51).values[0]
                                    Fit.error("2.706 52")
                                    kT_norm = mo(52).values[0]
                                    kT_norm_err = mo(52).error
                                    kT_norm_err_under, kT_norm_err_over = kT_norm_err[0], kT_norm_err[1]

                                
                                # extract values from fit
                                # –––––––––––––––––––––––
                                const1 = mo(1).values[0]
                                const2 = m2(1).values[0]

                                # normalinsation of powerlaw
                                Fit.error("2.706 46")  #cstat 2.706=1σ 
                                norm = mo(46).values[0]
                                norm_err = mo(46).error
                                err_low_norm, err_up_norm = norm_err[0], norm_err[1]

                                # cstatistics info
                                cstat = xspec.Fit.statistic  # cstat
                                chi2 = xspec.Fit.testStatistic  # chi2 stat.
                                dof = xspec.Fit.dof  # degrees of freedom = number of data points-number of free parameters

                                # write out the results
                                f.write(str(round(gamma,4))+'\t')
                                f.write(f'±{round(err_up_gamma,4)}\t')
                                f.write(format(norm,'.4e')+'\t')
                                f.write(f"±{format(err_up_norm,'.4e')}\t")
                                if detected=='yes':
                                    f.write(str(round(kT,4))+'\t')
                                    f.write(f"({format(kT_err_under,'.4e')},{format(kT_err_over,'.4e')})\t")
                                    f.write(format(kT_norm,'.4e')+'\t')
                                    f.write(f"({format(kT_norm_err_under,'.4e')},{format(kT_norm_err_over,'.4e')})\t")
                            
                                else:
                                    f.write('0\t')
                                    f.write('(0.0,0.0)\t')
                                    f.write('0\t')
                                    f.write('(0.0,0.0)\t')
                                f.write(str(round(cstat,2))+'\t')
                                f.write(str(round(chi2,2))+'\t')
                                f.write(str(dof)+'\t')

                                try1 = True
                            except:
                                if i==len(obsid)-1 and j==2:
                                    f.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t")
                                else:
                                    f.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n")
                            if try1==True:
                                # save models and plott the fit
                                #model_name = f'model_{name}_{obsid[i]}_{line_of_sight_string[j]}.xcm'
                                """
                                if sing_double=='double':
                                    Plot.setGroup("1-2")
                                #Plot.add=False
                                xspec.Plot.device = '/xs'
                                print('testtest 2.7')
                                xspec.Plot.xAxis = "keV"
                                xspec.Plot('ldata','del')
                                # test to rebin just so it will be easier to see in the fit
                                xspec.Plot.setRebin(minSig=20.0,maxBins=20)
                                xspec.Plot.commands = ()
                                
                                
                                print('testtest 3')
                                save_data = f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp'
                                if (os.path.isfile(save_data)):
                                    os.remove(save_data)
                                
                                
                                print('testtest 4')
                                
                                if sing_double=='double':
                                    Plot.setGroup("1-2")
                                    xspec.Plot.add = True
                                xspec.Plot.device = '/null'
                                xspec.Plot.addCommand(f'wd {save_data}')
                                xspec.Plot("ld","del")
                                # test to rebin just so it will be easier to see in the fit
                                xspec.Plot.setRebin(minSig=20.0,maxBins=20)
                                names = ['e','de','rate','rate_err','total']  # energi(x-axel), step in energy, y-axel data, y-axel error data, totala model
                                ncomp = len(mo.componentNames)
                                for a in range(ncomp):
                                    names.append(f'model{a}')
                                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp',skiprows=3,names=names, delimiter=' ')
                                index = np.where(df.e=='NO')[0]
                                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp', skiprows=3, names=names, delimiter=' ',nrows=index[0])
      
                                
                                fig = plt.figure(1)
                                ax = fig.add_subplot(211)
                                # Plot using Matplotlib:
                                # plott first spectra
                                ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err,fmt='.',label='data', markersize=1, color='black')
                                ax.step(df.e, df.total, color='black',label='Total model',linewidth=2)
                                                            
                                # plot second spectra
                                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp',skiprows=index[0]+4, names=names, delimiter=' ',nrows=index[1]-index[0]-1)
                                ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err,fmt='.',label='data2', color='red')
                                ax.step(df.e, df.total, color='red',label='Total model2',linewidth=2)
                                ax.set_ylabel(r'count$~$s$^{-1}~$keV$^{-1}$')
                                ax.set_xscale("linear")
                                ax.set_yscale("log")
                                ax.legend()
                    
                                
                                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp', skiprows=index[2]+4, names=names, delimiter=' ',nrows=index[3]-index[2]-1)
                                ax = fig.add_subplot(2,1,2)
                                ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err, fmt='.', color='black')
                                
                                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp', skiprows=index[3]+4, names=names, delimiter=' ',nrows=index[4]-index[3]-1)
                                ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err, fmt='.', color='red')
                                ax.set_ylabel('(data-model)/error')
                                ax.set_xlabel('Energy (keV)')
                                ax.axhline(y=0, linestyle='--', linewidth=1, color='black')
                                
                                #plt.show()
                                
                                fig.savefig(f'/Users/juliaahlvind/Documents/projekt_1/NuSTAR_ny/fits/tbvarabs/fit_{name}_{obsid[i]}_{line_of_sight_string[j]}')
                                fig.clf()
                                
                                Plot.commands = ()
                                
                                """
                                
                                for m in range(5):
                                    mo = None
                                    if detected=='yes':
                                        mo = Model("const*tbabs*(tbvarabs*cflux*po+mekal)")
                                        mo(50).values = kT
                                        mo(51).values = kT_2
                                        mo(52).values = kT_3
                                        mo(53).values = kT_4
                                        mo(54).values = kT_5
                                        mo(55).values = kT_norm
                                        variabel = 'detected'
                                        mo(50).frozen = True
                                        mo(55).frozen = True
                                        
                                    else:
                                        mo = Model("const*tbabs*(tbvarabs*cflux*po)")
                                        
                                    mo(1).values = const1
                                    mo(1).frozen = True
                                    
                                    m2 = AllModels(2)
                                    m2(1).values = const2
                                    m2(1).frozen = True
                                          
                                    # tbabs
                                    mo(2).values = nH
                                    mo(2).frozen = True

                                    # tbvarabs
                                    tbvarabs = tbvarabs_all[j]
                                    k=0
                                    for k in range(len(tbvarabs)):
                                        s1, s2, s3, s4, s5, s6 = tbvarabs[k].split()  # each row of tbvarabs aka each element
                                        tbvarabs_line = [float(s1), float(s2), float(s3), float(s4), float(s5), float(s6)]
                                        mo(k+3).values = tbvarabs_line  # set the model parameter
                                        mo(k+3).frozen = True # freezes them

                                    # cflux
                                    mo(45).values = E_int[m][0]
                                    mo(45).frozen = True
                                    mo(46).values = E_int[m][1]
                                    mo(46).frozen = True
                                    mo(47).values = -13
                                    
                                    mo(48).values = gamma
                                    mo(48).frozen = True
                                    mo(49).values = norm
                                    mo(49).frozen = True

                                    try:

                                        Fit.perform()
                                        
                                        if mo(48).values[0]!=2.0:
                                            mo(48).frozen = False
                                        if detected=='yes':
                                            mo(50).frozen = False
                                            mo(55).frozen = False
                                        Fit.perform()
                                        
                                        Fit.error("2.706 47")
                                        cflux = mo(47).values[0]
                                        cflux_err = mo(47).error
                                        test=1
                                        cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                                        
                                        # if PL index is out of range:
                                        if mo(48).values[0]<=0.51 or mo(48).values[0]>2.99:
                                            mo(48).values = 2
                                            mo(48).frozen = True
                                            
                                            Fit.perform()
                                            
                                            cflux = mo(47).values[0]
                                            cflux_err = mo(47).error
                                            cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                                        
                                        # if flux is an even number, do steppar instead
                                        if (cflux % 1==0) or cflux<-20:
                                            print('steg 470')
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=2.706,paramNr=47)
                                            
                                        flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                            
                                    except:
                                        try:
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=2.706,paramNr=47)
                                            flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                                        except:
                                            one_sigma_success = 0
                                            flux, flux_under, flux_over = 0,0,0
                                        
                                    gamma2 = mo(48).values[0]
                                    gamma_err2 = mo(48).error
                                    gamma_low2, gamma_up2 = gamma_err2[0], gamma_err2[1]
                                    f.write(f'{round(gamma2,4)}\t')
                                    f.write(f'({round(gamma_low2,3)},{round(gamma_up2,3)})\t')
                                    f.write(format(mo(49).values[0],'.4e')+'\t')
                                    
                                    if detected=='yes':
                                        Fit.error("2.706 50") # 1σ error
                                        kT2 = mo(50).values[0]
                                        kT_err2 = mo(50).error
                                        kT_err_under2, kT_err_over2 = kT_err2[0], kT_err2[1]
                                        
                                        Fit.error("2.706 55") # 1σ error
                                        kT_norm2 = mo(55).values[0]
                                        kT_norm_err2 = mo(55).error
                                        kT_norm_err_under2, kT_norm_err_over2 = kT_norm_err2[0], kT_norm_err2[1]
                                        
                                        f.write(str(round(kT2,4))+'\t')
                                        f.write(f"({format(kT_err_under2,'.4e')},{format(kT_err_over2,'.4e')})\t")
                                        f.write(format(kT_norm2,'.4e')+'\t')
                                        f.write(f"({format(kT_norm_err_under2,'.4e')},{format(kT_norm_err_over2,'.4e')})\t")
                                
                                    else:
                                        f.write('0\t')
                                        f.write('(0.0,0.0)\t')
                                        f.write('0\t')
                                        f.write('(0.0,0.0)\t')
                                    
                                    
                                    L, L_min, L_max = lum_funk(flux,flux_under,flux_over,dist[i])
                                                                                
                                    f.write(format(flux,'.4e')+'\t')
                                    f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                                    f.write(format(L,'.4e')+'\t')
                                    f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")
                                    one_sigma_success = 1
                                    try:
                                        Fit.error("9.0 47")
                                        cflux = mo(47).values[0]
                                        cflux_err = mo(47).error
                                        cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                                        test2=1
                                        if cflux_under == 0:
                                            Fit.error("7.74 47")
                                            cflux = mo(47).values[0]
                                            cflux_err = mo(47).error
                                            cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                                            
                                        if (cflux % 1==0) or cflux<-20:
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=47)
                                            
                                            if cflux_under==1:
                                                cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=47)

                                        flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                                    except:
                                        try:
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=47)

                                            if cflux_under==1:
                                                cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=47)
                                            flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                                        except:
                                            one_sigma_success = 0
                                            flux, flux_under, flux_over = 0,0,0
                                    L, L_min, L_max = lum_funk(flux,flux_under,flux_over,dist[i])
                                            
                                        
                                    # write out the results
                                    f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                                    if m==4:
                                        f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\n")
                                    else:
                                        f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")
                                        
    mo = None
    xspec.AllData.clear()
    s = None
                    
os.chdir(base)


f.close()


