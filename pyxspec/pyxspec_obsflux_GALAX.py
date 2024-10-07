
import os,sys
from scipy.integrate import quad
import numpy.ma as ma


# INPUT PARAMETERS
#telescope = 'Chandra'
telescope = 'XMM_ny'
folder = '22feb'

# Set up and files to read
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should give a string like '/Users/juliaahlvind'
wdir = os.path.join(home, 'Documents', 'projekt_1', telescope, 'data') # path to mother directory that includes all obsid folders
base = os.path.join(home, 'Documents', 'projekt_1')
start_folder = os.getcwd()

sys.path.insert(0,f'{base}/pipelines')
from pyxspec_functions import *

# create folder for results if it does not exist
#–––––––––––––––––––––––––––––––––––––––––––––––
folders_in_tele = os.listdir(f'{base}/{telescope}')
if folder not in folders_in_tele:
    os.mkdir(f'{base}/{telescope}/{folder}')


# setting up cstatistics and limited iterations
#––––––––––––––––––––––––––––––––––––––––––––––
Fit.nIterations = 100
Fit.statMethod = "cstat"
Fit.query = "yes" # run fitting until this limit:Fit.nIterations is reached
#Xset.chatter = 0 # less printing in the terminal
#Xset.logChatter = 0

# read the input file
# –––––––––––––––––––
obsid_file = ascii.read(f'{base}/{telescope}/test.txt') #main_list_EjGALAX
with open(f'{base}/{telescope}/test.txt') as p:
    lines = p.readlines()
obsid = []
for line in lines:
    obsid.append(line.split('\t')[1])
obsid.pop(0)
    
name_file = obsid_file['Name']
dist = obsid_file['dist(Mpc)']
epoch = obsid_file['epoch']
efft = obsid_file['efft(ks)']
type = obsid_file['type']
RA = obsid_file['ra']
DEC = obsid_file['dec']
dateSN = obsid_file['dateSN']
dateXray = obsid_file['dateXray']

# nH file to cross check with
nH_file = ascii.read(f"{base}/nH_60Mpc.csv", delimiter=';')
nH_all = nH_file['NH_tot_(Mean)_[atoms_cm-2]']
nH_names = nH_file['Name']
string_nH=[]
for k in nH_names:
    string_nH.append(k)
    
Mpc_cm = 3.08567758*1e24  # mult with distance
kev_erg = 1/624150964.7  # mult with integral output keV

# energy intervals to loop over
E_int = [[0.5,10.0],[2.0,10.0],[7.0,10.0]]

#f = open(f'{base}/{telescope}/{folder}/obsflux_GALAX.csv','+w')  # creates a csv document with all outputs
f = open(f'{base}/{telescope}/{folder}/testtesttest.csv','+w')  # creates a csv document with all outputs
f.write('"Name"\t"obsid"\t"dist(Mpc)"\t"type"\t"epoch"\t"efft(ks)"\t"RA"\t"DEC"\t"mean_nH(1e22atoms_cm-2)"\t"dateSN"\t"dateXray"\t"net_count_rate"\t"galax_kT"\t"galax_kT_norm"\t"galax_kT_norm_confInt(1.6σ)"\t"kT"\t"confInt(1.6σ)"\t"norm"\t"confInt(1.6)"\t"cstat"\t"chi2"\t"dof"\t"obsflux_05_10"\t"flux_confInt_05_10(1.6σ)"\t"Lum_05_10"\t"lum_confInt_05_10(1.6σ)"\t"flux_confInt_05_10(3σ)"\t"lum_confInt_05_10(3σ)"\t"obsflux_2_10"\t"flux_confInt_2_10(1.6σ)"\t"Lum_2_10"\t"lum_confInt_2_10(1.6σ)"\t"flux_confInt_2_10(3σ)"\t"lum_confInt_2_10(3σ)"\n')

# Start of the actuall program. Loop over each obsid
# ––––––––––––––––––––––––––––––––––––––––––––––––––
for i in range(len(obsid)):
    
    name = name_file[i]
    
    if telescope=='Chandra':
        # go into each repro file of the corresponding obsid
        list_files_in_repro = os.listdir(f'{wdir}/{obsid[i]}/repro')
        spec_dir = f'{wdir}/{obsid[i]}/repro'
        spec_file_name = f'src_{name}_acis_grp.pi'
        
    elif telescope=='XMM_ny':
        obsid[i]=obsid[i].replace('s','')
        list_files_in_repro = os.listdir(f'{wdir}/{obsid[i]}/pn')
        spec_dir = f'{wdir}/{obsid[i]}/pn'
        spec_file_name = f'{name}_pn_grp.pha'
    
    # match nH value with SN
    if name_file[i].startswith('i'):
        index_nH = np.where(np.array(string_nH)==name_file[i])[0]
    else:
        index_nH = np.where(np.array(string_nH)=='SN'+name_file[i])[0]
    nH = float(np.array(nH_all)[index_nH][0])/1e22
    
    # find the spectral file
    spec_file = f'{spec_dir}/{spec_file_name}'
    os.chdir(spec_dir) # same as cd but for "python environment"
        
    
    if (spec_file_name in list_files_in_repro):
        # clear all previous data
        xspec.AllData.clear()
        
        # load spectrum
        s = Spectrum(spec_file)
        
        # write out info to the list
        f.write(f"{name}\t")
        if telescope=='XMM_ny':
            f.write(f"s{obsid[i]}\t")
        else:
            f.write(f"{obsid[i]}\t")
        f.write(f"{dist[i]}\t")
        type[i] = type[i].replace('/','_')
        f.write(f"{type[i]}\t")
        f.write(f"{epoch[i]}\t")
        if telescope=='SWIFT':
            exp = s.exposure
            f.write(f"{round(exp/1e3,2)}\t")
        else:
            f.write(f"{efft[i]}\t")
        f.write(f"{RA[i]}\t")
        f.write(f"{DEC[i]}\t")
        f.write(f"{nH}\t")
        f.write(f"{dateSN[i]}\t")
        f.write(f"{dateXray[i]}\t")
        
        # restrain to the desired energy interval
        s.notice('all')
        s.ignore("**-0.5 10.0-**")
        AllData.ignore("bad")
        
        AllData.show()
        
        counts = float(s.exposure)*float(s.rate[0])
        net_count_rate = AllData(1).rate[0]
        net_c_s_format = "{:e}".format(net_count_rate)
        f.write(f"{net_c_s_format}\t")
        
        
        # If the spectrum gives zero counts, print nothing
        if counts<1:
            if i==len(obsid)-1:
                f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t ")
            else:
                f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n")
        else:
            
            # clear previous models
            mo = None

            # apply model
            mo = Model("tbabs*(mekal+mekal)")


            # define model parameters
            # –––––––––––––––––––––––
            # tbabs
            mo(1).values = nH
            mo(1).frozen = True
            
            # interaction mekal
            mo(2).values = [0.2, 0.01, 0.1, 0.1, 10.0, 10.0]
            mo(3).values = 1
            mo(4).values = 1
            mo(5).values = 0
            mo(6).values = 1
            mo(7).values = [1e-6, 1e-9, 1e-20, 1e-20, 1e2,1e2]
            
            # galaxy
            galax_mekal_kT, galax_mekal_norm = galax_func(name, obsid[i], telescope)
            mo(8).values = galax_mekal_kT  # set the model parameter
            mo(8).frozen= True
            mo(9).values = 1
            mo(10).values = 1
            mo(11).values = 0
            mo(12).values = 1
            mo(13).values = galax_mekal_norm  # set the model parameter
            
            # Fitting of the model first step without cflux
            # ––––––––––––––––––––––––––––––––––––––––––––––
            try:
                # fit the model to the spectrum
                Fit.perform()
                
                # extract values from fit
                # –––––––––––––––––––––––
                Fit.error("2.706 2") # 1σ error

                kT = mo(2).values[0]
                kT_err = mo(2).error
                kT_low, kT_up = kT_err[0], kT_err[1]
                kT_2= mo(3).values[0]
                kT_3 = mo(4).values[0]
                kT_4 = mo(5).values[0]
                kT_5 = mo(6).values[0]
                
                Fit.error("2.706 7") # 1σ error
                kT_norm = mo(7).values[0]
                kT_norm_err = mo(7).error
                kT_norm_low, kT_norm_up = kT_norm_err[0], kT_norm_err[1]
                    
                # cstatistics info
                cstat = xspec.Fit.statistic  # cstat
                chi2 = xspec.Fit.testStatistic  # chi2 stat.
                dof = xspec.Fit.dof  # degrees of freedom = number of data points-number of free parameters
                
                Fit.error("2.706 13")
                kT_galax = mo(8).values[0]
                kT_norm_galax = mo(13).values[0]
                kT_galax_norm_err = mo(13).error
                kT_galax_norm_low, kT_galax_norm_up = kT_galax_norm_err[0], kT_galax_norm_err[1]
                    
                # write out the results
                f.write(f"{round(kT_galax,4)}\t")
                f.write(f"{format(kT_norm_galax,'.4e')}\t")
                f.write(f"({format(kT_galax_norm_low,'.4e')},{format(kT_galax_norm_up,'.4e')})\t")
                f.write(f"{round(kT,4)}\t")
                f.write(f"({round(kT_low,3)},{round(kT_up,3)})\t")
                f.write(f"{format(kT_norm,'.4e')}\t")
                f.write(f"({format(kT_norm_low,'.4e')},{format(kT_norm_up,'.4e')})\t")
                f.write(f"{round(cstat,2)}\t")
                f.write(f"{round(chi2,2)}\t")
                f.write(f"{dof}\t")
                
                
                # save a parameter that allowes for plotting and saving the plot & model
                try1 = True
                
            except Exception as e:
                print(f"Error in Fit.perform() without cflux : {e}")
                # If the inital fitting does not work, we print blank
                if i==len(obsid)-1:
                    f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t")
                else:
                    f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n")
                    try1 = False
                
            if try1 == True:
            
                # save the models for each line of sight within the repro folder of each obsid
                model_name = f'model_{name}_{obsid[i]}_obsflux.xcm'
                if (os.path.isfile(model_name)):
                     os.remove(model_name)
                Xset.save(f'model_{name}_{obsid[i]}_obsflux.xcm', info='m')
                
                # plott
                Plot.device = '/xs'
                Plot.xAxis = "keV"
                xspec.Plot('ldata','del')
                xspec.Plot.commands = ()
                # save the corresponding plot in the folder fit_plots
                save_data = f'specfit_delchi_{name}_{obsid[i]}_obsflux.qdp' # data for plot
               
                if (os.path.isfile(save_data)):
                    os.remove(save_data)
               
                xspec.Plot.device = '/null'
                xspec.Plot.add = True
                xspec.Plot.addCommand(f'wd {save_data}')
                xspec.Plot("ld","del")
                names = ['e','de','rate','rate_err','total']  # energi(x-axel), step in energy, y-axel data, y-axel error data, totala model
                ncomp = len(mo.componentNames)
                for a in range(ncomp):
                    names.append(f'model{a}')
               
                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsflux.qdp',skiprows=3,names=names, delimiter=' ')
                df_original = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsflux.qdp',skiprows=3,names=names, delimiter=' ')
                index = np.where(df.e=='NO')[0][0]
                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsflux.qdp', skiprows=3, names=names, delimiter=' ',nrows=index)
                
                fig = plt.figure(1)
                ax = fig.add_subplot(211)
               
                # Plot using Matplotlib:
                ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err,fmt='.',label='data')
                ax.step(df.e, df.total, color='black',label='Total model',linewidth=2)
                ax.set_ylabel(r'count$~$s$^{-1}~$keV$^{-1}$')
                ax.set_xscale("linear")
                ax.set_yscale("log")
                ax.legend()

                df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_obsflux.qdp', skiprows=index+1+3, names=names, delimiter=' ')
                ax = fig.add_subplot(2,1,2)
                ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err, fmt='.')
                ax.set_ylabel('(data-model)/error')
                ax.set_xlabel('Energy (keV)')
                ax.axhline(y=0, linestyle='--', linewidth=1, color='black')
                figure_name = f'{base}/{telescope}/fits/obsflux/{name}_{obsid[i]}.png'
                if os.path.isfile(figure_name):
                    os.remove(figure_name)
                fig.savefig(figure_name, format='png')
                #plt.show()
                fig.clf()
                Plot.commands = ()

                # Loop over each energy interval for cstat and do a new fitting
                # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
                for m in range(2):
                
                    # to check coutns in "new range"
                    s.notice('all')
                    s.ignore(f"**-{E_int[m][0]} {E_int[m][1]}-**")
                    counts2 = float(s.exposure)*float(s.rate[0])
                    s.notice('all')
                    s.ignore("**-0.5 10.0-**")
                    AllData.ignore("bad")
                    
                    if counts2<1:
                        
                        if m==1:
                            f.write("\t \t \t \t \t \n")
                        else:
                            f.write("\t \t \t \t \t \t")
                        
                    else:
                    
                        # set up the new model
                        # ––––––––––––––––––––
                        mo = None
                        mo = Model("cflux*(tbabs*(mekal+mekal))" )
                        
                        # cflux
                        mo(1).values = E_int[m][0]
                        mo(1).frozen = True
                        mo(2).values = E_int[m][1]
                        mo(2).frozen = True
                        mo(3).values = -13
                        
                        # tbabs
                        mo(4).values = nH
                        mo(4).frozen = True
                        
                        # interaction
                        mo(5).values = kT
                        mo(5).frozen = True
                        mo(6).values = kT_2
                        mo(7).values = kT_3
                        mo(8).values = kT_4
                        mo(9).values = kT_5
                        mo(10).values = kT_norm
                        mo(10).frozen = True
                        
                        # galaxy
                        mo(11).values = kT_galax  # set the model parameter
                        mo(11).frozen= True
                        mo(12).values = 1
                        mo(13).values = 1
                        mo(14).values = 0
                        mo(15).values = 1
                        mo(16).values = kT_norm_galax  # set the model parameter
                        mo(16).frozen = True
            
                        # Try to fit the new model with cflux for 1σ
                        # –––––––––––––––––––––––––––––––––––––––––
                        try:
                            Fit.perform()
                            mo(5).frozen = False
                            
                            # save flux values
                            Fit.error("2.706 3")
                            cflux = mo(3).values[0]
                            cflux_err = mo(3).error
                            cflux_under, cflux_over = cflux_err[0], cflux_err[1]

                            # if flux is an even number, do steppar instead
                            if (cflux % 1==0) or cflux<-20:
                                raise Exception("The result of Fit.Perform() is not sufficient, do steppar instead")

                            
                        except Exception as e: # fitting for Γ didn't work from the begining so we do steppar instead
                            print(f"Error in Fit.perform(), do steppar instead. Error: {e}")
                            try:
                                cflux, cflux_under, cflux_over = call_steppar(cstat=2.706,paramNr=3,step_start=-16, step_end=-10)
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
                            
                            flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                            
                            # if flux is an even number, do steppar instead
                            if (cflux % 1==0) or cflux<-20:
                                raise Exception("The result is not sufficient")
                        
                        # if fit not possible to steppar directly
                        except Exception as e:
                            print(f"Error in Fit.perform() for 3σ, do steppar instead: {e}")
                            try:
                                cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=3,step_start=-16, step_end=-10)

                                if cflux_under==0 or cflux_under==1:
                                    cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=3,step_start=-16, step_end=-10)
                                    cflux_under = 0
                                flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                            except Exception as e:
                                print(f"Error in steppar for 2.706, cflux=0: {e}")
                                flux, flux_under, flux_over = 0,0,0
                                    
                        # calculate luminosity
                        L, L_min, L_max = lum_funk(10**cflux,10**cflux_under,10**cflux_over,dist[i])
                            
                        # write out the results
                        if flux==1 or flux==0:
                            if m==1:
                                f.write("\t \n")
                            else:
                                f.write("\t \t")
                        else:
                            f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                            if m==1:
                                f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\n")
                            else:
                                f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")

                        # save model with cflux
                        model_name = f'model_{name}_{obsid[i]}_obsflux.xcm'
                        if (os.path.isfile(model_name)):
                            os.remove(model_name)
                        Xset.save(f'model_{name}_{obsid[i]}_obsflux.xcm', info='m')
                
        mo = None
        xspec.AllData.clear()

f.close()
os.chdir(start_folder)


