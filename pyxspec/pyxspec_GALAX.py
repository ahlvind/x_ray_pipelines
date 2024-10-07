
import os, sys
import numpy.ma as ma


# INPUT PARAMETERS
#telescope = 'XMM_ny' #'Chandra' or 'XMM_ny'
telescope = 'XMM_ny'

folder = '22feb'
# Set up and files to read
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should giv a string like '/Users/juliaahlvind'
wdir = os.path.join(home, 'Documents', 'projekt_1', telescope, 'data') # path to mother directory that includes all obsid folders
base = os.path.join(home, 'Documents', 'projekt_1')
start_folder = os.getcwd()

sys.path.insert(0,f'{base}/pipelines')
from pyxspec_functions import *

# setting up cstatistics and limited iterations
Fit.nIterations = 100
Fit.statMethod = "cstat"
Fit.query = "yes"
Xset.chatter = 0 # less printing in the terminal
Xset.logChatter = 0

# Chandra file with SN names and obsid
obsid_file = read_csv(f'{base}/{telescope}/{folder}/testtesttest.csv')#detected_obsflux_GALAX.csv')

obsid = obsid_file['obsid']
name_file = obsid_file['Name']
dist = obsid_file['dist(Mpc)']
type = obsid_file['type']
epoch = obsid_file['epoch']
efft = obsid_file['efft(ks)']
dateSN = obsid_file['dateSN']
dateXray = obsid_file['dateXray']
RA = obsid_file['RA']
DEC = obsid_file['DEC']

# nH file to cross check with
nH_file = ascii.read(f"{base}/nH_60Mpc.csv", delimiter=';')
nH_all = nH_file['NH_tot_(Mean)_[atoms_cm-2]']
nH_names = nH_file['Name']
string_nH=[]
for k in nH_names:
    string_nH.append(k)
    
Mpc_cm = 3.08567758*1e24  # mult with distance
kev_erg = 1/624150964.7  # mult with integral output keV

line_of_sight_string = ["10th", "median", "90th"]

# energy intervals to loop over
E_int = [[0.5,10.0],[2.0,10.0],[7.0,10.0]]

f = open(f'{base}/{telescope}/{folder}/resultat/testtesttest_tbvarabs.csv','+w')  #detected_pulsar_GALAX creates a csv document with all outputs
f.write('"Name"\t"obsid"\t"dist[Mpc]"\t"type"\t"epoch"\t"efft(ks)"\t"RA"\t"DEC"\t"mean_nH(1e22atoms_cm-2)"\t"dateSN"\t"dateXray"\t"line_of_sight"\t"net_count_rate"\t"galax_kT"\t"galax_kT_norm"\t"galax_kT_norm_confInt(1.6σ)"\t"PhoIndex"\t"confInt(1.6σ)"\t"norm"\t"confInt(1.6σ)"\t"kT"\t"confInt(1.6σ)"\t"norm"\t"confInt(1.6)"\t"cstat"\t"chi2"\t"dof"\t"intflux05_10"\t"PhoInd_05_10"\t"confInt(1.6σ)_05_10"\t"PL_norm_05_10"\t"kT_05_10"\t"confInt(1.6σ)_05_10"\t"kT_norm_05_10"\t"confInt_05_10(1.6σ)"\t"kT_galax_norm_05_10"\t"kT_galax_norm_05_10_confInt(1.6σ)"\t"absflux_05_10"\t"flux_confInt_05_10(1.6σ)"\t"Lum_05_10"\t"lum_confInt_05_10(1.6σ)"\t"flux_confInt_05_10(3σ)"\t"lum_confInt_05_10(3σ)"\t"PhoInd_2_10"\t"confInt(1.6σ)_2_10"\t"PL_norm_2_10"\t"kT_2_10"\t"confInt(1.6σ)_2_10"\t"kT_norm_2_10"\t"confInt_2_10(1.6σ)"\t"kT_galax_norm_2_10"\t"kT_galax_norm_2_10_confInt(1.6σ)"\t"absflux_2_10"\t"flux_confInt_2_10(1.6σ)"\t"Lum_2_10"\t"lum_confInt_2_10(1.6σ)"\t"flux_confInt_2_10(3σ)"\t"lum_confInt_2_10(3σ)"\n')
# Start of the actuall program. Loop over each obsid
# ––––––––––––––––––––––––––––––––––––––––––––––––––
for i in range(len(obsid)):
    
    name = name_file[i]
    name = name.replace(' ','')
    
    if telescope=='XMM_ny':
        obsid[i]=obsid[i].replace('s','')
        
        
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
        index_nH = np.where(np.array(string_nH)==name)[0]
    else:
        index_nH = np.where(np.array(string_nH)=='SN'+name)[0]
    nH = float(np.array(nH_all)[index_nH][0])/1e22

    # find the spectral file
    spec_file = f'{spec_dir}/{spec_file_name}'
    current = os.getcwd()  # current folder we are in
    os.chdir(spec_dir) # same as cd but for "python environment"

    # if spectral file does exist
    if (spec_file_name in list_files_in_repro):
    
        # write out on the list
        f.write(f"{name}\t")
        if telescope=='XMM_ny':
            f.write(f"s{obsid[i]}\t")
        else:
            f.write(f"{obsid[i]}\t")
        f.write(f"{dist[i]}\t")
        type[i] = type[i].replace('/','_')
        f.write(f"{type[i]}\t")
        f.write(f"{epoch[i]}\t")
        f.write(f"{efft[i]}\t")
        f.write(f"{RA[i]}\t")
        f.write(f"{DEC[i]}\t")
        f.write(f"{nH}\t")
        f.write(f"{dateSN[i]}\t")
        f.write(f"{dateXray[i]}\t")
        
        # load spectrum
        s = Spectrum(spec_file)
        
        # restrain to the desired energy interval
        s.notice('all')
        s.ignore("**-0.5 10.0-**")
        AllData.ignore('bad')
        
        # tbvarabs
        tbvarabs_all = tbvarabs_funk(name,obsid[i],type[i],epoch[i])  # = [ten_percentile, median, ninty_percentile]
        
        # loops over all line of sights
        # –––––––––––––––––––––––––––––
        for j in range(3):  # 3 if all 3 line of sights
            if j==0:
                f.write(f"{line_of_sight_string[j]}\t")
            else:
                f.write(f"\t \t \t \t \t \t \t \t \t \t \t {line_of_sight_string[j]}\t")
            
            counts = float(s.exposure)*float(s.rate[0])
            net_count_rate = AllData(1).rate[0]
            net_c_s_format = "{:e}".format(net_count_rate)
            f.write(f"{net_c_s_format} \t")
            
            if ma.is_masked(obsid_file['obsflux_2_10'].data[i])  or obsid_file['obsflux_2_10'].data[i]==' ':
                if i==len(obsid)-1 and j==2:
                    f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t ")
                else:
                    f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n")
                
            else:
                # clear previous model
                mo = None
                
                # apply model
                mo = Model("tbabs*(tbvarabs*po+mekal+mekal)") # galax and interaction mekal
                
                # interaction mekal
                mo(52).values = [0.6, 0.01, 0.1, 0.1, 10.0, 10.0]
                mo(53).values = 1
                mo(54).values = 1
                mo(55).values = 0
                mo(56).values = 1
                mo(57).values = 1e-4

                # define model parameters
                # –––––––––––––––––––––––
                # tbabs
                mo(1).values = nH
                mo(1).frozen = True
                
                # powerlaw
                mo(44).values = [2.0, 0.01, 0.5, 0.5, 3.0, 3.0]   # value, fit delta, min, bot, top, max
                mo(45).values = [1e-6, 1e-9, 1e-20, 1e-20, 1e2,1e2]
                 
                tbvarabs = tbvarabs_all[j]
                # loops over all inputs for tbvarabs
                k=0
                for k in range(len(tbvarabs)):
                    s1, s2, s3, s4, s5, s6 = tbvarabs[k].split()  # each row of tbvarabs aka each element
                    tbvarabs_line = [float(s1), float(s2), float(s3), float(s4), float(s5), float(s6)]
                    mo(k+2).values = tbvarabs_line  # set the model parameter
                    mo(k+2).frozen = True # freezes them

                
                # galaxy
                galax_mekal_kT, galax_mekal_norm = galax_func(name, obsid[i], telescope)
                mo(46).values = galax_mekal_kT  # set the model parameter
                mo(46).frozen= True
                mo(47).values = 1
                mo(48).values = 1
                mo(49).values = 0
                mo(50).values = 1
                mo(51).values = galax_mekal_norm  # set the model parameter
                   

                # Fitting of the model first step without cflux
                # ––––––––––––––––––––––––––––––––––––––––––––––
                try:
                    # fit the model to the spectrum
                    Fit.perform()
                    gamma = mo(44).values[0]
                    Fit.error("2.706 44")  # 1.6σ on parameter 2
                    gamma_err = mo(44).error
                    gamma_low, gamma_up = gamma_err[0], gamma_err[1]
                    
                    # constrains so that photon index is not skyhigh or deepsea low
                    if mo(44).values[0]<0.51 or mo(44).values[0]>2.99 or gamma_up==0:
                        mo(44).values = 2
                        mo(44).values[0] = 2
                        mo(44).frozen = True
                        gamma_low, gamma_up = 0,0
                        Fit.perform()

                    gamma = mo(44).values[0]
                    
                    Fit.error("2.706 52")
                    kT = mo(52).values[0]
                    kT_err = mo(52).error
                    kT_low, kT_up = kT_err[0], kT_err[1]
                                        
                    # constrains so that kT is deepsea low
                    if kT<=0.10:
                        mo(52).values = 0.1
                        mo(52).frozen = True
                        kT_low, kT_up = 0,0
                        
                    kT_2= mo(53).values[0]
                    kT_3 = mo(54).values[0]
                    kT_4 = mo(55).values[0]
                    kT_5 = mo(56).values[0]
                    
                    Fit.error("2.706 57")
                    kT_norm = mo(57).values[0]
                    kT_norm_err = mo(57).error
                    kT_norm_low, kT_norm_up = kT_norm_err[0], kT_norm_err[1]

                    # extract values from fit
                    # –––––––––––––––––––––––
                    # normalinsation of powerlaw
                    Fit.error("2.706 45")  #cstat 2.706=1σ doubble limit on the second parameter
                    norm = mo(45).values[0]
                    norm_err = mo(45).error
                    gamma_norm_low, gamma_norm_up = norm_err[0], norm_err[1]

                    # mekal info
                    # mekal galaxy
                    Fit.error("2.706 51")
                    kT_galax = mo(46).values[0]
                    kT_2_galax= mo(47).values[0]
                    kT_3_galax = mo(48).values[0]
                    kT_4_galax = mo(49).values[0]
                    kT_5_galax = mo(50).values[0]
                    kT_norm_galax = mo(51).values[0]
                    kT_galax_norm_err = mo(51).error
                    kT_galax_norm_low, kT_galax_norm_up = kT_galax_norm_err[0], kT_galax_norm_err[1]
                    
                    # cstatistics info
                    cstat = xspec.Fit.statistic  # cstat
                    chi2 = xspec.Fit.testStatistic  # chi2 stat.
                    dof = xspec.Fit.dof  # degrees of freedom = number of data points-number of free parameters
                            
                    # calculate the absorbed flux
                    I, err = quad(PL,0.5,10.0, args = (norm,gamma))
                    abs_flux = I*kev_erg

                    # write out the results
                    f.write(f"{round(kT_galax,4)}\t")
                    f.write(f"{format(kT_norm_galax,'.4e')}\t")
                    f.write(f"({format(kT_galax_norm_low,'.4e')},{format(kT_galax_norm_up,'.4e')})\t")
                    f.write(f"{round(gamma,4)}\t")
                    f.write(f"({round(gamma_low,3)},{round(gamma_up,3)})\t")
                    f.write(f"{format(norm,'.4e')}\t")
                    f.write(f"({format(gamma_norm_low,'.4e')},{format(gamma_norm_up,'.4e')})\t")
                    f.write(f"{round(kT,4)}\t")
                    f.write(f"({round(kT_low,3)},{round(kT_up,3)})\t")
                    f.write(f"{format(kT_norm,'.4e')}\t")
                    f.write(f"({format(kT_norm_low,'.4e')},{format(kT_norm_up,'.4e')})\t")
                    f.write(f"{round(cstat,2)}\t")
                    f.write(f"{round(chi2,2)}\t")
                    f.write(f"{dof}\t")
                    f.write(f"{abs_flux}\t")
    
                    try1 = True
                except Exception as e:
                    print(f"Error in Fit.perform() without cflux, print empty row. Error: {e}")
                    if i==len(obsid)-1 and j==2:
                        f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t ")
                    else:
                        f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \n")
                        try1 = False
                
                if try1 == True:

                    # save the models for each line of sight within the repro folder of each obsid
                    # Galax
                    Plot.device = '/xs'
                    Plot.xAxis = "keV"
                    #figure = Plot.addCommand('ldata','del')
                    xspec.Plot('ldata','del')
                    # save the corresponding plot in the folder fit_plots
                    save_data = f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp' # data for plot
                    
                    if (os.path.isfile(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp')):
                        os.remove(save_data)
                    
                    xspec.Plot.device = '/null'
                    xspec.Plot.add = True
                    xspec.Plot.addCommand(f'wd {save_data}')
                    xspec.Plot("ld","del")
                    names = ['e','de','rate','rate_err','total']  # energi(x-axel), step in energy, y-axel data, y-axel error data, totala model
                    ncomp = len(mo.componentNames)
                    for a in range(ncomp):
                        names.append(f'model{a}')
                        
                    df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp',skiprows=3,names=names, delimiter=' ')
                    index = np.where(df.e=='NO')[0][0]
                    df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp', skiprows=3, names=names, delimiter=' ',nrows=index)

                    fig = plt.figure(1)
                    ax = fig.add_subplot(211)
                    print(df.e)
                    
                    # Plot using Matplotlib:
                    ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err,fmt='.',label='data')
                    for b in range(ncomp):
                        ax.step(df.e, df[f'model{b}'],label=f'{mo.componentNames[b]}',linewidth=1)
                    ax.step(df.e, df.total, color='black',label='Total model',linewidth=2)
                    ax.set_ylabel(r'count$~$s$^{-1}~$keV$^{-1}$')
                    ax.set_xscale("linear")
                    ax.set_yscale("log")
                    ax.legend()

                    df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_{line_of_sight_string[j]}.qdp', skiprows=index+1+3, names=names, delimiter=' ')
                    ax = fig.add_subplot(2,1,2)
                    ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err, fmt='.')
                    ax.set_ylabel('(data-model)/error')
                    ax.set_xlabel('Energy (keV)')
                    ax.axhline(y=0, linestyle='--', linewidth=1, color='black')
                    #plt.show()

                    fig.savefig(f'{base}/{telescope}/fits/tbvarabs/fit{name}_{obsid[i]}_{line_of_sight_string[j]}_galax')
                    fig.clf()
                    Plot.commands = ()


                    # Loop over each energy interval for cstat and do a new fitting
                    # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
                    for m in range(2):  # Loop over each energy interval for cstat
                    
                        # to check coutns in "new range"
                        s.notice('all')
                        s.ignore(f"**-{E_int[m][0]} {E_int[m][1]}-**")
                        counts2 = float(s.exposure)*float(s.rate[0])
                        s.notice('all')
                        s.ignore("**-0.5 10.0-**")
                        AllData.ignore("bad")
                        
                        if counts2<1:
                            if m==1:
                                f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \n")
                            else:
                                f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t")
                            
                        else:
                            
                            mo = None
                            mo = Model("tbabs*(tbvarabs*cflux*powerlaw + mekal + mekal)" )
                            
                            # Powerlaw
                            mo(47).values = [gamma, 0.01, 0.5, 0.5, 3.0, 3.0]
                            mo(47).frozen = True
                            mo(48).values = norm
                            mo(48).frozen = True

                            # tbabs
                            mo(1).values = nH
                            mo(1).frozen = True
                            
                            # cflux
                            mo(44).values = E_int[m][0]
                            mo(44).frozen = True
                            mo(45).values = E_int[m][1]
                            mo(45).frozen = True
                            mo(46).values = np.log10(abs_flux)
                            
                            # tbvarabs
                            for k in range(len(tbvarabs)):
                                s1, s2, s3, s4, s5, s6 = tbvarabs[k].split()  # each row of tbvarabs aka each element
                                tbvarabs_line = [float(s1), float(s2), float(s3), float(s4), float(s5), float(s6)]
                                mo(k+2).values = tbvarabs_line  # set the model parameter
                                mo(k+2).frozen = True # freezes them

                            # mekal
                            # galaxy
                            mo(49).values = kT_galax
                            mo(49).frozen = True
                            mo(50).values = kT_2_galax
                            mo(51).values = kT_3_galax
                            mo(52).values = kT_4_galax
                            mo(53).values = kT_5_galax
                            mo(54).values = kT_norm_galax
                            mo(54).frozen = True
                            
                            # interaction
                            mo(55).values = [kT, 0.01, 0.1, 0.1, 10.0, 10.0]
                            mo(55).frozen = True
                            mo(56).values = kT_2
                            mo(57).values = kT_3
                            mo(58).values = kT_4
                            mo(59).values = kT_5
                            mo(60).values = kT_norm
                            mo(60).frozen = True
                            
                            # first fit with all paran except cflux frozen
                            Fit.perform()
                            if mo(47).values[0]!=2.0:
                                mo(47).frozen = False
                                
                            mo(54).frozen = False
                            mo(60).frozen = False
                            
                            if mo(55).values[0]>0.1:
                                mo(55).frozen = False
                            
                            
                            try:
                                
                                Fit.perform()
                                if mo(47).values[0]!=2.0:
                                    Fit.error("2.706 47")
                                    gamma_err2 = mo(47).error
                                    gamma_low2 , gamma_up2 = gamma_err[0], gamma_err[1]
                                else:
                                    gamma_low2 , gamma_up2 = 0,0
                                
                                # if PL index is out of range:
                                if mo(47).values[0]<=0.51 or mo(47).values[0]>2.99 or gamma_up2==0:
                                    mo(47).values = 2
                                    mo(47).frozen = True
                                    
                                    Fit.perform()
                                    gamma2 = mo(47).values[0]
                                    gamma_low2, gamma_up2 = 0,0
                                    
                                gamma_norm2 = mo(48).values[0]
                                                                    
                                # constrains so that kT is deepsea low
                                if mo(55).frozen == False:
                                    if mo(55).values[0]<=0.10:
                                        mo(55).values = 0.1
                                        mo(55).frozen = True
                                        Fit.perform()
                                        
                                # save flux values
                                Fit.error("2.706 46")
                                cflux = mo(46).values[0]
                                cflux_err = mo(46).error
                                cflux_under, cflux_over = cflux_err[0], cflux_err[1]

                                # if flux is an even number, do steppar instead
                                if (cflux % 1==0) or cflux<-20:
                                    raise Exception("The result is not sufficient, do steppar instead")
                                            
                                try2 = True
                                
                            except Exception as e:
                                print(f"Error in Fit.perform(), do steppar instead: {e}")
                                mo(47).values = 2
                                mo(47).frozen = True
                                up_lim2, low_lim2 = 0.0, 0.0
                                                                    
                                # constrains so that kT is deepsea low
                                if mo(55).frozen == False:
                                    if mo(55).values[0]<=0.10:
                                        mo(55).values = 0.1
                                        mo(55).frozen = True
                                        
                                try: # make sure PL index is frozen at 2 before steppar, so make a default failed fit again
                                    Fit.perform()
                                except:
                                    try:
                                        cflux, cflux_under, cflux_over = call_steppar(cstat=2.706, paramNr=46,step_start=-16, step_end=-10)
                                    except Exception as e:
                                        print(f"Error steppar for 2.706: {e}. Printing cflux=0 ")
                                        cflux, cflux_under, cflux_over = 0,0,0
                                        crach_par = 0
                                        try2 = False

                            gamma_norm2 = mo(48).values[0]
                            gamma2 = mo(47).values[0]
                            flux, flux_under, flux_over = 10**cflux, 10**cflux_under, 10**cflux_over
                            # calculate luminosity
                            L, L_min, L_max = lum_funk(flux,flux_under,flux_over,dist[i])
                            
                            if mo(55).frozen==False:
                                Fit.error("2.706 55") # 1σ error
                                kT2 = mo(55).values[0]
                                kT_err2 = mo(55).error
                                kT2_low2, kT2_up2 = kT_err2[0], kT_err2[1]
                            else:
                                kT2_low2, kT2_up2 = 0,0
                            
                            Fit.error("2.706 60") # 1σ error
                            kT_norm2 = mo(60).values[0]
                            kT_norm_err2 = mo(60).error
                            kT_norm_low2, kT_norm_up2 = kT_norm_err2[0], kT_norm_err2[1]

                            # normalisation for glaxy mekal
                            Fit.error("2.706 54")
                            kT_norm_galax2 = mo(54).values[0]
                            kT_err2 = mo(54).error
                            kT_norm_galax2_low, kT_norm_galax2_up = kT_err2[0], kT_err2[1]
                            
                            if flux==1 or flux==0:
                                f.write("\t \t \t \t \t \t \t \t \t \t \t \t \t")
                            else:
                                f.write(f"{round(gamma2,4)}\t")
                                f.write(f"({round(gamma_low2,3)},{round(gamma_up2,3)})\t")
                                f.write(f"{format(gamma_norm2,'.4e')}\t")
                                f.write(f"{round(kT2,4)}\t")
                                f.write(f"({round(kT2_low2,3)},{round(kT2_up2,3)})\t")
                                f.write(f"{format(kT_norm2,'.4e')}\t")
                                f.write(f"({format(kT_norm_low2,'.4e')},{format(kT_norm_up2,'.4e')})\t")
                                f.write(f"{format(kT_norm_galax2,'.4e')}\t")
                                f.write(f"({round(kT_norm_galax2_low,3)},{round(kT_norm_galax2_up,3)})\t")
                                f.write(f"{format(flux,'.4e')}\t")
                                f.write(f"({format(flux_under,'.4e')},{format(flux_over,'.4e')})\t")
                                f.write(f"{format(L,'.4e')}\t")
                                f.write(f"({format(L_min,'.4e')},{format(L_max,'.4e')})\t")
                            
                        
                            # Try fitting for cflux with 3σ
                            if try2 == False: # not worth trying fitting for 3σ if 1σ didn't work
                                if m==1:
                                    f.write("\t \n")
                                else:
                                    f.write("\t \t")
                            else:
                                try:
                                    Fit.error("9.0 46")
                                    cflux = mo(46).values[0]
                                    cflux_err = mo(46).error
                                    cflux_under, cflux_over = cflux_err[0], cflux_err[1]
                                    

                                    if cflux_under == 0 or cflux_under == 1:
                                        Fit.error("7.74 46")
                                        cflux = mo(46).values[0]
                                        cflux_err = mo(46).error
                                        cflux_under, cflux_over = 0, cflux_err[1]
                                        
                                    flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                                    # if flux is an even number, do steppar instead
                                    if (cflux % 1==0) or cflux<-20 or cflux_under==0 or cflux_under == 1:
                                        raise Exception("Fit.perform did not yield in valid result, doing steppar instead")
                                    
                                    
                                # if fit not possible to steppar
                                except Exception as e:
                                    print(f"Fitting failed, doing steppar instead. Error: {e}")
                                    try:
                                        cflux, cflux_under, cflux_over = call_steppar(cstat=9.0,paramNr=46,step_start=-16, step_end=-10) # first try fitting for 3σ with upper and lower bound

                                        if cflux_under==0 or cflux_under==1:
                                            cflux, cflux_under, cflux_over = call_steppar(cstat=7.74,paramNr=46,step_start=-16, step_end=-10) # if no lower bound found, new upper bound criteria
                                            cflux_under = 0
                                        flux, flux_under, flux_over = 10**cflux,10**cflux_under,10**cflux_over
                                    except Exception as e:
                                        print(f"Steppar 3σ failed, error: {e}")
                                        flux, flux_under, flux_over = 0,0,0
                                    
                                # calculate luminosity
                                L, L_min, L_max = lum_funk(flux,flux_under,flux_over,dist[i])
                                    
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

                    model_name = 'model_'+name+'_'+str(obsid[i])+'_'+line_of_sight_string[j]+'.xcm'
                    if (os.path.isfile(model_name)):
                        os.remove(model_name)
                    Xset.save('model_'+name+'_'+str(obsid[i])+'_'+line_of_sight_string[j], info='m')
        mo = None
        xspec.AllData.clear()

f.close()
os.chdir(start_folder)
