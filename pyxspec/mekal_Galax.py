
import os, sys


telescope = 'XMM_ny' #'Chandra'

# Set up and files to read
#––––––––––––––––––––––––––
home = os.path.expanduser('~') # should giv a string like '/Users/juliaahlvind'
wdir = f'{home}/Documents/projekt_1/{telescope}/data'  # path tp mother directory that includes all obsid folders
base = f'{home}/Documents/projekt_1/'
start_folder = os.getcwd()

sys.path.insert(0,f'{base}/pipelines')
from pyxspec_functions import *

# setting up cstatistics and limited iterations
Fit.nIterations = 100
Fit.statMethod = "cstat"

# Chandra file with SN names and obsid
#obsid_file = ascii.read(f'{base}/{telescope}/main_list_GALAX.txt')
obsid_file = ascii.read(f'{base}/{telescope}/test.txt')

# Chandra file with SN names and obsid
#with open(f'{base}/{telescope}/main_list_GALAX.txt') as p:
with open(f'{base}/{telescope}/test.txt') as p:
    lines = p.readlines()
obsid = []
for line in lines:
    obsid.append(line.split('\t')[1])
obsid.pop(0)
name_file = obsid_file['Name']
epoch = obsid_file['epoch']

# nH file to cross check with
nH_file = ascii.read(f"{base}/nH_60Mpc.csv", delimiter=';')
nH_all = nH_file['NH_tot_(Mean)_[atoms_cm-2]']
nH_names = nH_file['Name']
string_nH=[]
for k in nH_names:
    string_nH.append(k)
    
Mpc_cm = 3.08567758*1e24  # mult with distance
kev_erg = 1/624150964.7  # mult with integral output keV

# Start of the actuall program. Loop over each obsid
# ––––––––––––––––––––––––––––––––––––––––––––––––––
saknar_spec_o = []
saknar_spec_n = []
for i in range(len(obsid)):

    # remove SN in name
    name = name_file[i]
    
    if telescope=='Chandra':
        # go into each repro file of the corresponding obsid
        list_files_in_repro = os.listdir(wdir+'/'+str(obsid[i])+"/repro")
        spec_dir = f'{wdir}/'+str(obsid[i])+'/repro'
        spec_file_name = f'galax_{name}_acis_grp.pi'
        
    elif telescope=='XMM_ny':
        list_files_in_repro = os.listdir(wdir+'/'+str(obsid[i])+"/pn")
        spec_dir = f'{wdir}/'+str(obsid[i])+'/pn'
        spec_file_name = f'{name}_pn_galax_grp.pha'

    # match nH value with SN
    if name_file[i].startswith('i'):
        index_nH = np.where(np.array(string_nH)==name_file[i])[0]
    else:
        index_nH = np.where(np.array(string_nH)=='SN'+name_file[i])[0]
    nH = float(np.array(nH_all)[index_nH][0])/1e22

    # find the spectral file
    spec_file = f'{spec_dir}/{spec_file_name}'
    
    
    current = os.getcwd()  # current folder we are in
    os.chdir(spec_dir) # same as cd but for "python environment"

    # if spectral file does exist
    if (spec_file_name in list_files_in_repro):
        print('test')

        # load spectrum
        s = Spectrum(spec_file)
        
        # restrain to the desired energy interval
        s.notice('all')
        s.ignore("**-0.5 10.0-**")
        s.ignore("bad")
        
        # apply model
        mo = Model("tbabs*(mekal)")
        
        # define model parameters
        # –––––––––––––––––––––––
        # tbabs
        mo(1).values = nH
        mo(1).frozen = True
        
        # mekal
        mo(2).values = [0.1, 0.01, 0.1, 0.1, 2.0, 2.0] #0.1
        mo(3).values = 1
        mo(4).values = 1
        mo(5).values = 0
        mo(6).values = 1
        mo(7).values = 1e-4

        # here we start the fit of the model for 90% confidence
        # –––––––––––––––––––––––––––––––––––––––––––––––––––––
        #try:
        # fit the model to the spectrum
        Fit.query = "yes"
        Fit.perform()
        Fit.error("2.706 2") # value not used, but to ensure that the model finds the global minimum and not local

        Fit.perform()
        # save the models for each line of sight within the repro folder of each obsid
        model_name = 'mekal_'+name+'.xcm'
        
        if (os.path.isfile(model_name)):
            os.remove(model_name)
        Xset.save('mekal_'+name+'.xcm', info='m')
        print('model should have been saved')

        
        # plott
        Plot.device = '/xs'
        Plot.xAxis = "keV"
        xspec.Plot('ldata','del')
        xspec.Plot.commands = ()
        # save the corresponding plot in the folder fit_plots
        save_data = f'specfit_delchi_{name}_{obsid[i]}_mekal_galax.qdp' # data for plot
        
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
            
        df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_mekal_galax.qdp',skiprows=3,names=names, delimiter=' ')
        index = np.where(df.e=='NO')[0][0]
        df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_mekal_galax.qdp', skiprows=3, names=names, delimiter=' ',nrows=index)

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
        

        df = pd.read_table(f'specfit_delchi_{name}_{obsid[i]}_mekal_galax.qdp', skiprows=index+1+3, names=names, delimiter=' ')
        ax = fig.add_subplot(2,1,2)
        ax.errorbar(df.e, df.rate, xerr=df.de ,yerr=df.rate_err, fmt='.')
        ax.set_ylabel('(data-model)/error')
        ax.set_xlabel('Energy (keV)')
        ax.axhline(y=0, linestyle='--', linewidth=1, color='black')
        #plt.show()

        fig.savefig(f'{base}/{telescope}/fits/mekalGalax/{name}_{obsid[i]}.png')
        fig.clf()
        Plot.commands = ()
        
        
        f = open(spec_dir+'/mekal_'+str(name)+'.txt','+w')  # creates a txt document with all outputs
        f.write(str(mo(2).values[0])+'\n')
        f.write(str(mo(3).values[0])+'\n')
        f.write(str(mo(4).values[0])+'\n')
        f.write(str(mo(5).values[0])+'\n')
        f.write(str(mo(6).values[0])+'\n')
        f.write(str(mo(7).values[0])+'\n')
        f.close()
        
        mo = None
        xspec.AllData.clear()
        xspec.AllModels.clear()


    else:
        saknar_spec_o.append(obsid[i])
        saknar_spec_n.append(name)
    os.chdir(wdir)


os.chdir(start_folder)








