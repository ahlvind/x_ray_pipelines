
import numpy as np
from astropy.io import ascii
import glob
import os
from scipy.integrate import quad
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy.ma as ma


#9
f_1σ = [   -13.5201      ,   -13.8131   ,  -13.3103  ]
f_3σ = [     -13.5201       , -14.212  ,   -13.1708   ]



f_1σ_med = [  0,0,0 ]
f_3σ_med = [    0,0,0 ]
f_1σ_90 = [   0,0,0  ]
f_3σ_90 = [   0,0,0]

f_1σ_2_10 = [    -13.6369  ,  -13.9299  ,   -13.4271  ]
f_3σ_2_10 = [     -13.6369  ,   -14.3288  ,   -13.2877   ]
f_1σ_2_10_med = [   0,0,0  ]
f_3σ_2_10_med = [   0,0,0 ]
f_1σ_2_10_90 = [     0,0,0 ]
f_3σ_2_10_90 = [   0,0,0 ]


gamma, gamma_med, gamma_90 = 2.0, 2.0 , 2.0
g_l, g_u = 0,0
g_l_med, g_u_med = 0,0
g_l_90, g_u_90 = 0,0

norm, norm_med, norm_90 =  7.61133E-07, 0,0
norm_l, norm_u =   0 , 1.17583e-06 
norm_l_med, norm_u_med =   0 ,     0
norm_l_90, norm_u_90 =   0  ,   0

kT, kT_med, kT_90 =     5.16593 , 0,0
kT_l, kT_u =    1.73075   ,         0
kT_l_med, kT_u_med = 0,0
kT_l_90, kT_u_90 =0,0

kT_norm, kT_norm_med, kT_norm_90 =       2.10783E-05 ,  0,0
kT_norm_l, kT_norm_u =  1.07401e-05 , 3.42031e-05
kT_norm_l_med, kT_norm_u_med =  0,0
kT_norm_l_90, kT_norm_u_90 = 0,0
dist = 4.7


# 5.66122e-06
# (0.0,6.89515e-06)
Mpc_cm = 3.08567758*1e24

# luminosity funktion: derives the luminosities and errors
def lum_funk(flux,err_low_flux,err_up_flux,dist):
    L = flux*4*np.pi*(float(dist)*Mpc_cm)**2
    L_min = err_low_flux*4*np.pi*(float(dist)*Mpc_cm)**2
    L_max = err_up_flux*4*np.pi*(float(dist)*Mpc_cm)**2
    return L, L_min, L_max
    
flux_1σ = [10**f_1σ[0],10**f_1σ[1],10**f_1σ[2]] # best, low, up
flux_1σ_med = [10**f_1σ_med[0],10**f_1σ_med[1],10**f_1σ_med[2]] # best, low, up
flux_1σ_90 = [10**f_1σ_90[0],10**f_1σ_90[1],10**f_1σ_90[2]] # best, low, up

flux_3σ = [10**f_3σ[0],10**f_3σ[1],10**f_3σ[2]] # best, low, up
flux_3σ_med = [10**f_3σ_med[0],10**f_3σ_med[1],10**f_3σ_med[2]] # best, low, up
flux_3σ_90 = [10**f_3σ_90[0],10**f_3σ_90[1],10**f_3σ_90[2]] # best, low, up

flux_1σ_2_10 = [10**f_1σ_2_10[0],10**f_1σ_2_10[1],10**f_1σ_2_10[2]] # best, low, up
flux_1σ_2_10_med = [10**f_1σ_2_10_med[0],10**f_1σ_2_10_med[1],10**f_1σ_2_10_med[2]] # best, low, up
flux_1σ_2_10_90 = [10**f_1σ_2_10_90[0],10**f_1σ_2_10_90[1],10**f_1σ_2_10_90[2]] # best, low, up

flux_3σ_2_10 = [10**f_3σ_2_10[0],10**f_3σ_2_10[1],10**f_3σ_2_10[2]] # best, low, up
flux_3σ_2_10_med = [10**f_3σ_2_10_med[0],10**f_3σ_2_10_med[1],10**f_3σ_2_10_med[2]] # best, low, up
flux_3σ_2_10_90 = [10**f_3σ_2_10_90[0],10**f_3σ_2_10_90[1],10**f_3σ_2_10_90[2]] # best, low, up


Lf_1σ, L_minf_1σ, L_maxf_1σ = lum_funk(flux_1σ[0],flux_1σ[1], flux_1σ[2], dist)
Lf_1σ_med, L_minf_1σ_med, L_maxf_1σ_med = lum_funk(flux_1σ_med[0],flux_1σ_med[1], flux_1σ_med[2], dist)
Lf_1σ_90, L_minf_1σ_90, L_maxf_1σ_90 = lum_funk(flux_1σ_90[0],flux_1σ_90[1], flux_1σ_90[2], dist)

Lf_3σ, L_minf_3σ, L_maxf_3σ = lum_funk(flux_3σ[0],flux_3σ[1], flux_3σ[2], dist)
Lf_3σ_med, L_minf_3σ_med, L_maxf_3σ_med = lum_funk(flux_3σ_med[0],flux_3σ_med[1], flux_3σ_med[2], dist)
Lf_3σ_90, L_minf_3σ_90, L_maxf_3σ_90 = lum_funk(flux_3σ_90[0],flux_3σ_90[1], flux_3σ_90[2], dist)

Lf_1σ_2_10, L_minf_1σ_2_10, L_maxf_1σ_2_10 = lum_funk(flux_1σ_2_10[0],flux_1σ_2_10[1], flux_1σ_2_10[2], dist)
Lf_1σ_2_10_med, L_minf_1σ_2_10_med, L_maxf_1σ_2_10_med = lum_funk(flux_1σ_2_10_med[0],flux_1σ_2_10_med[1], flux_1σ_2_10_med[2], dist)
Lf_1σ_2_10_90, L_minf_1σ_2_10_90, L_maxf_1σ_2_10_90 = lum_funk(flux_1σ_2_10_90[0],flux_1σ_2_10_90[1], flux_1σ_2_10_90[2], dist)

Lf_3σ_2_10, L_minf_3σ_2_10, L_maxf_3σ_2_10 = lum_funk(flux_3σ_2_10[0],flux_3σ_2_10[1], flux_3σ_2_10[2], dist)
Lf_3σ_2_10_med, L_minf_3σ_2_10_med, L_maxf_3σ_2_10_med = lum_funk(flux_3σ_2_10_med[0],flux_3σ_2_10_med[1], flux_3σ_2_10_med[2], dist)
Lf_3σ_2_10_90, L_minf_3σ_2_10_90, L_maxf_3σ_2_10_90 = lum_funk(flux_3σ_2_10_90[0],flux_3σ_2_10_90[1], flux_3σ_2_10_90[2], dist)

f = open('manuell_output.csv', '+w')
f.write('"PhoInd_05_10"\t"confInt(1.6σ)_05_10"\t"kT_norm_05_10"\t"confInt_05_10(1.6σ)"\t"PL_norm_05_10"\t"PL_norm_confint"\t"kT_05_10"\t"confInt(1.6σ)_05_10"\t\t"obsflux_05_10(1.6σ)"\t"confInt_05_10(1.6σ)"\t"Lum_05_10(1.6σ)"\t"confInt_05_10(1.6σ)"\t\t"confInt_05_10(3σ)"\t"confInt_05_10(3σ)"\t\t"obsflux_2_10(1.6σ)"\t"confInt_2_10(1.6σ)"\t"Lum_2_10(1.6σ)"\t"confInt_2_10(1.6σ)"\t\t"confInt_2_10(3σ)"\t"confInt_2_10(3σ)"\n')
f.write(str(kT)+'\t')
f.write(f'({round(kT_l,2)},{round(kT_u,2)})\t')
f.write(format(kT_norm,'.4e')+'\t')
f.write('('+format(kT_norm_l,'.4e')+','+format(kT_norm_u,'.4e')+')\t')
f.write(str(round(gamma,2))+'\t')
f.write(f'({round(g_l,2)},{round(g_u,2)})\t')
f.write(format(norm,'.4e')+'\t')
f.write('('+format(norm_l,'.4e')+','+format(norm_u,'.4e')+')\t\t')

f.write(format(flux_1σ[0],'.4e')+'\t')
f.write('('+format(flux_1σ[1],'.4e')+','+format(flux_1σ[2],'.4e')+')\t')
f.write(format(Lf_1σ,'.4e')+'\t')
f.write('('+format(L_minf_1σ,'.4e')+','+format(L_maxf_1σ,'.4e')+')\t')
f.write('\t')
f.write('('+format(flux_3σ[1],'.4e')+','+format(flux_3σ[2],'.4e')+')\t')
f.write('('+format(L_minf_3σ,'.4e')+','+format(L_maxf_3σ,'.4e')+')\t\t')

f.write(format(flux_1σ_2_10[0],'.4e')+'\t')
f.write('('+format(flux_1σ_2_10[1],'.4e')+','+format(flux_1σ_2_10[2],'.4e')+')\t')
f.write(format(Lf_1σ_2_10,'.4e')+'\t')
f.write('('+format(L_minf_1σ_2_10,'.4e')+','+format(L_maxf_1σ_2_10,'.4e')+')\t')
f.write('\t')
f.write('('+format(flux_3σ_2_10[1],'.4e')+','+format(flux_3σ_2_10[2],'.4e')+')\t')
f.write('('+format(L_minf_3σ_2_10,'.4e')+','+format(L_maxf_3σ_2_10,'.4e')+')\n')


# median
f.write(str(kT_med)+'\t')
f.write(f'({round(kT_l_med,2)},{round(kT_u_med,2)})\t')
f.write(format(kT_norm_med,'.4e')+'\t')
f.write('('+format(kT_norm_l_med,'.4e')+','+format(kT_norm_u_med,'.4e')+')\t')
f.write(str(round(gamma_med,2))+'\t')
f.write(f'({round(g_l_med,2)},{round(g_u_med,2)})\t')
f.write(format(norm_med,'.4e')+'\t')
f.write('('+format(norm_l_med,'.4e')+','+format(norm_u_med,'.4e')+')\t\t')

f.write(format(flux_1σ_med[0],'.4e')+'\t')
f.write('('+format(flux_1σ_med[1],'.4e')+','+format(flux_1σ_med[2],'.4e')+')\t')
f.write(format(Lf_1σ_med,'.4e')+'\t')
f.write('('+format(L_minf_1σ_med,'.4e')+','+format(L_maxf_1σ_med,'.4e')+')\t')
f.write('\t')
f.write('('+format(flux_3σ_med[1],'.4e')+','+format(flux_3σ_med[2],'.4e')+')\t')
f.write('('+format(L_minf_3σ_med,'.4e')+','+format(L_maxf_3σ_med,'.4e')+')\t\t')

f.write(format(flux_1σ_2_10_med[0],'.4e')+'\t')
f.write('('+format(flux_1σ_2_10_med[1],'.4e')+','+format(flux_1σ_2_10_med[2],'.4e')+')\t')
f.write(format(Lf_1σ_2_10_med,'.4e')+'\t')
f.write('('+format(L_minf_1σ_2_10_med,'.4e')+','+format(L_maxf_1σ_2_10_med,'.4e')+')\t')
f.write('\t')
f.write('('+format(flux_3σ_2_10_med[1],'.4e')+','+format(flux_3σ_2_10_med[2],'.4e')+')\t')
f.write('('+format(L_minf_3σ_2_10_med,'.4e')+','+format(L_maxf_3σ_2_10_med,'.4e')+')\n')

# 90th
f.write(str(kT_90)+'\t')
f.write(f'({round(kT_l_90,2)},{round(kT_u_90,2)})\t')
f.write(format(kT_norm_90,'.4e')+'\t')
f.write('('+format(kT_norm_l_90,'.4e')+','+format(kT_norm_u_90,'.4e')+')\t')
f.write(str(round(gamma_90,2))+'\t')
f.write(f'({round(g_l_90,2)},{round(g_u_90,2)})\t')
f.write(format(norm_90,'.4e')+'\t')
f.write('('+format(norm_l_90,'.4e')+','+format(norm_u_90,'.4e')+')\t\t')

f.write(format(flux_1σ_90[0],'.4e')+'\t')
f.write('('+format(flux_1σ_90[1],'.4e')+','+format(flux_1σ_90[2],'.4e')+')\t')
f.write(format(Lf_1σ_90,'.4e')+'\t')
f.write('('+format(L_minf_1σ_90,'.4e')+','+format(L_maxf_1σ_90,'.4e')+')\t')
f.write('\t')
f.write('('+format(flux_3σ_90[1],'.4e')+','+format(flux_3σ_90[2],'.4e')+')\t')
f.write('('+format(L_minf_3σ_90,'.4e')+','+format(L_maxf_3σ_90,'.4e')+')\t\t')

f.write(format(flux_1σ_2_10_90[0],'.4e')+'\t')
f.write('('+format(flux_1σ_2_10_90[1],'.4e')+','+format(flux_1σ_2_10_90[2],'.4e')+')\t')
f.write(format(Lf_1σ_2_10_90,'.4e')+'\t')
f.write('('+format(L_minf_1σ_2_10_90,'.4e')+','+format(L_maxf_1σ_2_10_90,'.4e')+')\t')
f.write('\t')
f.write('('+format(flux_3σ_2_10_90[1],'.4e')+','+format(flux_3σ_2_10_90[2],'.4e')+')\t')
f.write('('+format(L_minf_3σ_2_10_90,'.4e')+','+format(L_maxf_3σ_2_10_90,'.4e')+')\t')


f.close()


def PL(E,Norm,Gamma):
    """ Mannually calculate flux component over powerlaw component """
    return Norm*E**(-Gamma)
    
kev_erg = 1/624150964.7
I, err = quad(PL,0.5,10.0, args = (norm,gamma))
abs_flux = I*kev_erg
print(format(abs_flux,'.4e'))

"""
f_1σ = []
f_3σ = [  ]
f_1σ_med = []
f_3σ_med = [ ]
f_1σ_90 = [  ]
f_3σ_90 = [  ]

f_1σ_2_10 = [  ]
f_3σ_2_10 = [  ]
f_1σ_2_10_med = [  ]
f_3σ_2_10_med = [  ]
f_1σ_2_10_90 = [  ]
f_3σ_2_10_90 = [  ]


gamma, gamma_med, gamma_90 = 2.0, 2.0 , 2.0
g_l, g_u = 0,0
g_l_med, g_u_med = 0,0
g_l_90, g_u_90 = 0,0
norm, norm_med, norm_90 =
norm_l, norm_u =
norm_l_med, norm_u_med =
norm_l_90, norm_u_90 =


kT, kT_med, kT_90 =
kT_l, kT_u =
kT_l_med, kT_u_med =
kT_l_90, kT_u_90 =

kT_norm, kT_norm_med, kT_norm_90 =
kT_norm_l, kT_norm_u =
kT_norm_l_med, kT_norm_u_med =
kT_norm_l_90, kT_norm_u_90 =
"""
