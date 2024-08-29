import xspec
from xspec import *
import numpy as np
from astropy.io import ascii
import glob
import os
from scipy.integrate import quad
import re
import matplotlib.pyplot as plt
import pandas as pd


def read_csv(fileName):
    """ Reading csv files without worrying about the format """
    file = ascii.read(fileName, delimiter=';')
    if '\t' in file[0][0]:
        file = ascii.read(fileName, delimiter='\t')
        
    return file
    

def find_closest_indices(array, target_value):
    """ This function identifies the index of the position of the value in an array that is closest to the target value """
    min_index, array_zero_min = min(enumerate(array), key=lambda x: x[1] if 0 < x[0] < len(array) - 1 else float('inf'))
   
    min_index = np.argmin(array)
    
    if min_index == 0 or min_index == len(array) - 1:
        target_indices = np.where(array == target_value)[0]
        if target_indices.size > 0:
            min_index = target_indices[0], array_zero_min = array[target_indices[0]]
    
    array_zero_min = array[min_index]
    
    k = min_index
    
    
    # finds the upper bound
    # –––––––––––––––––––––
    while array[k] < target_value:
        k+=1
        if array[k] and array[k]>target_value:
            break
        
    if abs(array[k]-target_value) < abs(array[k-1]-target_value):
        index_up = k
    else:
        index_up = k-1
    
    # also save the nearest upper bound that is still positive
    new_upper_index = index_up+1
    print('index, övre nya flux gräns',new_upper_index, array[new_upper_index])
    
    k = min_index
    # finds the lower bound
    # –––––––––––––––––––––
    # if a lower limits does not exist (can be that we start at negative cstat for best fit)

    if any(value > 0 for value in array[0:k])==True: # no lower limit
        while array[k] < target_value:
            k -= 1
            if array[k] and array[k]>target_value:
                break

        
        if k==min_index: # do nothing, there is no lower limit
            index_low = 0
        
        else:
            if abs(array[k]-target_value) < abs(array[k+1]-target_value):
                index_low = k
            else:
                index_low = k+1
        
        new_lower_index = k-1
    else:
        index_low = 0
        
        new_lower_index = index_up-5
    
    
    print(f'Target value: {target_value}')
    print(f'Original lower bound index: {index_low}, New lower bound index: {new_lower_index}, Value: {array[new_lower_index]}')
    print(f'Original upper bound index: {index_up}, New upper bound index: {new_upper_index}, Value: {array[new_upper_index]}')
    
    return min_index, index_up, index_low, new_lower_index, new_upper_index




def steppar_func(step_start, step_end, tolerance, cstat, steps, paramNr, max_call_limit):
    ## TODO: ***XSPEC Error:  No variable parameters for fit ***Warning: Error search cannot find bracket around delta fit-stat, will automatically exit from trials loop. parameter upper  bound is INVALID.
    ## TODO: xspec exit loop without try, error/warning not cauth in except
    #https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/extended.html
    """ A complete function that makes steppar for a given parameter number  """
    
    # Max number of itterations of steppar. This narrows down the interval of the steppar limits
    if max_call_limit<=0:
        raise ValueError(f"Exceeded the maximum number of function calls ({max_call_limit}).")
        return None

    try:
        Fit.steppar(f"{paramNr} {step_start} {step_end} {steps}")  # steppar through cflux
    except Exception as e:
        print(f"Error in Fit.steppar: {e}")
        
    power = Fit.stepparResults(paramNr)  # log flux output array
    power = np.array(power)
    delcstat = Fit.stepparResults('delstat') # delcstat output array
    delcstat = np.array(delcstat)
  

    min_index, index_up, index_low, new_lower_index, new_upper_index = find_closest_indices(delcstat,cstat)  # index of array where best fit value
    best_power = power[min_index] # best fit value of cflux
    upper_power = power[index_up]

    if index_low==0:
        lower_power = 0
    else:
        lower_power = power[index_low]
    
    step_start, step_end = power[new_lower_index], power[new_upper_index]
    
    # check tolerance if it is close enough to cstat
    if abs(delcstat[index_up]-cstat)>tolerance:
        print('här igen')
        # Do steppar again
        return steppar_func(step_start, step_end, tolerance, cstat, 100, paramNr, max_call_limit-1)
        
    else:
        if lower_power==0:
            best_power = upper_power
        # Return the requested values
        return best_power, lower_power, upper_power

def call_steppar(cstat,paramNr, step_start=-18, step_end=-8, **kwargs):
    """ Collecting function of steppar for an easier call of steppar """
    print('step_start_is',step_start)
    tolerance = 0.1
    steps = 100
    cflux, cflux_under, cflux_over = steppar_func(step_start, step_end, tolerance, cstat, steps,paramNr,max_call_limit=20)
    return cflux, cflux_under, cflux_over
def PL(E,Norm,Gamma):
    """ Mannually calculate flux component over powerlaw component """
    return Norm*E**(-Gamma)
    
def lum_funk(flux,err_low_flux,err_up_flux,dist):
    """ Calculates luminosities from given flux values and distance """
    Mpc_cm = 3.08567758*1e24  # mult with distance

    L = flux*4*np.pi*(float(dist)*Mpc_cm)**2
    L_min = err_low_flux*4*np.pi*(float(dist)*Mpc_cm)**2
    L_max = err_up_flux*4*np.pi*(float(dist)*Mpc_cm)**2
    return L, L_min, L_max
    
def galax_func(name,obsid,telescope):
    """ Function to extract the mekal component values from previous galaxy fitting """
    if telescope=='Chandra':
        with open("/Users/juliaahlvind/Documents/projekt_1/"+telescope+"/data/"+str(obsid)+"/repro/mekal_"+name+".xcm") as f:
            contents = f.readlines()
    else: #telescope=='XMM_ny':
        with open("/Users/juliaahlvind/Documents/projekt_1/XMM_ny/data/"+str(obsid)+"/pn/mekal_"+name+".xcm") as f:
            contents = f.readlines()
    i = 0
    for c in contents:
        if c.startswith('model'):
            if 'TBabs' in c:
                kT_index = i+2
            else:
                kT_index = i+1

        i += 1
    norm_index = kT_index+5
    kT, s2,s3,s4,s5,s6 = contents[kT_index].split()
    norm, s2,s3,s4,s5,s6 = contents[norm_index].split()
    return float(kT), float(norm)
    
def tbvarabs_funk(name,obsid,type,epoch):
    """ Function to extract values for tbvarabs from txt files """
    with open(f"/Users/juliaahlvind/Documents/projekt_1/tbvarabs/scaled_values/tbvarabs_{name}_{obsid}_{type}_{epoch}d.txt") as f:
        contents = f.readlines()
        #mean = contents[2:43]
    ten_percentile = contents[44:85]
    median = contents[86:(86+41)]
    ninty_percentile = contents[87+41:87+2*41]
    tbvarabs_all = [ten_percentile, median, ninty_percentile]
    return tbvarabs_all

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx
