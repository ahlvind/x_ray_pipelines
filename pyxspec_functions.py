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
    
    return min_index, index_up, index_low, new_lower_index, new_upper_index

def steppar_func(step_start, step_end, steps, paramNr, max_call_limit, force_upper_limit):
    """
    Execute the `Fit.steppar` function with the specified input parameters and handle any errors.

    This function performs a parameter step search using `pyxspec`'s `Fit.steppar` for a given parameter,
    and returns the results as arrays. It limits the number of calls to prevent excessive iteration and
    handles errors that may occur during the step operation.

    Parameters
    ----------
    step_start : float
        The starting value for the steppar search.
    step_end : float
        The ending value for the steppar search.
    steps : int
        The number of steps between `step_start` and `step_end`.
    paramNr : int
        The parameter number to be stepped.
    max_call_limit : int
        Maximum number of function calls allowed to avoid excessive iterations.
    force_upper_limit : bool
        A flag indicating whether to force the upper limit during the steppar operation.
        
    Returns
    -------
    cflux_array : numpy.ndarray or None
        Array of `cflux` values from the `Fit.stepparResults` if successful; otherwise `None`.
    delcstat_array : numpy.ndarray or None
        Array of `delstat` values from the `Fit.stepparResults` if successful; otherwise `None`.
        
    Raises
    ------
    ValueError
        If `max_call_limit` is less than or equal to 0, indicating that the maximum number of calls has been reached.

    Notes
    -----
    The function logs the number of iterations using `max_call_limit` and limits it to avoid excessive recursive calls.
    It returns `None` for both arrays if an error occurs in the `Fit.steppar` call.

    Examples
    --------
    >>> cflux, delstat = steppar_func(0.1, 10.0, 50, 2, 20, False)
    >>> if cflux is not None:
    ...     print(f"Cflux array: {cflux}")
    ... else:
    ...     print("Steppar operation failed.")
    """
    print(f'itteration nr: {max_call_limit}')
    if max_call_limit <= 0:
        raise ValueError(f"Exceeded the maximum number of function calls ({max_call_limit}).")
    
    try:
        Fit.steppar(f"{paramNr} {step_start} {step_end} {steps}")
    except Exception as e:
        print(f"Error in Fit.steppar: {e}")
        return None, None

    cflux_array = np.array(Fit.stepparResults(paramNr))
    delcstat_array = np.array(Fit.stepparResults('delstat'))
    
    return cflux_array, delcstat_array



def find_indices_for_steppar(cflux_array, delcstat_array, tolerance, target, steps, paramNr, max_call_limit, force_upper_limit=False):
    """
    Determine the best flux, (lower if exists) and upper flux limits that correspond to a given target Δcstat value using `Fit.steppar()` results.

    This function identifies the index positions within the `cflux_array` and `delcstat_array` that correspond to a specified target value for Δcstat, within a given tolerance. If the lower or upper limit is not found within the current data range, the function recursively zooms in on the region around the best Δcstat value, performing additional `Fit.steppar()` operations.

    Parameters
    ----------
    cflux_array : numpy.ndarray
        Array of flux values corresponding to a given parameter from the `steppar` search.
    delcstat_array : numpy.ndarray
        Array of Δcstat values corresponding to each Δcstatistics value.
    tolerance : float
        The acceptable tolerance of |Δcstat-target|<tolerance.
    target : float
        The target Δcstat value to find within the `delcstat_array`.
    steps : int
        Number of steps to use in each `steppar` operation.
    paramNr : int
        Which model parameter is used for the `Fit.steppar` function call.
    max_call_limit : int
        Maximum number of recursive calls of the whole steppar funciton allowed to prevent excessive iterations.
    force_upper_limit : bool, optional
        A flag indicating whether to ignore the search for a lower limit if it does not exist (default is `False`). This parameter becomes important for when searching upper values in steppar called a second time when array has zoomed in on upper limit region.

    Returns
    -------
    cflux : float
        The best-fit flux value corresponding to the minimum Δcstat value.
    lower_cflux : float
        The flux value corresponding to the lower limit for the target Δcstat value, or 0 if not found.
    upper_cflux : float
        The flux value corresponding to the upper limit for the target Δcstat value.

    Raises
    ------
    ValueError
        If `max_call_limit` is less than or equal to 0, indicating the maximum number of recursive calls has been reached.
    TypeError
        If either `cflux_array` or `delcstat_array` is not a NumPy array.

    Notes
    -----
    - If the best-fit index is at position 0 or if the lower limit does not exist (based on `force_upper_limit`),
      the function only searches for the upper limit.
    - If the tolerance condition is not met, the function recursively refines the search region until the target value is found or the `max_call_limit` is reached.
    - The `steppar_func` is called to perform additional `Fit.steppar` operations in the refined range during recursion.

    Examples
    --------
    >>> cflux_array = np.array([0.1, 0.5, 1.0, 1.5, 2.0])
    >>> delcstat_array = np.array([12.0, 10.5, 8.0, 9.5, 11.0])
    >>> cflux, lower, upper = find_indices_for_steppar(cflux_array, delcstat_array, tolerance=0.1, target=9.0, steps=5, paramNr=2, max_call_limit=10)
    >>> print(f"Best-fit flux: {cflux}, Lower limit: {lower}, Upper limit: {upper}")
    Best-fit flux: 1.0, Lower limit: 0, Upper limit: 1.5
    """
    
    print('Δcstat target value: ', target)
    best_index = np.where(delcstat_array == np.min(delcstat_array))[0][0]  # Find the index of min Δcstat
    print(f'Forced upper limit: {force_upper_limit}')
    
    if best_index == 0 or delcstat_array[0]<target or force_upper_limit:
        print('The lower limit does not exist, searching for only upper limit.')
        
        if target==9.0:
            print('Exiting function, no search for only upper value of 9.0!')
            return 1, 1, 1
            
        # Only search for the upper limit
        index_up = best_index + min(range(len(delcstat_array[best_index:])),
                                    key=lambda i: abs(delcstat_array[best_index + i] - target))
        print('tolerance check:',abs(delcstat_array[index_up] - target))
        if abs(delcstat_array[index_up] - target) <= tolerance:
            upper_cflux = cflux_array[index_up]
            print('|Δcstat-target|< tolerance')
            return upper_cflux, 0, upper_cflux  # No lower limit
        else:
            print('|Δcstat-target|> tolerance')
            # Recursive call if tolerance not met. Zoom in on the region around max to find only upper limit
            new_step_end = cflux_array[min(index_up + 3, len(cflux_array) - 1)]
            new_step_start = cflux_array[min(index_up -3, len(cflux_array) - 1)]
            cflux_array, delcstat_array = steppar_func(new_step_start, new_step_end, steps, paramNr, max_call_limit - 1, force_upper_limit=True)
            return find_indices_for_steppar(cflux_array, delcstat_array, tolerance, target, steps, paramNr, max_call_limit-1, force_upper_limit=True)

    else:
        # Find the best flux
        cflux = cflux_array[best_index]

        # Find the upper limit
        index_up = best_index + min(range(len(delcstat_array[best_index:])),
                                    key=lambda i: abs(delcstat_array[best_index + i] - target))
        
        if abs(delcstat_array[index_up] - target) <= tolerance:
            upper_cflux = cflux_array[index_up]
        else:
            new_step_start = cflux_array[best_index]
            new_step_end = cflux_array[min(index_up + 3, len(cflux_array) - 1)]
            return steppar_func(new_step_start, new_step_end, steps, paramNr, max_call_limit - 1)

        # Find the lower limit
        index_low = min(range(best_index), key=lambda i: abs(delcstat_array[i] - target))
        
        if index_low == 0:  # No lower limit
            return upper_cflux, 0, upper_cflux

        if abs(delcstat_array[index_low] - target) <= tolerance:
            lower_cflux = cflux_array[index_low]
            return cflux, lower_cflux, upper_cflux
        
        # Recursive call if tolerance not met for lower limit
        return steppar_func(cflux_array[best_index], cflux_array[index_up + 3], steps, paramNr, max_call_limit - 1)





def call_steppar(cstat, paramNr, step_start, step_end, **kwargs):
    """
    Call the `steppar_func` to perform a parameter step search and find the flux values corresponding to a target Δcstat.
    This function performs a `steppar` search for a given parameter over a specified range and identifies the best-fit flux, as well as the flux values corresponding to the lower and upper limits of a target Δcstat value. It uses the `steppar_func` to get the initial results and `find_indices_for_steppar` to analyze the results and refine the search if necessary.

    Parameters
    ----------
    cstat : float
        The target Δcstat value for which the flux limits should be found.
    paramNr : int
        The parameter number to be used in the `steppar` function.
    step_start : float
        The starting value for the parameter in the step search.
    step_end : float
        The ending value for the parameter in the step search.
    **kwargs : dict, optional
        Additional keyword arguments to be passed to other functions, if necessary.

    Returns
    -------
    cflux : float
        The best-fit flux value corresponding to the minimum Δcstat value.
    cflux_under : float
        The flux value corresponding to the lower limit for the target Δcstat value.
    cflux_over : float
        The flux value corresponding to the upper limit for the target Δcstat value.

    Notes
    -----
    - The function uses a fixed tolerance value of `0.1` and a default number of steps (`100`) for the initial `steppar` search.
    - The `steppar_func` is called with the maximum number of allowed iterations set to `20` to avoid excessive calls.
    - This function is only called once, further refinements in steppar procedure is incorporated in functions called by this funciton.

    Examples
    --------
    >>> cstat = 9.0
    >>> paramNr = 2
    >>> step_start = 0.1
    >>> step_end = 10.0
    >>> cflux, lower_limit, upper_limit = call_steppar(cstat, paramNr, step_start, step_end)
    >>> print(f"Best-fit flux: {cflux}, Lower limit: {lower_limit}, Upper limit: {upper_limit}")
    Best-fit flux: 1.0, Lower limit: 0.5, Upper limit: 2.0
    """
    
    tolerance = 0.1
    steps = 100
    cflux_array, delcstat_array = steppar_func(step_start, step_end, steps, paramNr, 20, False)
    print(cflux_array, delcstat_array)
    cflux, cflux_under, cflux_over = find_indices_for_steppar(cflux_array, delcstat_array, tolerance, cstat, steps, paramNr, 20)
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
    #with open(f"/Users/juliaahlvind/Documents/projekt_1/tbvarabs/scaled_values/tbvarabs_{name}_{obsid}_{type}_{epoch}d.txt") as f:
    with open(f"/Users/juliaahlvind/Documents/projekt_1/tbvarabs/scaled_values/tbvarabs_{name}_{type}_{epoch}d_G{obsid}.txt") as f:
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
