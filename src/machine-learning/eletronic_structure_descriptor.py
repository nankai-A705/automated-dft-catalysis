# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 20:10:47 2022

@author: hongyeqin

"""
import numpy as np


def smooth(y, box_pts):
    """Smooth the noisy density of state distribution using convolution operater.
    Args:
       y (array): one-dimensional input array
       box_pts (int): average box size parameter
    Returns:
       array: one-dimentional array of smoothed density of state values
    """
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def moment(x, y, n):
    """moments function to calculate the porbability distribution characteristics of density of states.
    Args:
       x (array): energy values
       y (array): density of states
       n (int): order parameter of moments function
    Returns:
       float: moment descriptor 
    """
    p = x**n * y
    return np.trapz(p, x)/np.trapz(y, x)

def density_moments(energies, dos, mode=None):
    """Calculate the moment descriptors for the density of states distributions.
    Args:
       energies (array): energy values with respect to fermi energy
       dos (array): density of states
    returns:
       list: moment characteristics including filling, center, sigma, skewness, kurtosis
    """
    # smooth the noisy density state 
    dos_rev = smooth(dos[:], 15)# smoothing function
    # determine the index of first non-positive value  
    Ind = np.argmax(energies>0)
    # calculate the moment descriptors for the density of states 
    if mode is None:
       filling = np.trapz(dos_rev[0:Ind], energies[0:Ind])/np.trapz(dos_rev, energies)
       center = moment(energies[:], dos_rev[:], 1)             
       sigma_c = np.sqrt(moment(energies[:]-center, dos_rev[:], 2))
       skewness = moment(energies[:]-center, dos_rev[:], 3)/sigma_c**3 
       kurtosis = moment(energies[:]-center, dos_rev[:], 4)/sigma_c**4 
       return [filling, center, sigma_c, skewness, kurtosis]
    if mode is center:
        return [moment(energies[:], dos_rev[:], 1)]
    