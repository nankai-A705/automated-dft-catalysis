#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 19:09:11 2022

@author: hongyeqin
@email:290720931@qq.com
"""
from ase.db import connect
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from xml_parser import get_pdos
import numpy as np
import time
import os
"""
                                                    
"""
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
    if mode == 'center':
      return [moment(energies[:], dos_rev[:], 1)]
    
    
def get_atom(name, idx):
    
    """
    Get Atom object according to id 
    
    """
    db = connect(name)
    atoms = db.get_atoms(id=idx)
    return atoms
    
def get_id(file):
    
    bulk = connect(file)
    id_list = []
    for row in bulk.select():
        id_list.append(row.id)
    return id_list, len(id_list)

file = 'surface-RuNi-bulk.db'
db = connect(file)
relaxed_db = connect('RuNi-relaxed.db')
vasp = Vasp(ibrion=2,
            nsw=500,
            isif=2,
            ispin=2,
            isym=-1,
            ediff=1e-5,
            ediffg=-0.02,
            algo='Fast',
            lasph=True,
            ismear=1,
            lreal='Auto',
            sigma=0.2,
            kpts=(5, 5, 1),
            encut=600,
            nelm=60,
            ivdw=11,
            lorbit=11,
            lwave=False,
            lcharg=False,
            ldipol=True,
            idipol=3,
            xc='pbe',
            ncore=16,
            lmaxmix=4)
number_of_atoms = get_id(file)[1] 
for i in range(number_of_atoms):

    atoms = get_atom(file,i+1)
    row = db.get(id=i+1)
    atom_index = [atom.index+1 for atom in atoms]
    if os.path.exists('vasprun.xml'):
        os.remove('vasprun.xml')
        os.remove('WAVECAR')
        os.remove('CHGCAR')
        os.remove('CHG')
        atoms.calc = vasp
        time.sleep(30)
        energy = atoms.get_potential_energy()
        d_band = []
        for j in atom_index:
          energy = get_pdos('./','vasprun.xml',j,['dxy'])[0]
          dos = get_pdos('./','vasprun.xml',j,['dxy','dyz','dxz','x2-y2','dz2'])[1]
          d = density_moments(energy, dos,'center')
          d_band.append(d)
        relaxed_db.write(atoms,relaxed=True, data={'d_band':d_band})

    time.sleep(30)
# atoms = get_atom(file, 1)
# atoms.calc = vasp
# dyn = BFGS(atoms, trajectory='H2O.traj')
# dyn.run(fmax=0.05)
