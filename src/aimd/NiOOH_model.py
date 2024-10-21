# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:41:05 2024

@author: win
"""

from ase.build.tools import cut, stack
from ase.io import read, write
import numpy as np
from ase import Atoms
from ase.build import hcp0001, add_adsorbate, add_vacuum
from ase.constraints import FixBondLengths
from ase.calculators.tip3p import TIP3P, rOH, angleHOH
from ase.md import Langevin
import ase.units as units

from ase.visualize import view
from ase.build import molecule

atoms = read('POSCAR')
niooh_001 = cut(atoms, (1,0,0),(0,1,0),nlayers=10,)



print(niooh_001.get_cell()[2])


def water_volume(atoms, number):
    """
    a*b*np.sqrt(3) /2 * c = total_volume
    total_volum = number_water * vol
    n * a_grid = a
    m * b_grid = b
    l * c_grid = c
    
    n * m * l = number
    
    n = np.sqrt(number / l)
    """
    
    vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))
    total_volume = number * vol
    a,b,c,angle_bc,angle_ac,angle_ab = atoms.get_cell_lengths_and_angles()
    lattice_a = atoms.get_cell()[0]
    lattice_b = atoms.get_cell()[1]
    surface_area = np.linalg.norm(np.cross(lattice_a,lattice_b))
    water_layer_c = total_volume / (surface_area)
    water_length = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
    
    Ga = lattice_a
    Gb = lattice_b
    Gc = np.array([0, 0, water_layer_c])
    n1 = int(a/water_length)
    n2 = int(b/water_length)
    n3 = int(number /(n1*n2))
    # assert n1*n2*n3 == number
    e1 = Ga/n1
    e2 = Gb/n2
    e3 = Gc/n3
    pos = np.zeros((n1,n2,n3,3))
    lpos = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                pos[i][j][k][:] = i*e1 + j*e2 + k*e3
                lpos.append(pos[i][j][k][:])
    
    position = [p.tolist() for p in lpos]
    x = angleHOH * np.pi / 180 / 2
    pos = [[0, 0, 0],
       [0, rOH * np.cos(x), rOH * np.sin(x)],
       [0, rOH * np.cos(x), -rOH * np.sin(x)]]
     
    water = molecule('H2O')
    
    for i in position:
        
        water.rotate(np.random.randint(0,360),[np.random.randint(1,5),np.random.randint(1,5),np.random.randint(1,5)])   
        
        add_adsorbate(atoms, water, height=i[2]+2, position = (i[0]+0.2,i[1]+0.2))
        atoms.center(vacuum=10, axis=2)
        atoms.write('test1.cif')
    return position

a = water_volume(niooh_001, 160)