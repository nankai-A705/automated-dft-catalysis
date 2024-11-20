# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:25:12 2024

@author: Win
"""

from ase import Atom, Atoms
from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms, Hookean
from ase.optimize.minimahopping import MinimaHopping
from ase import Atoms
from ase.build import hcp0001, add_adsorbate, add_vacuum
from ase.constraints import FixBondLengths
from ase.calculators.tip3p import TIP3P, rOH, angleHOH
from ase.md import Langevin
import ase.units as units
from ase.io import write, read
from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.build import molecule
import numpy as np
from ase.cluster.cubic import FaceCenteredCubic
from ase.visualize import view
surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [2, 2, -1]

cluster = FaceCenteredCubic('Ag', surfaces, layers)


atoms = fcc111('Cu',(6,6,4),a=3.6346, vacuum=7.5)




def water_volume(layer):
    """
    a*b*np.sqrt(3) /2 * c = total_volume
    total_volum = number_water * vol
    n * a_grid = a
    m * b_grid = b
    l * c_grid = c
    
    n * m * l = number
    
    n = np.sqrt(number / l)
    """
    Ru = fcc111('Cu',(6,6,4),a=3.6346)
    cell_a = Ru.get_cell()[0][0]
    print(cell_a)
    vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))
    water_length = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
    n1 = int(cell_a/water_length)
    n2 = n1 
    total_volume = layer * vol * n1 * n2
    c = total_volume / (cell_a**2 * np.sqrt(3) /2 )

    
    Ga = np.array([cell_a, 0, 0])
    Gb = np.array([cell_a/2, (cell_a * np.sqrt(3))/2, 0])
    Gc = np.array([0, 0, c])

    n3 = layer
    print(n1)
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
    # Ru = hcp0001('Ru', (6,6,4),a=2.728355)
#    add_adsorbate(Ru, 'H',height=1, position='hcp')

    for i in position:
        
        water.rotate(np.random.randint(0,360),[np.random.randint(1,5),np.random.randint(1,5),np.random.randint(1,5)])   
        
        add_adsorbate(Ru, water, height=i[2]+2.5, position = (i[0]+(np.random.rand()),i[1]+np.random.rand()))
        Ru.center(vacuum=10, axis=2)
        Ru.write('test884.cif')
    return position

water_volume(9)