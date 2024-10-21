# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 21:27:41 2021

@author: win
"""

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

Ru = hcp0001('Ru', (4,4,4),a=10.91342)
a = Ru.cell[0]
b= np.sqrt(np.sum(a ** 2))
print(Ru)

def get_water_layer(d):
    x = angleHOH * np.pi / 180 / 2
    pos = [[0, 0, 0],
       [0, rOH * np.cos(x), rOH * np.sin(x)],
       [0, rOH * np.cos(x), -rOH * np.sin(x)]]
    atoms = Atoms('OH2', positions=pos)
    vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.) #single water side length 0.998mol/L
    layer = int(d / vol)
    #total = int(d**3 / vol**3)
    z_new = d**3 / (d**2 * np.sqrt(3) /2)
    atoms.set_cell((vol, vol, vol))
    atoms.center()
    atoms = atoms.repeat((layer, layer, int(z_new/vol)))
    atoms.set_cell((d,d,z_new,90,90,60))
    #atoms.set_cell((d,d,d))
    #atoms.set_pbc(True)
    return atoms

    

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
    Ru = hcp0001('Ru', (8,8,4),a=2.728355)
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


a = water_volume(6)
print(len(a))

# x = angleHOH * np.pi / 180 / 2
# pos = [[0, 0, 0],
#        [0, rOH * np.cos(x), rOH * np.sin(x)],
#        [0, rOH * np.cos(x), -rOH * np.sin(x)]]
# water = Atoms('OH2', positions=pos)
# water.rotate(np.random.randint(0,90),'z')
# view(water)
# from ase.build import molecule
# w = molecule('H2O')
# w.rotate(np.random.randint(0,90),'z')
# view(w)
# view(c)
# print(c)
# print(np.random.randint(0,90))
# a = 12

# Ga = np.array([a, 0, 0])
# Gb = np.array([a/2, (a * np.sqrt(3))/2, 0])
# Gc = np.array([0, 0, 9])

# n1 = 3
# n2 = 3 
# n3 = 2
# e1 = Ga/n1
# e2 = Gb/n2
# e3 = Gc/n3
# lpos = []
# pos = np.zeros((n1,n2,n3,3))
# for i in range(n1):
#       for j in range(n2):
#           for k in range(n3):
#               pos[i][j][k]= i*e1 + j*e2 + k*e3
#               lpos.append(pos[i][j][k])
# lpos = pos.tolist()
# print(lpos)
# print(e1)

# water = get_water_layer(11.08743432) 
# water.write('water.cif')   
# from ase.build import stack  
# from ase.build import fcc111
# from ase.build import add_adsorbate
# slab1 = fcc111('Pt', size=(4,4,6))
# add_adsorbate(slab1,'H',1.5,'bridge')
# #slab1.center(vacuum=0, axis =2)

# a = slab1.get_cell()[0]
# b = slab1.get_cell()[1]
# slab1.write('slab1.cif')
# cosangle = a.dot(b)/(np.linalg.norm(a) * np.linalg.norm(b))
# inv = np.arccos(cosangle)
# angle = np.degrees(inv)  
# print(angle)
# slab2 = fcc111('Al', size=(3,3,4))
# interface = stack(slab1, water, maxstrain=1, distance=2.3)
# print(slab1)
# interface.write('heterostructure.cif')    
