# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 19:06:15 2022

@author: win
"""
from ase import Atoms
from ase.io import read
d = 1.1
co = Atoms('COO', positions=[(0, 0, 0), (0, 0, d),(0, d, 0)])
print(co[1].symbol)
a = read('POSCAR')
def magmon_array(atoms, **kwargs):
    """
    array = magmon_array(Atoms, element1=magmon1, element2=magmon2,....)
    
    """
    file = {}
    for key, value in kwargs.items():
        file[key] = value
    atoms_symbols = [atom.symbol for atom in atoms]
    index = [atom.index for atom in atoms]
    array = []
    for i in atoms_symbols:
        
        for k, v in file.items():
            
            if i == k:
                array.append(file[k])
                
        
        
        
    return atoms_symbols, array

b = magmon_array(a, O=1, Mn=3, Co=5)
print(b)
# a.set_initial_magnetic_moments(b)
# print(a.get_initial_magnetic_moments())