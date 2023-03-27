# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 09:02:15 2023

@author: win
"""

from mayavi import mlab
import numpy as np
from ase.io import read


# d_bands = np.loadtxt("dbands.txt", dtype='int32', delimiter=',')
file = open('dbands.txt')
value = file.readlines()[0]

v = eval(value)

# file = pd.read_table('dbands.txt', sep=',')
# print(file)
# d = [float(x) for x in value]


atoms = read('CONTCAR')
atom_index = [atom.index for atom in atoms]
print(atom_index)
postions = atoms.positions
x = [pos[0] for pos in postions]
y = [pos[1] for pos in postions]
z = [pos[2] for pos in postions]
# print(res)
# #symbol = [surface.symbols[i] for i in atom_index]
# print(x)

mlab.points3d(x,y,z,v,line_width=4, scale_mode='vector',resolution=20)
mlab.show()