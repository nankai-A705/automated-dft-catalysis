# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 11:55:58 2024

@author: hongye qin
"""
from ase.visualize import view

from ase.cluster.cubic import FaceCenteredCubic
from ase.visualize import view
surfaces = [(1, 1, 1), (1, 1, 1), (1, -1, 1)]
layers = [2, 2, -1]

cluster = FaceCenteredCubic('Ag', surfaces, layers, vacuum=10)
cluster.center(vacuum=10, axis=(0,1,2))
cluster.set_pbc([True])
# cluster.add_vacuum(10)
cluster.write('Ag-cluster-111-111-1-11.cif')
view(cluster)