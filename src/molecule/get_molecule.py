# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 22:52:21 2023

@author: win
"""

from ase.data.pubchem import pubchem_atoms_conformer_search
from ase.visualize import view

hmf = pubchem_atoms_conformer_search(cid=69980)
hmf = hmf[1]
hmf.set_pbc(True)
hmf.set_cell([30,30,30])
hmf.center
print(hmf)
hmf.write('POSCAR')
view(hmf)