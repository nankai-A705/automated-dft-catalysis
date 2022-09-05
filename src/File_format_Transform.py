#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import os
import numpy as np
from ase.build import fcc111
from ase.io import read
from pymatgen import vis
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.operations import SymmOp
from pymatgen.core.surface import generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import in_coord_list_pbc
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor


# In[2]:


class File_Format_Transform():
    
    def __init__(self,filename):
        self.filename = filename
        
        pass
    def Transform_to_Atoms(self,filename):
        sturcture = Structure.from_file(filename)
        atoms = AseAtomsAdaptor.get_atoms(sturcture)
        return atoms
    
    def Transform_to_Structure(self,filename):
        atoms = read(filename)
        structure = AseAtomsAdaptor.get_structure(atoms)
        return structure
        


# In[ ]:




