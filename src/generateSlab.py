# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 21:11:23 2021

@author: win
"""
import os
from ase import Atoms
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor


def Generater_Slab_with_Pymatgen(filename,miller_index, min_slab_size, min_vacuum_size, surface_symmetry = None, transform_primitive_to_unit = None):
    
    if isinstance(filename, Structure):
            
       structure = Structure.from_file(filename)
              
       
       if transform_primitive_to_unit:
          unit_cell = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
       
          slab = SlabGenerator(unit_cell, miller_index, min_slab_size, min_vacuum_size,center_slab=False, primitive=False)
          slabs = slab.get_slabs()
       else:
           slab = SlabGenerator(structure, miller_index, min_slab_size, min_vacuum_size, center_slab=False,primitive=False)
           slabs = slab.get_slabs()
           
           #slabs = [slab for slab in slabs if  slab.is_symmetric() is True]
       return slabs
       
    if isinstance(filename, Atoms):
        
       structure = AseAtomsAdaptor.get_structure(filename)
       
       if transform_primitive_to_unit:
          unit_cell = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
       
          slab = SlabGenerator(unit_cell, miller_index, min_slab_size, min_vacuum_size,center_slab=False, primitive=False)
          slabs = slab.get_slabs()
       else:
           slab = SlabGenerator(structure, miller_index, min_slab_size, min_vacuum_size, center_slab=False,primitive=False)
           slabs = slab.get_slabs()
           #slabs = [slab for slab in slabs if  slab.is_symmetric() is True]
       return slabs 
   

