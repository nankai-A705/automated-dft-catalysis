# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 15:39:43 2021

@author: win
"""
from pymatgen.core.surface import Slab
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.ase import AseAtomsAdaptor


try: 
    from ase.atoms import Atoms
    from ase.atom import Atom
    ase_loaded = True
    
except ImportError:
    ase_loaded = False
    
def Atoms_to_Slab(atoms, miller_index,shift,scale_factors,cls = None,):
    """
    This function could transform Atoms type from ase.build.fcc111  module to 
    pymatgen Slab type
    Usage:
        miller_index:array
        shift:float
        scale_factors: array [1,1,1] which means a*1,b*1,c*1
        miller_index:array [1,1,1]
    """
    
    atoms.center(vacuum=0, axis=2)

    oriented_unit_cell = AseAtomsAdaptor().get_structure(atoms)
    
    #try:
    #    atoms.get_cell()
        
    #except:
    #    atoms.center(vacuum=0, axis=2)

    
    symbols = oriented_unit_cell.species
    positions = atoms.get_positions()
    lattice = oriented_unit_cell.lattice
    
    oriented_unit_cell = AseAtomsAdaptor().get_structure(atoms)
    cls = Slab if cls is None else cls
    return cls(lattice, symbols, positions, miller_index, oriented_unit_cell, shift, scale_factors,coords_are_cartesian=True,)
    
        








