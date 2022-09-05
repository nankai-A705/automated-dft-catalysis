#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 21:48:54 2022

@author: hongyeqin
"""


#rom ase.qhyio import Atoms_to_slab
from ase.db import connect
from generateSlab import Generater_Slab_with_Pymatgen
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from ase.build import surface
def get_file_name(**kwargs):
    
    """
    input  : get_file_name(Cu=0.2, Ni=0.3)
    return : {'Cu': 0.2, 'Ni': 0.3}-1.db

    """
    file = {} # key : atoms name value : fraction ['Pt'=0.1, 'Pd'=0.2]
    for key, value in kwargs.items():
        file[key] = value    
    png_name = str(file)
  
    name = str('{}.db'.format(png_name))
    
    return name

def get_id(file):
    
    bulk = connect(file)
    id_list = []
    for row in bulk.select():
        id_list.append(row.id)
    return id_list, len(id_list)

def check_state(atoms):
    
    """
       Return slab state surface with or without adsorbates
    """
    tags = atoms.get_tags().tolist()    
    b = min([tags.count(j) for j in tags])
    c = max([tags.count(j) for j in tags])
    t = [tags.count(j) for  j in tags]
    if c - b >= 2:
        tag = tags[t.index(min(t))] # return adsorbates tag 
        return str('adsorbate'), tag
    else:
        return str('surface')
    
def get_label(atoms):
    
    state = check_state(atoms)[0]
    if state == 'surface':
        
       atom_index = [atom.index for atom in atoms]
       symbol = [atoms.symbols[i] for i in atom_index]
       return symbol
    else:
        tag = check_state(atoms)[1]
        atom_index = [atom.index for atom in atoms if atom.tag != tag]
        symbol = [atoms.symbols[i] for i in atom_index]
        return symbol

def get_atom(name, idx):
    
    """
    Get Atom object according to id 
    
    """
    db = connect(name)
    atoms = db.get_atoms(id=idx)
    
    return atoms

def get_slab(atoms,miller_index,min_slab_size,min_vacuum_size):
    
    slabs = Generater_Slab_with_Pymatgen(atoms, miller_index, min_slab_size, min_vacuum_size, surface_symmetry = None)
    
    return slabs

def get_slab_db(atoms,miller_index,min_slab_size,min_vacuum_size,filename=None, **kwargs):
    
    if filename is None:
       file_name = get_file_name(**kwargs)
       state = check_state(atoms)
       if file_name.isspace == True:
          print('You should input atoms species and ratio like {Ru=0.5, Ni=0.5}')
    else: 
        file_name = filename
        state = check_state(atoms)
        if state == 'surface':
           surface_db = connect('{}-{}.db'.format(state,file_name))
           slab =  get_slab(atoms,miller_index,min_slab_size,min_vacuum_size)
           for sla in slab:
               slab = AseAtomsAdaptor().get_atoms(sla) # Tranform to Atoms type
               symbol = get_label(slab)
               surface_db.write(slab,relaxed=False,data={'symbols':symbol})
           return surface_db
        elif state =='adsorbate':
            surface_db = connect('{}-{}.db'.format(state,file_name))
            slab =  get_slab(atoms,miller_index,min_slab_size,min_vacuum_size)
            for sla in slab:
                slab = AseAtomsAdaptor().get_atoms(sla)
                symbol = get_label(slab)
                surface_db.write(slab,relaxed=False, data={'symbols':symbol})
            return surface_db
        else:
            print('Please make sure with right crystal Structure')
        


def get_surface_db(atoms,lattice,layer,vacuum):
    s = surface(atoms,lattice,layer)
    s.center(vacuum=vacuum, axis=2)
    
    return s

filename = 'RuNi-bulk.db'
for i in range(get_id(filename)[1]):
    #atoms = read('{}@id={}'.format(filename, i))
    atoms = get_atom(filename, i+1)
    label = get_label(atoms)
    slab = get_slab_db(atoms, (0,0,1), 4, 15,filename) 
# atoms = get_atom(filename, 3)
# atoms.write('p.cif')
# #     label = get_label(atoms)
# slab = get_surface_db(atoms, (0,0,1),1,10) 
# print(slab)


# slab.write('POSCAR')

       




