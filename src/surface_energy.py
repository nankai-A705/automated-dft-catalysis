# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:09:49 2021

@author: win
"""
import os
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.io.ase import  AseAtomsAdaptor
from ase.db import connect
from ase.build import bulk
from ase.io import read
from ase.calculators.vasp import Vasp
import subprocess
from ase import io


non_metal_atoms_energy= {'H':-3.3574, 'O':-4.9480}
metal_atoms_energy = {}
None_metal_elements = ['H','O','N','C']
Metal_elements = ['Fe','Co','Ni','Mn']

def make_dir_name(atoms):
    name = atoms.get_chemical_formula()
    os.mkdir('./{}'.format(name))
    
def run_vasp_bulk(cell):
    name = cell.get_chemical_formula()
    calc = Vasp(prec='Normal',
                xc='PBE',
                isif=2,
                encut=520,
                lreal=False,
                ispin=2,
                kts=(9,9,9),
                directory='./{}_bulk'.format(name))
    cell.calc = calc
    energies = cell.get_potential_energy()
    return energies
def run_vasp_surface(cell, miller_index):
    #surface_db = connect('surface.db')
    name = cell.get_chemical_formula()
    calc = Vasp(prec='Normal',
                xc='PBE',
                isif=2,
                encut=520,
                potim=0.1,
                lreal=False,
                ispin=2,
                kpts=(get_kpoints(cell)[0],get_kpoints(cell)[1],1),
                ldipol=True,
                idipo=3,
                directory='./{}_{}surface'.format(name, miller_index))
    cell.calc = calc
    cell.get_potential_energy()
    #surface_db.write(cell, miller='{}'.format(miller_index))
    #os.chdir('../')
    
def get_kpoints(atoms):
    """
    Parameters
    ----------
    atoms : Atoms
    Ka * Vector_a = 30 A
    """
    cell = atoms.get_cell()
    ka = int(30.0 / np.linalg.norm(cell[0]))
    kb = int(30.0 / np.linalg.norm(cell[1]))
    kc = int(30.0 / np.linalg.norm(cell[2]))
    return ka, kb, kc

def calculate_bulk_chemical_potential(atoms):
    """
    atoms type : Atoms ase
    
    """
    if isinstance(atoms, str):
        atoms = bulk(atoms)
    atoms_type = parse_formula(atoms)[0]
    if len(atoms_type) > 1:
        print('You should know what are you doing, chemical potential from single phase generally')
        
    energy = run_vasp_bulk(atoms)
    subprocess.call('cp -rf OUTCAR OUTCAR_$(date +%s)', shell=True)
    io.write('out.cif', atoms)
    atoms_number = len(atoms.get_atomic_numbers())
    
    return energy/atoms_number

#def parse_formula(atoms):
#    """
#    atoms type : Atoms ase
#    warnning: NiO2 Will not be recognized
#    """
#    formula = atoms.get_chemical_formula()
#    name = re.split(r'(\d+)', formula)
#    while '' in name:
#        name.remove('')
#    chemical_species = []
#    species_number = []
#    for i in range(((len(name))//2)):
#        chemical_species.append(name[2*i])
#        species_number.append(name[2*i + 1])
#    return chemical_species, species_number

def parse_formula(atoms):
    """
    atoms type : Atoms ase
    return : chemical species ['Fe','Co']
             species number   ['2','1']
    
    """
    formula = atoms.get_chemical_symbols()
    formula_dict = {}
    for i in formula:
        formula_dict[i] = formula.count(i)
    chemical_species = []
    species_number = []   
    for key in formula_dict:
        chemical_species.append(key)
        species_number.append(formula_dict[key])
    return chemical_species, species_number

def surface_area(atoms):
    """
    Calculates the surface area of the slab
    atoms type : Atoms ase
    """
    cell = atoms.get_cell()
    return np.linalg.norm(np.cross(cell[0], cell[1]))
    

def get_index(file, name):
    """
    
    get the same atoms index in metal_dict and file
    
    Parameters
    ----------
    file : list
    name : elements in list

    """
    index = [x for (x,y) in enumerate(file) if y==name ]
    
    return index

def bulk_atoms_potential(atoms, vasp=True,*arg):
    """
    gama = (Eslab - N*Ebulk) / 2* surface_Area
    return N*Ebulk
    usage: bulk_atoms_potential(Atoms,'Fe')
    
    """

    if isinstance(arg, tuple):
       for i in range(len(arg)):
           new_bulk = bulk(arg[i])
           print(new_bulk)
           #atoms_energy.append(calculate_bulk_chemical_potential(new_bulk))
    chemical_species, species_number = parse_formula(atoms)
    non_metal_element = [element for element in chemical_species if element in None_metal_elements]
    metal_element = [element for element in chemical_species if element in Metal_elements]
    
    non_metal_chemical_potential=0 #calculate non-metal chemical potential. All the data from Material Project
    for i in range(len(non_metal_element)):
        
        non_metal_chemical_potential += non_metal_atoms_energy[non_metal_element[i]] * species_number[get_index(chemical_species, non_metal_element[i])[0]]
        
    if vasp:  # weather use vasp to get chemical potential or not! 
       metal_chemical_potential = 0
       for i in range(len(metal_element)):
           metal_chemical_potential += calculate_bulk_chemical_potential(metal_element[i]) * species_number[get_index(chemical_species, metal_element[i])[0]]
    else:    # this module is not available now
        metal_chemical_potential = 0
        for i in range(len(metal_element)):
            metal_chemical_potential += metal_atoms_energy(metal_element[i]) * species_number[get_index(chemical_species, metal_element[i])[0]]
    
    return metal_chemical_potential + non_metal_chemical_potential
        
def get_slab_with_miller(structure, max_index, min_slab_size, min_vacuum_size,primitive=False,symmetrize=True):
    """
    Get surface slab and miller index

    """
    #surface_db = connect('surface.db')
    if isinstance(structure, Structure):
       slabs = generate_all_slabs(structure, max_index, min_slab_size, min_vacuum_size,primitive=False,symmetrize=True)
    ase_slab = []
    miller_index = []
    for slab in slabs:
        miller_index.append(slab.miller_index)
        ase_slab.append(AseAtomsAdaptor.get_atoms(slab))
    return ase_slab, miller_index

def get_slab_energy(cell, mill_index):
    run_vasp_surface(cell, miller_index)
    
    
        
        
       
surface_db = connect('surface.db')       
         
cell1 = Structure.from_file('NiO_mp-19009_conventional_standard.cif') 

ase_slab, miller_index = get_slab_with_miller(cell1, 1, 5, 15,primitive=False,symmetrize=True)
    

print(ase_slab, miller_index)

for i in range(len(ase_slab)):
    run_vasp_surface(ase_slab[i], miller_index[i])
    surface_db.write(ase_slab[i], miller='{}'.format(miller_index[i]))





