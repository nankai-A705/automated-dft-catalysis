# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 22:40:24 2022

@author: hongye
"""
import numpy as np
from ase.io.trajectory import Trajectory
from ase.cluster import Octahedron
from ase.calculators.vasp import Vasp
from ase.io import read
from ase.visualize import view
from ase import Atoms
import time


from ase.db import connect

# scog = SCOG(atoms, elements=['Co', 'Mn', 'Fe', 'Cu', 'Mo','Sn'],
#             symmetry='spherical',
#             composition={'Co':0.166, 'Mn':0.166, 'Fe':0.166, 'Cu':0.166, 'Mo':0.166,'Sn':0.17})
# scog.run(max_gen=25, mode='stochastic', verbose=True)

# images = read('orderings.traj', index=':')
#view(atoms)
# atoms = Atoms(atoms)

def get_random_hcp_alloy(num_structure, host,  **kwargs):
    """
    generate high entropy alloy
    output type: type1 db
                 type2 POSCAR
    this module only works for hcp
    """
    
    
    file = {} # key : atoms name value : fraction ['Pt'=0.1, 'Pd'=0.2]
    for key, value in kwargs.items():
        file[key] = value
      
    key_dict = []  # type : dict ['Pt','Pd']
    value_dict = []# type : dict [0.1, 0.2]
    
    
    
    
    #assert(sum(value_dict) == 1)
    for key in file:
        key_dict.append(key)
        value_dict.append(file[key])
    png = key_dict + value_dict
    png_name = str(png)
    dopant_species = len(key_dict)

    db = connect('{}.db'.format(png_name))
    traj = Trajectory('{}.traj'.format(png_name),'w')
    for i in range(num_structure):
        
        
        seed=1120210312+i*100
        atoms =Octahedron(host, 5,2)
        atoms.center(vacuum=10, axis=(0,1,2))
        atoms.set_pbc([True])

        super_cell = atoms
        positions = super_cell.get_positions()
        cell = super_cell.cell
        angles = cell.angles()
        
        
            
        atoms_ordinal = super_cell.get_atomic_numbers() # get supercell atoms number
        atoms_num = len(atoms_ordinal)
        elements = super_cell.get_chemical_symbols() #[Al, Al, Al ,etc]
        impurity = [] # dict for dopants numbers 
        np.random.seed(seed) # for repeat
        
        for j in range(dopant_species):
            impurity.append(np.round(atoms_num * value_dict[j]))
        for k in range(len(impurity)):    
            l=0
            while l < int(impurity[k]):
                  r = np.random.rand()
                  n = int(np.ceil(r*atoms_num-1))
                  if elements[n] == host:
                      elements[n]=key_dict[k]
                      super_cell.set_chemical_symbols(elements)
                      l=l+1
        super_cell.set_positions(positions)

                  
            

        k0=1
        k1=1
        k2=1
        calc = Vasp(algo='Fast',
                            xc='PBE',
                            encut=500,
                            ediff=0.00001,
                            ediffg=-0.05,
                            isif=2,
                            nsw=300,
                            ibrion=2,
                            lreal='Auto',
                            ivdw=11,
                            ldipol=True,
                            idipol=4,
                            ispin=1,
                            lwave=False,
                            lcharg=False,
                            kpts=(k0,k1,k2),
                            directory='.')
        super_cell.calc = calc
        super_cell.get_potential_energy()
        db.write(super_cell)
        traj.write(super_cell)
        time.sleep(30)
                


            
        
                
# view(atoms)
# print(atoms.cell)
# a = Atoms(atoms)   
# print(a.get_positions())
a = get_random_hcp_alloy(50,  'Cu', Sn=0.166, Cu=0.166, Mn=0.166,Mo=0.166,Fe=0.166, Co=0.17) # run with this command

# print(a)
#get_random_hcp_alloy(2,'Ru', 2, 0.05, Ru=0.5, Fe=0.5)