#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 10:07:00 2022

@author: hxps
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:59:34 2022

@author: hxps
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 21:03:26 2021
@author: hongye
@email:290720931@qq.com
"""
import numpy as np
from ase.build import bulk
from ase.db import connect
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from ase.constraints import StrainFilter
from ase.io.trajectory import Trajectory
import time


hcp_lattice = {
               'Ru' : [2.733,4.314],
               'Ni' : [2.474,4.070],
               'Fe' : [2.466,3.900]
               }
#bulk_setting = Vasp_default.bulk_settings()
#bulk_setting = bulk_setting['vasp']
#bulk_setting['isif'] = 2


db = connect('RuNi-bulk.db')
            
def get_random_hcp_alloy(num_structure, host, repeat, eps, **kwargs):
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
    png_name = str(file)
    
    #assert(sum(value_dict) == 1)
    for key in file:
        key_dict.append(key)
        value_dict.append(file[key])
    dopant_species = len(key_dict)
    a0 = 0
    c0 = 0
    for i, key in enumerate(key_dict):
        a0 += hcp_lattice[key][0] * value_dict[i]
        c0 += hcp_lattice[key][1] * value_dict[i]
    
    for i in range(num_structure):
        
        seed=1120210312+i*100
        model = bulk(host, 'hcp', a=a0, c=c0) # get cell
        super_cell = model.repeat(repeat)# generate supercell
        positions = super_cell.get_positions()
        cell = super_cell.cell
        angles = cell.angles()
        
        
            
        atoms_ordinal = super_cell.get_atomic_numbers() # get supercell atoms number
        atoms_num = len(atoms_ordinal)
        elements = super_cell.get_chemical_symbols() #[Al, Al, Al ,etc]
        impurity = [] # dict for dopants numbers 
        np.random.seed(seed) # for repeat
        #traj = Trajectory('{}-{}.traj'.format(png_name, i),'w')
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
        k0 = np.ceil(30.0/(repeat*a0))  
        k1 = np.ceil(30.0/(repeat*a0))
        k2 = np.ceil(30.0/(repeat*c0)) 
                  
            
        super_cell.write('POSCAR-{}-{}'.format(png_name,i))
        

        calc = Vasp(algo='Normal',
                            xc='PBE',
                            encut=520,
                            ediff=0.00001,
                            ediffg=-0.02,
                            isif=3,
                            nsw=200,
                            ibrion=2,
                            lreal='Auto',
                            ispin=2,
                            lcharg=False,
                            lwave=False,
                            kpts=(k0,k1,k2),
                            directory='.')
        super_cell.calc = calc
        super_cell.get_potential_energy()
        db.write(super_cell)
        time.sleep(30)
                


            
        
                

       
Ru = [0.125,0.25,0.5,0.75,0.875]
for i in Ru:      
    get_random_hcp_alloy(2,'Ru', 2, 0.1, Ru=i, Ni=1-i) # run with this command
#get_random_hcp_alloy(2,'Ru', 2, 0.05, Ru=0.5, Fe=0.5)
