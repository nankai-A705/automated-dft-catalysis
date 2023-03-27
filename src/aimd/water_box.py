# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 10:19:11 2022

@author: win
"""

from ase import Atoms
from ase.constraints import FixBondLengths
from ase.calculators.tip3p import TIP3P, rOH, angleHOH
from ase.calculators.tip4p import TIP4P
from ase.md import Langevin
import ase.units as units
from ase.io.trajectory import Trajectory
import numpy as np
from ase.build import molecule





def water_tip3_box(target_a, number, offset):
    """
    Get water box phi(a,b) not equal to 90 degree.
    
    target_a : lattice a (slab)
    number: water molecule in box
    offset： tune this parameters to wrap molecule in box
    
    """
    x = angleHOH * np.pi / 180 / 2

    pos = [[0, 0, 0],
       [0, rOH * np.cos(x), rOH * np.sin(x)],
       [0, rOH * np.cos(x), -rOH * np.sin(x)]]
    atoms = Atoms('OH2', positions=pos)
    total_volume = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))
    target_volume = number * total_volume
    target_c = target_volume / (target_a**2 * np.sqrt(3) /2 )
    a = ((18.01528 / 6.022140857e23) / (1 / 1e24))**(1 / 3.) # 1g/ml 25 degree
    c = total_volume / (a**2 * np.sqrt(3) /2 )
    atoms.set_cell([a,a,c, 90,90,60])
    atoms.center()
    x_n = int(target_a / a)
    z_n = int(number / (x_n**2))
    # atoms.repeat((x_n,x_n,z_n))
    # atoms.set_pbc(True)
    Ga = np.array([target_a, 0, 0])
    Gb = np.array([target_a/2, (target_a * np.sqrt(3))/2, 0])
    Gc = np.array([0, 0, target_c])
    n1 = x_n
    n2 = x_n
    n3 = z_n
    assert n1*n2*n3 == number
    e1 = Ga/n1
    e2 = Gb/n2
    e3 = Gc/n3
    pos = np.zeros((n1,n2,n3,3))
    lpos = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                pos[i][j][k][:] = i*e1 + j*e2 + k*e3
                lpos.append(pos[i][j][k][:])
    
    position = [p.tolist() for p in lpos]
    box = Atoms()
    box.set_cell([target_a,target_a, target_c, 90,90,60])
    for pos in position:
        newp = [[pos[0]+offset,pos[1]+offset,pos[2]+offset],
               [pos[0]+offset, pos[1]+rOH * np.cos(x)+offset, pos[2]+rOH * np.sin(x)+offset],
               [pos[0]+offset, pos[1]+rOH * np.cos(x)+offset, pos[2]-rOH * np.sin(x)+offset]]
        water =  Atoms('OH2', positions=newp)
        
        box.extend(water)
        box.set_pbc(True)
        
    
    return box

atoms = water_tip3_box(10.91, 27,1)
atoms.write('a.cif')

    
from ase.visualize import view


# RATTLE-type constraints on O-H1, O-H2, H1-H2.
# atoms.constraints = FixBondLengths([(3 * i + j, 3 * i + (j + 1) % 3)
#                                     for i in range(3**3)
#                                     for j in [0, 1, 2]])

# tag = 'tip3p_27mol_equil'
# atoms.calc = TIP4P(rc=4.5)
# md = Langevin(atoms, 1 * units.fs, temperature=300 * units.kB,
#               friction=0.01, logfile=tag + '.log')

# traj = Trajectory(tag + '.traj', 'w', atoms)
# md.attach(traj.write, interval=1)
# md.run(4000)