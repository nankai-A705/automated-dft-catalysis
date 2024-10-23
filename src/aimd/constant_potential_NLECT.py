#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 16:26:09 2024

@author: hxps
"""

from utils import get_GCDFT, get_fmax, get_nelect_neu_pp, get_nelect_incar, get_nelect_neu 
from ase.calculators.vasp import Vasp
from ase.io import read, write
from ase.db import connect
import os
import numpy as np
from ase import units
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.vasp import Vasp
from ase.constraints import FixInternals
import re

def update_incar(incar_file, new_nsw_value,inputs=None):
    """
    Update the NSW value in the INCAR file.
    
    Parameters:
    incar_file (str): Path to the INCAR file.
    new_nsw_value (int): New value for the NSW parameter.
    """
    # Read the contents of the INCAR file
    with open(incar_file, 'r') as file:
        lines = file.readlines()

    # Update the NSW value
    with open(incar_file, 'w') as file:
        for line in lines:
            if inputs in line:
                file.write(f'{inputs} = {new_nsw_value}\n')
            else:
                file.write(line)

def NVT_NoseHoover(atoms=None,
    timestep_fs=1,
    temperature_K=300,
    nsteps = 1000,
    label='',
    pot_target=0.0,
    pot_ref=4.44,
    pot_step=None,
    potentiostat_nsteps=10,
    command=None,
    ):

    # Read geometry if trajectory file exists
    


    nelect_neu = get_nelect_neu()

    print(f'\nFIX-POTENTIAL BOMD (pot_target= {pot_target}, timestep={timestep_fs} fs, potentialstat per {potentiostat_nsteps} fs)')
    if pot_step is None:
        pot_step = 0.0005 * nelect_neu
        print(f'Potentiostat step set to {pot_step:.6f} |e|/V (NELECT_neu={nelect_neu:8.3f} |e|)')

    traj_nelect = []
    traj_pot = []
    traj_gcfe = []
    nsteps_done = 0
    update_incar('INCAR', potentiostat_nsteps,'NSW')
    while nsteps_done < nsteps:
         
        
        
        os.system('mpirun -np 64 /opt/soft/vasp6.3.0/vasp.6.3.0/bin/vasp_gam > vasp.out')
        geom_md = read('CONTCAR')
        os.system('cp CONTCAR POSCAR')
        
        
        print(os.getcwd())
        nsteps_done += potentiostat_nsteps
        new_loop = nsteps_done + potentiostat_nsteps
        os.system(f'ase convert XDATCAR {new_loop}.traj')
        os.system(f'cp OUTCAR OUTCAR{new_loop}')
        
        nelect_net_now, pot_now, energy_gcdft_now = get_GCDFT(pot_ref, dirName=label)
        traj_nelect.append(nelect_net_now)
        traj_pot.append(pot_now) #U_SHE
        traj_gcfe.append(energy_gcdft_now)
        np.savetxt(f'{label}/log_nelect.txt', np.array(traj_nelect))
        np.savetxt(f'{label}/log_potential.txt', np.array(traj_pot))
        np.savetxt(f'{label}/log_gcfe.txt', np.array(traj_gcfe))
        energy = eval([l for l in open('OUTCAR').readlines()
                    if 'free energy' in l][-1].split()[4])
        temperature= eval([l for l in open('OSZICAR').readlines()
                    if 'EK' in l][-1].split()[2])
        print(f'Step {len(traj_nelect):4}: NELECT_net= {traj_nelect[-1]:8.4f} |e|;  U_she= {traj_pot[-1]:8.4f} V; GCFE_el= {traj_gcfe[-1]:12.4f} eV')
        with connect(f'{label}_traj.db', append=True) as db:
            db.write(
                geom_md,
                charge_net=-traj_nelect[-1],
                nelect_net = traj_nelect[-1],
                pot=traj_pot[-1],
                gcfe=traj_gcfe[-1],
                temperature = temperature,
                potene=energy
            )
        if nsteps_done >= nsteps:
            break
        else:
            nelect_net_new = nelect_net_now - (pot_target-pot_now)*pot_step
            with open('INCAR', 'r') as file:
                data = file.read()
            
            # 定义新的NELECT变量值
                new_value = nelect_neu + nelect_net_new
                
                # 使用正则表达式替换NELECT的值
                data = re.sub(r'NELECT = [0-9.]+', f'NELECT = {new_value}', data)
                
                # 将修改后的内容写回INCAR文件
                with open('INCAR', 'w') as file:
                    file.write(data)
           

atoms = read('POSCAR0')
NVT_NoseHoover(atoms,timestep_fs=1,
    temperature_K=300,
    nsteps = 5000,
    label='./',
    pot_target=-0.2,
    pot_ref=4.44,
    pot_step=None,
    potentiostat_nsteps=5,
    command=None,  )    
