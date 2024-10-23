#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:21:35 2024

@author: Hongye Qin
@email: 290720931@qq.com
"""



from utils import get_GCDFT, get_fmax, get_nelect_neu_pp, get_nelect_incar, get_nelect_neu 
from ase.calculators.vasp import Vasp
from ase.io import read, write
from ase.db import connect
import os
import subprocess
import numpy as np
from ase import units
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.vasp import Vasp
from ase.constraints import FixInternals
import re



# Define H pseudopotential
max_target_numbers = [1, 1.25, 1.33, 1.5, 1.66, 1.75]
min_target_numbers = [0.25, 0.33, 0.42, 0.5, 0.58, 0.66, 0.75, 1]

def closest_combinations(input_number):
    closest_pair = []
    min_diff = float('inf')

    if input_number > 0:
        # 从 max_target_numbers 中减去 1
        adjusted_numbers = [num - 1 for num in max_target_numbers]
        target_sum = input_number
    else:
        # 从 min_target_numbers 中减去 1
        adjusted_numbers = [num - 1 for num in min_target_numbers]
        target_sum = input_number  # 目标和直接为输入值

    # 遍历所有可能的组合
    for i in range(len(adjusted_numbers)):
        for j in range(i, len(adjusted_numbers)):
            num1 = adjusted_numbers[i]
            num2 = adjusted_numbers[j]
            current_sum = num1 + num2
            
            # 确保当前和与目标和的差最小
            current_diff = abs(current_sum - target_sum)
            if current_diff < min_diff:
                min_diff = current_diff
                # 返回原始值为列表
                if input_number > 0:
                    closest_pair = [max_target_numbers[i], max_target_numbers[j]]
                else:
                    closest_pair = [min_target_numbers[i], min_target_numbers[j]]

    return closest_pair

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
        nelect_net_new = pot_target-pot_now

        with connect(f'{label}_traj.db', append=True) as db:
            db.write(
                geom_md,
                charge_net=-traj_nelect[-1],
                nelect_net = traj_nelect[-1],
                pot=traj_pot[-1],
                gcfe=traj_gcfe[-1],
                temperature = temperature,
                potene=energy,
                H_net = nelect_net_new
            )
        if nsteps_done >= nsteps:
            break
        else:
          #  os.remove('POTCAR')
            nelect_net_new = pot_now-pot_target
            print(nelect_net_new)
  
            if  -1.5 <= nelect_net_new <= 1.5:
                psepotentional_comination = closest_combinations(nelect_net_new)
                np.savetxt(f'{label}/log_pseH.txt', np.array(psepotentional_comination))
                for key,value in enumerate(psepotentional_comination):
                    filename1 = f'POTCAR_H{psepotentional_comination[0]}'
                    filename2 = f'POTCAR_H{psepotentional_comination[1]}'
                    
                    # 执行 cat 命令
                    try:
                        command = f'cat POTCAR_Ru POTCAR_O POTCAR_H1 {filename1} {filename2}> POTCAR'
                        
                        subprocess.run(command, shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"执行命令失败: {e},please check the code")



           

atoms = read('POSCAR')
NVT_NoseHoover(atoms,timestep_fs=1,
    temperature_K=300,
    nsteps = 5000,
    label='./',
    pot_target=-0.2,
    pot_ref=4.44,
    pot_step=None,
    potentiostat_nsteps=3,
    command=None,  )    
