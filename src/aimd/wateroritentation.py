# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 22:13:01 2024

@author: Hongye Qin
@email: hongyechin@qq.com
"""
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import subprocess
import os
from post_processing_for_aimd1 import get_water_according_to_atom, assert_distance, get_water_according_to_space
from ase.visualize import view

class WaterOrientationalRelaxation():
    
    """
    This module forked from 
    https://docs.mdanalysis.org/1.0.1/_modules/MDAnalysis/analysis/waterdynamics.html#WaterOrientationalRelaxation
    """
    
    def __init__(self, file, dt, start_frame=None, end_frame=None):
        
        if file.endswith(".traj"):
            
             self.file = Trajectory(file)
             
        else:
            
            print('Please input AIMD Traj file')
        self.dt = dt
        
        
        if start_frame is None:
           self.start_frame = 1
           self.ref_id = 1
           # we define the first frame as water oritentaion reference
        else:
            self.start_frame = start_frame
            self.ref_id = start_frame
            
        if end_frame is None:
            
            self.end_frame = len(self.file)
        
        else:
            self.end_frame = end_frame
            
        pass
    
    def get_objetctive_water(self):
        
        
        
        
        pass
    def get_ref_Cvect(self, min_posposition_z, max_position_z, atom_symbol):
        
        atoms = self.file[self.ref_id]
        
        h, o = get_water_according_to_space(atoms, min_posposition_z,max_position_z, atom_symbol)
        water_group, distance_group, angle_list, position_group = assert_distance(atoms, o, h, 1.2)
        self.ref_bisetion_vect = np.zeros([len(angle_list),3])
        self.ref_water_group = water_group
        self.min_posposition_z = min_posposition_z
        self.max_position_z = max_position_z
        # self.ref_bisetion_vect = {}
        self.ref_O = []
        self.ref_H = []
        for index,(k, v) in enumerate(water_group.items()):
            self.ref_O.append(k)
            self.ref_H.extend(v)
            O_pos = atoms.get_positions()[k]
            H1_index = water_group[k][0]
            H2_index = water_group[k][1]
            H1_pos = atoms.get_positions()[H1_index]
            H2_pos = atoms.get_positions()[H2_index]
            vec1 = (H1_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            vec2 = (H2_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            normal_water_bisetion_vector = (vec1 + vec2)/np.linalg.norm(vec1 + vec2)
            self.ref_bisetion_vect[index] =  normal_water_bisetion_vector

    def get_Cvect(self, atoms):
        water_group, distance_group, angle_list, position_group = assert_distance(atoms, self.ref_O, self.ref_H, 1.2)
        
        water_bisetion_vector = np.zeros([len(water_group),3])
        for index, (k, v) in enumerate(water_group.items()):

            O_pos = atoms.get_positions()[k]
            H1_index = water_group[k][0]
            H2_index = water_group[k][1]
            H1_pos = atoms.get_positions()[H1_index]
            H2_pos = atoms.get_positions()[H2_index]
            vec1 = (H1_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            vec2 = (H2_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            bisetion_vector = (vec1 + vec2)/np.linalg.norm(vec1 + vec2)
            water_bisetion_vector[index] = bisetion_vector
        
        
        val_bisection  = self.lg2(np.sum(water_bisetion_vector*self.ref_bisetion_vect, axis=1))
        # print(val_bisection.shape)
        # print(val_bisection)
        
        
    
        return val_bisection
    
    def run(self, mode=None):


            
        dt = [t for t in range(self.start_frame+1, self.end_frame, self.dt)] 
        p2 = []
        

        t = []   
        for j in dt:
            n=0
            val = np.zeros(len(self.ref_water_group))
            for i in range(self.start_frame,int(j), self.dt):
                atoms = self.file[i]
                val += self.get_Cvect(atoms)
                n += 1
            print(n*self.dt)
            val = val / n
            p2.append(np.average(val))
            t.append(j/1000*0.5)
            
        with open("output.txt", "w") as file:
            # 将时间和数据同时写入文件，每行一个时间和对应的数据
            for time, value in zip(t, p2):
                file.write(f"{time} {value}\n")   
        print(p2)
        return t , p2 
        
    @staticmethod
    def lg2(x):
        """Second Legendre polynomial"""
        return (3*x*x - 1)/2
    
    def plot(self):
        dt2, p22 = self.run()
        
        fig = plt.figure( figsize = ( 3.5, 3.5 ) )
        # ax1 = fig.add_subplot( 2, 1, 1 )
        # ax1.plot( dt1, p21, label = '1 distribution' )
        # ax1.set_xlabel( 'time (A)' )
        # ax1.set_ylabel( 'P2 (a.u.)' )
        ax1 = fig.add_subplot( 2, 1, 2 )
        ax1.plot( dt2, p22, label = '2 distribution' )
        ax1.set_xlabel( 'time (A)' )
        ax1.set_ylabel( 'P2 (a.u.)' )


        
        
    pass
# a = WaterOrientationalRelaxation('RuNi-1k.traj', 100,1000,5000)
# a.get_ref_Cvect(3.5, 'Ru', 2)

# dt1, p21 = a.run()
# print(a.ref_water_group)

b = WaterOrientationalRelaxation('Ru-1k.traj', 250,2000,12000)
b.get_ref_Cvect(0,3.5, 'Ru',)
print(b.ref_bisetion_vect)
b.plot()
