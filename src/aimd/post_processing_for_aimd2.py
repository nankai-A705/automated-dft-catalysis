# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 20:02:46 2023

@author: Hongye Qin
@email: 290720931@qq.com
"""

import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
import numpy as np
import pandas as pd
import subprocess
import os

from ase.visualize import view


# traj = Trajectory( 'Ru.traj', 'r' )
# total_frame = len(traj)
# tag_group = []
# first_frame = traj[2551]

#view(first_frame)
# print(total_frame)
class aimd_post_process:
    def __init__(self, file=None, start_frame=None):
        
        if file is None:
            if os.path.exists('XDATCAR'):
               subprocess.run(["ase", "convert", "XDATCAR", "aimd.traj"])
               file = 'aimd.traj'
        
        elif file.endswith(".traj"):
            
             file = Trajectory(file)
        
        else:
            raise ValueError('please input right file format, This version only support XDATCAR and traj format')
        
        self.start_frame = 0
        self.dpoints = 61 # define space interval point toward Z vectors
        self.file = file
        self.OH_threshold = 1.2 # define O-H bonds length less than 1.2 belongs to one water molecule
    def get_water_according_to_distance(self, atoms, thickness, atom_symbol):
        
        top_metal_index = atoms[[ at.index for at in atoms if at.symbol == atom_symbol]]
        maxz_atom_symbol = np.max( top_metal_index.positions[ :, 2 ] )
        hydrogen_index = [i.index for i in atoms if i.symbol == 'H']
        oxygen_index = [i.index for i in atoms if i.symbol == 'O' and i.position[2] - maxz_atom_symbol <= thickness]
        
        if len(oxygen_index) == 0:
            
            raise ValueError('You should choose right frame for analysis or change the thickness of edl')
            
            
        # else:
            
        #     water_group = {}
        #     for oxygen in oxygen_index:
        #         single_water_group = []
        #         single_hydrogen_group = []
        #         distances = atoms.get_distances(oxygen, hydrogen_index)
        #         bottom_k_index = ArrayBottomK(2, distances)
        #         for i in bottom_k_index:
        #             single_hydrogen_group.append(hydrogen_index[i])
                    
           
        #     #single_hydrogen_group.append()
        #         water_group[oxygen] = single_hydrogen_group #{O index: H1 index, H2 index}
            
        return  hydrogen_index, oxygen_index
    
    def assert_distance(self,atom, oxygen_index, hydrogen_index):
        
        """
         threshold :  O-H bond length to make sure belong to the water molecule
         
        """
        
        cell = atom.get_cell()
        water_group = {}
        position_group = {}
        distance_group = {}
        total_angle = []
        for i in oxygen_index:
            hydrogen = []
            distance = []
            H_position = []
            O_pos = atom.get_positions()[i]
            
            for j in hydrogen_index:
                if atom.get_distance(i,j) < self.OH_threshold:
                    hydrogen.append(j)
                    distance.append(atom.get_distance(i, j))
                    H_position.append(atom.get_positions()[j])
                    
                
                    
                elif atom.get_distance(i,j) >= self.OH_threshold:
                    
                    
                    new_maxxpos = atom.get_positions()[j] + cell[0]
                    dis = np.linalg.norm(O_pos - new_maxxpos)
                    if dis < self.OH_threshold:
                       hydrogen.append(j)
                       atom[j].position = new_maxxpos
                       distance.append(dis)
                       H_position.append(new_maxxpos)
                       
                    else:
                        new_minxpos = atom.get_positions()[j]  - cell[0]
                        
                        new_dis = np.linalg.norm(O_pos - new_minxpos)
                        if new_dis < self.OH_threshold:
                            hydrogen.append(j)
                            atom[j].position = new_minxpos
                            distance.append(new_dis)
                            H_position.append(new_minxpos)
                            
                            
                        
                        else: 
                            new_maxypos = atom.get_positions()[j]  + cell[1]
                            new_maxydis = np.linalg.norm(O_pos - new_maxypos)
                            
                            if new_maxydis < self.OH_threshold:
                                hydrogen.append(j)
                                atom[j].position = new_maxypos
                                distance.append(new_maxydis)
                                H_position.append(new_maxypos)
                                
                            else:
                                new_minypos = atom.get_positions()[j]  - cell[1]
                                new_minydis = np.linalg.norm(O_pos - new_minypos)
                                if new_minydis < self.OH_threshold:
                                    hydrogen.append(j)
                                    atom[j].position = new_minxpos
                                    distance.append(new_minydis)
                                    H_position.append(new_minypos)
                                
                                else:
                                    new_maxxypos = atom.get_positions()[j]  + cell[0] + cell[1]
                                    new_maxxydis = np.linalg.norm(O_pos - new_maxxypos)
                                    if new_maxxydis < self.OH_threshold:
                                        hydrogen.append(j)
                                        atom[j].position = new_maxxypos
                                        distance.append(new_maxxydis)
                                        H_position.append(new_maxxypos)
                                    
                                    else:
                                         new_minxypos = atom.get_positions()[j]  - cell[0] - cell[1]
                                         new_minxydis = np.linalg.norm(O_pos - new_minxypos)
                                         if new_minxydis < self.OH_threshold:
                                             hydrogen.append(j)
                                             atom[j].position = new_minxypos
                                             distance.append(new_minxydis)
                                             H_position.append(new_minxypos)
                                             
                                         else:
                                             new_max_xypos = atom.get_positions()[j]  - cell[0] + cell[1]
                                             new_max_xypos = np.linalg.norm(O_pos - new_max_xypos)
                                             if new_max_xypos < self.OH_threshold:
                                                 hydrogen.append(j)
                                                 atom[j].position = new_max_xypos
                                                 distance.append(new_max_xypos)
                                                 H_position.append(new_max_xypos)
                                             else:
                                                 new_min_xypos = atom.get_positions()[j]  + cell[0] - cell[1]
                                                 new_min_xypos = np.linalg.norm(O_pos - new_min_xypos)
                                                 if new_min_xypos < self.OH_threshold:
                                                     hydrogen.append(j)
                                                     atom[j].position = new_min_xypos
                                                     distance.append(new_min_xypos)
                                                     H_position.append(new_min_xypos)
                                                 
                                         
                                         
                                    
                                
                                
                                
            if len(hydrogen) != 2:
                print(hydrogen)
            water_group[i] = hydrogen
            distance_group[i] = distance
            position_group[i] = H_position
            angle_list = []
            for k,v in water_group.items():
                
                O_position = atom.get_positions()[k]
                vec1 = v[0] - O_position
                vec2 = v[1] - O_position
                ag = get_angle(vec1,vec2)
                angle_list.append(ag)
            
            
            
                                
                            
                            
                            
        
        
        return water_group, distance_group, angle_list, position_group
        
        
        
        
 
        
        
        
             
                 
       
    
    
    
def return_error_water_group(atoms, groups):
    oxygen_index = []
    for k, v in groups.items():
        if abs(v[0]-v[1]) >= 1:
            oxygen_index.append(k)
            
        
        
    return oxygen_index

def get_angle(vec1, vec2):
    
    data1 = np.sqrt(np.sum(vec1*vec1))
    data2 = np.sqrt(np.sum(vec2*vec2))
    cos_theta = np.sum(vec1*vec2)/(data1*data2)
    
    
    angle = np.degrees(np.arccos(cos_theta))
    return angle

def assert_distance(atom, oxygen_index, hydrogen_index,threshold):
    
    """
     threshold :  O-H bond length to make sure belong to the water molecule
     
    """
    
    cell = atom.get_cell()
    water_group = {}
    position_group = {}
    distance_group = {}
    
    for i in oxygen_index:
        hydrogen = []
        distance = []
        H_position = []
        O_pos = atom.get_positions()[i]
        
        for j in hydrogen_index:
            if atom.get_distance(i,j) < threshold:
                hydrogen.append(j)
                distance.append(atom.get_distance(i, j))
                H_position.append(atom.get_positions()[j])
                
            
                
            elif atom.get_distance(i,j) >= threshold:
                
                
                new_maxxpos = atom.get_positions()[j] + cell[0]
                dis = np.linalg.norm(O_pos - new_maxxpos)
                if dis < threshold:
                   hydrogen.append(j)
                   atom[j].position = new_maxxpos
                   distance.append(dis)
                   H_position.append(new_maxxpos)
                   
                else:
                    new_minxpos = atom.get_positions()[j]  - cell[0]
                    
                    new_dis = np.linalg.norm(O_pos - new_minxpos)
                    if new_dis < threshold:
                        hydrogen.append(j)
                        atom[j].position = new_minxpos
                        distance.append(new_dis)
                        H_position.append(new_minxpos)
                        
                        
                    
                    else: 
                        new_maxypos = atom.get_positions()[j]  + cell[1]
                        new_maxydis = np.linalg.norm(O_pos - new_maxypos)
                        
                        if new_maxydis < threshold:
                            hydrogen.append(j)
                            atom[j].position = new_maxypos
                            distance.append(new_maxydis)
                            H_position.append(new_maxypos)
                            
                        else:
                            new_minypos = atom.get_positions()[j]  - cell[1]
                            new_minydis = np.linalg.norm(O_pos - new_minypos)
                            if new_minydis < threshold:
                                hydrogen.append(j)
                                atom[j].position = new_minypos
                                distance.append(new_minydis)
                                H_position.append(new_minypos)
                            
                            else:
                                new_maxxypos = atom.get_positions()[j]  + cell[0] + cell[1]
                                new_maxxydis = np.linalg.norm(O_pos - new_maxxypos)
                                if new_maxxydis < threshold:
                                    hydrogen.append(j)
                                    atom[j].position = new_maxxypos
                                    distance.append(new_maxxydis)
                                    H_position.append(new_maxxypos)
                                
                                else:
                                     new_minxypos = atom.get_positions()[j]  - cell[0] - cell[1]
                                     new_minxydis = np.linalg.norm(O_pos - new_minxypos)
                                     if new_minxydis < threshold:
                                         hydrogen.append(j)
                                         atom[j].position = new_minxypos
                                         distance.append(new_minxydis)
                                         H_position.append(new_minxypos)
                                         
                                     else:
                                          # print(j)
                                          
                                           new_max_xypos = atom.get_positions()[j]  - cell[0] + cell[1]
                                           new_max_xydis = np.linalg.norm(O_pos - new_max_xypos)
                                           if new_max_xydis < threshold:
                                               # print(atom[j].position)
                                               hydrogen.append(j)
                                               atom[j].position = new_max_xypos
                                               distance.append(new_max_xydis)
                                               H_position.append(new_max_xypos)
                                               # print(j)
                                               # print(atom[j].position)
                                           else:
                                               new_min_xypos = atom.get_positions()[j]  + cell[0] - cell[1]
                                               new_min_xydis = np.linalg.norm(O_pos - new_min_xypos)
                                               if new_min_xydis < threshold:
                                                   hydrogen.append(j)
                                                   atom[j].position = new_min_xypos
                                                   distance.append(new_min_xydis)
                                                   H_position.append(new_min_xypos)
                                             
                                     
                                     
                                
                            
                            
                            
        if len(hydrogen) != 2:
            print(hydrogen)
        water_group[i] = hydrogen
        distance_group[i] = distance
        position_group[i] = H_position
        
        total_angle = []
        
        for k,v in water_group.items():
            
            O_position = atom.get_positions()[k]
            vec1 = atom.get_positions()[v[0]] - O_position
            vec2 = atom.get_positions()[v[1]] - O_position
            ag = get_angle(vec1,vec2)
            total_angle.append(ag)
        
        
        
                            
                        
                        
                        
    
    
    return water_group, distance_group, total_angle, position_group

           
def ArrayTopK(top_k,arr):
    top_k_idx=arr.argsort()[::-1][0:top_k]
    return top_k_idx

def ArrayBottomK(top_k,arr):
    bottom_k_idx=arr.argsort()[0:top_k]
    return bottom_k_idx      
    
# a = np.array([10,2,1,5])
# b = ArrayBottomK(1, a)
# print(b)

def get_water_according_to_distance(atoms, thickness, atom_symbol):
    
    top_metal_index = atoms[[ at.index for at in atoms if at.symbol == atom_symbol ]]
    maxz_atom_symbol = np.max( top_metal_index.positions[ :, 2 ] )
    hydrogen_index = [i.index for i in atoms if i.symbol == 'H']
    oxygen_index = [i.index for i in atoms if i.symbol == 'O' and i.position[2] - maxz_atom_symbol <= thickness]
    
    if len(oxygen_index) == 0:
        
        raise ValueError('You should choose right frame for analysis or change the thickness of edl')
        
        
    # else:
        
    #     water_group = {}
    #     for oxygen in oxygen_index:
    #         single_water_group = []
    #         single_hydrogen_group = []
    #         distances = atoms.get_distances(oxygen, hydrogen_index)
    #         bottom_k_index = ArrayBottomK(2, distances)
    #         for i in bottom_k_index:
    #             single_hydrogen_group.append(hydrogen_index[i])
                
       
    #     #single_hydrogen_group.append()
    #         water_group[oxygen] = single_hydrogen_group #{O index: H1 index, H2 index}
        
    return  hydrogen_index, oxygen_index


def get_water_according_to_space(atoms, min_pos, max_pos, atom_symbol):
    
    top_metal_index = atoms[[ at.index for at in atoms if at.symbol == atom_symbol ]]
    maxz_atom_symbol = np.max( top_metal_index.positions[ :, 2 ] )
    hydrogen_index = [i.index for i in atoms if i.symbol == 'H'] # get all H index
    oxygen_index = [i.index for i in atoms if i.symbol == 'O' and min_pos <= i.position[2] - maxz_atom_symbol <= max_pos]
    
    # if len(oxygen_index) == 0:
        
    #     raise ValueError('You should choose right frame for analysis or change the thickness of edl')
        
        
    # else:
        
    #     water_group = {}
    #     for oxygen in oxygen_index:
    #         single_water_group = []
    #         single_hydrogen_group = []
    #         distances = atoms.get_distances(oxygen, hydrogen_index)
    #         bottom_k_index = ArrayBottomK(2, distances)
    #         for i in bottom_k_index:
    #             single_hydrogen_group.append(hydrogen_index[i])
                
       
    #     #single_hydrogen_group.append()
    #         water_group[oxygen] = single_hydrogen_group #{O index: H1 index, H2 index}
        
    return  hydrogen_index, oxygen_index

def get_water_according_to_atom(atoms, position_z, atom_symbol, number):
    
    specific_atoms_index = [atom.index for atom in atoms if atom.symbol == atom_symbol and atom.position[2] >= position_z]
    hydrogen_index = [i.index for i in atoms if i.symbol == 'H']
    oxygen_index = [i.index for i in atoms if i.symbol == 'O']
    nearest_oxygen_index = []
    for i in specific_atoms_index:
        distances = atoms.get_distances(i, oxygen_index, mic=True)
        near_id = ArrayBottomK(number, distances)
        for j in range(len(near_id)):
            nearest_oxygen_index.append(oxygen_index[near_id[j]])
    
    return hydrogen_index, nearest_oxygen_index
######
"""test"""
#######
# traj = Trajectory('RuNi.traj', 'r' )
# atoms = traj[10550]
# # view(atoms)
# # print(atoms.get_cell())
# # cell = atoms.get_cell()
# # print(atoms.get_positions()[127])

# h,o = get_water_according_to_atom(atoms,7, 'Ni',  2 )
# print(o)

# water_group, distance_group, angle_list, position = assert_distance(atoms, o, h, 1.3)
# # print(atoms.get_positions()[127])
# print(water_group)
# print(angle_list)
# atoms.write('test.cif')


def get_nearest_Oxygen(atom, thickness,atom_symbol, mthreshold,threshold):
    """
    Return the number of oxygen atoms near the surface according to thikness
         。。。     max_z
           |       thikness
    。。。。。。。。surface
    。。。。。。。。
    threshold: if the distance of nearest oxygen atoms large than that of threshold, then the oxygen will not be consider.
    """
    
    cell = atom.get_cell()
    
    hydogen_list, oxygen_list = get_water_according_to_distance(atom, thickness, atom_symbol)
    total_oxygen_list = [i.index for i in atom if i.symbol == 'O' ]
    primitive = {}
    primitive_maxx = {}
    primitive_minx = {}
    primitive_maxy = {}
    primitive_miny = {}
    primitive_maxxy = {}
    primitive_minxy = {}
    
    for i in oxygen_list:
        distance = []
        near_distance_index = []
        for j in total_oxygen_list:
            distance = atom.get_distances(i,j)
            if distance != 0 and mthreshold < distance <= threshold:
                near_distance_index.append(j)
                
                
        primitive[i] = near_distance_index
    
    for i in oxygen_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in total_oxygen_list:
            new_maxxpos = atom.get_positions()[j]  + cell[0]
            distance = np.linalg.norm(O_pos - new_maxxpos)
            if distance != 0 and mthreshold < distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_maxx[i] = near_distance_index
        
    for i in oxygen_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in total_oxygen_list:
            new_minxpos = atom.get_positions()[j]  - cell[0]
            distance = np.linalg.norm(O_pos - new_minxpos)
            if distance != 0 and mthreshold < distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_minx[i] = near_distance_index

    for i in oxygen_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in total_oxygen_list:
            new_maxypos = atom.get_positions()[j]  + cell[1]
            distance = np.linalg.norm(O_pos - new_maxypos)
            if distance != 0 and mthreshold < distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_maxy[i] = near_distance_index
        
    for i in oxygen_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in total_oxygen_list:
            new_minypos = atom.get_positions()[j]  - cell[1]
            distance = np.linalg.norm(O_pos - new_minypos)
            if distance != 0 and mthreshold < distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_miny[i] = near_distance_index

    for i in oxygen_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in total_oxygen_list:
            new_maxxypos = atom.get_positions()[j]  + cell[1] + cell[0]
            distance = np.linalg.norm(O_pos - new_maxxypos)
            if distance != 0 and mthreshold < distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_maxxy[i] = near_distance_index

    for i in oxygen_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in total_oxygen_list:
            new_minxypos = atom.get_positions()[j]  - cell[1] - cell[0]
            distance = np.linalg.norm(O_pos - new_minxypos)
            if distance != 0 and mthreshold < distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_minxy[i] = near_distance_index
            
    links = {key:primitive[key]+primitive_minx[key]+primitive_maxy[key]+primitive_miny[key]+primitive_maxx[key]+primitive_maxxy[key]+primitive_minxy[key] for key in primitive} 
    # return primitive, primitive_minx, primitive_maxx,  primitive_maxy, primitive_miny
    return links

def get_global_nearest_Oxygen(atom, pair_A, pair_B, mthreshold,threshold):
    """
    Return the number of oxygen atoms near the surface according to thikness
         。。。     max_z
           |       thikness
    。。。。。。。。surface
    。。。。。。。。
    threshold: if the distance of nearest oxygen atoms large than that of threshold, then the oxygen will not be consider.
    """
    
    cell = atom.get_cell()
    
    B_atom_list = [i.index for i in atom if i.symbol == pair_B]
    A_atom_list = [i.index for i in atom if i.symbol == pair_A ]
    primitive = {}
    primitive_maxx = {}
    primitive_minx = {}
    primitive_maxy = {}
    primitive_miny = {}
    primitive_maxxy = {}
    primitive_minxy = {}
    
    for i in A_atom_list:
        distance = []
        near_distance_index = []
        for j in B_atom_list:
            distance = atom.get_distances(i,j)
            if distance != 0 and mthreshold < distance and distance <= threshold:
                near_distance_index.append(j)
                
                
        primitive[i] = near_distance_index
    
    for i in A_atom_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in B_atom_list:
            new_maxxpos = atom.get_positions()[j]  + cell[0]
            distance = np.linalg.norm(O_pos - new_maxxpos)
            if distance != 0 and mthreshold < distance and distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_maxx[i] = near_distance_index
        
    for i in A_atom_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in B_atom_list:
            new_minxpos = atom.get_positions()[j]  - cell[0]
            distance = np.linalg.norm(O_pos - new_minxpos)
            if distance != 0 and mthreshold < distance and distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_minx[i] = near_distance_index

    for i in A_atom_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in B_atom_list:
            new_maxypos = atom.get_positions()[j]  + cell[1]
            distance = np.linalg.norm(O_pos - new_maxypos)
            if distance != 0 and mthreshold < distance and distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_maxy[i] = near_distance_index
        
    for i in A_atom_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in B_atom_list:
            new_minypos = atom.get_positions()[j]  - cell[1]
            distance = np.linalg.norm(O_pos - new_minypos)
            if distance != 0 and mthreshold < distance and distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_miny[i] = near_distance_index

    for i in A_atom_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in B_atom_list:
            new_maxxypos = atom.get_positions()[j]  + cell[1] + cell[0]
            distance = np.linalg.norm(O_pos - new_maxxypos)
            if distance != 0 and mthreshold < distance and distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_maxxy[i] = near_distance_index

    for i in A_atom_list:
        distance = []
        near_distance_index = []
        O_pos = atom.get_positions()[i]
        for j in B_atom_list:
            new_minxypos = atom.get_positions()[j]  - cell[1] - cell[0]
            distance = np.linalg.norm(O_pos - new_minxypos)
            if distance != 0 and mthreshold < distance and distance <= threshold:
                near_distance_index.append(j)
                
        
        primitive_minxy[i] = near_distance_index
            
    links = {key:primitive[key]+primitive_minx[key]+primitive_maxy[key]+primitive_miny[key]+primitive_maxx[key]+primitive_maxxy[key]+primitive_minxy[key] for key in primitive} 
    # return primitive, primitive_minx, primitive_maxx,  primitive_maxy, primitive_miny
    return links
    
### test the number of H bonds     
# O_links = get_nearest_Oxygen(atoms, 5.0, 'Ru',  4.9)
# oxygen_index = [key for key,v in O_links.items()]
# h, o = get_water_according_to_distance(atoms, 5.0, 'Ru')
# water_group, b, c, d = assert_distance(atoms, o, h, 1.2)
# # print(water_group)
# print(O_links)
# degree = []
# for k, v in O_links.items():
#     deg1 = []
#     deg2 = []
#     count = 0
#     for i in v:
#         O_accepter = atoms.get_positions()[k]
#         O_doner = atoms.get_positions()[i]
        
#         H1_index = water_group[i][0]
#         H2_index = water_group[i][1]
#         H1 = atoms.get_positions()[H1_index]
#         H2 = atoms.get_positions()[H2_index]
        
#         angle1 = get_angle(O_accepter-O_doner, H1-O_doner)
#         angle2 = get_angle(O_accepter-O_doner, H2-O_doner)
#         dis1 = np.linalg.norm(O_accepter - H1)
#         dis2 = np.linalg.norm(O_accepter - H2)
        
#         if angle1 < 35 and dis1 <=3.5:
#             deg1.append(angle1)
        
#         if angle2 < 35 and dis2 <=3.5:
#             deg2.append(angle1)
            
#     counts = len(deg1+deg2)
#     print(deg1+deg2)
    
#     degree.append(counts)
# print(degree)

        
    
# print(b) 
# print(c)  
# print(d)  
# print(e)  
# f = {key:a[key]+b[key]+c[key]+d[key]+e[key] for key in a}   
# print(f)  
    
# oxygen_index = get_water_according_to_distance(first_frame, 3.5, 'Ru')[2]
# hydrogen_index = get_water_according_to_distance(first_frame, 3.5, 'Ru')[1]
# water_group = assert_distance(first_frame, oxygen_index, hydrogen_index, 1.3)



def aimd_plot(traj, min_pos, max_pos, threshold, start_frame, end_frame=None):
    
    
    traj = Trajectory(traj, 'r' )
    if end_frame is None:
        
       total_frame = len(traj)
    else:
        total_frame = end_frame
    total_angle = []
    total_distance = []
    for i in range(start_frame,total_frame):
        img = traj[i]
        img.set_pbc = [True]

    

        distance = []
        h,o = get_water_according_to_space(img, min_pos, max_pos, 'Ru')
        water_group, distance_group, angle_list, position_group = assert_distance(img, o, h, threshold)
    
    
        for k, v in distance_group.items():
           distance.append(v[0])

        
        mean_angle = np.mean(angle_list)
        mean_distance = np.mean(distance)
        total_angle.append(mean_angle)
        total_distance.append(mean_distance)
    

    angle = np.array(total_angle)
    distance = np.array(total_distance)
    time = np.arange(start_frame, total_frame )
    fig = plt.figure( figsize = ( 5, 5 ) )
    ax1 = fig.add_subplot( 2, 1, 1 )
    ax1.plot( time, angle, label = 'water angle' )
    ax1.plot(time, [np.mean(angle)]*len(time))
    ax1.set_xlabel( 'Time/ps' )
    ax1.set_ylabel( 'angle [K]' )
    # ax1.annotate(text='mean angle', xy=(time/2, np.mean(angle)), xytext=(time/2, np.mean(angle)+0.03))
    ax1.set_ylim(np.min(angle)-0.5, np.max(angle)+0.5)

    ax2 = fig.add_subplot( 2, 1, 2 )
    ax2.plot( time, distance, label = 'OH distance' )
    ax2.plot(time, [np.mean(distance)]*len(time))
    ax2.set_xlabel( 'Time/ps' )
    ax2.set_ylabel( 'distance [A]' )
    ax2.set_ylim(np.min(distance)-0.05, np.max(distance)+0.05)

# a = aimd_plot('Ru.traj', 0, 3.5, 1.2, 1)
def get_water_index(atoms, cdl, symbol):

    h, o = get_water_according_to_distance(atoms, cdl, symbol)
    water_group, b, c, d = assert_distance(atoms, o, h, 1.2)
    index = []

    for k,v in water_group.items():
        index.append(k)
        index.append(v[0])
        index.append(v[1])
    
    index = [i+1 for i in index]
        
    return index

def get_H_bonds_according_to_index(atoms, index, circle):
    
    H_index  =  [atom.index for atom in atoms if atom.symbol == 'H' ]
    O_index = [atom.index for atom in atoms if atom.symbol == 'O' ]
    total_index = H_index + O_index 
    nearest_atom_index = []
    water_H_index = []
    nearest_atom_distance = atoms.get_distances(index, total_index, mic=True)
    for  i  in range(len(nearest_atom_distance)):
        
        if 0.1 <= nearest_atom_distance[i]  <= circle:
            
            nearest_atom_index.append(total_index[i])
    for  i  in range(len(nearest_atom_distance)):
        
        if 0.1 <= nearest_atom_distance[i]  <= 1.2:
            
            water_H_index.append(total_index[i])
    assert len(water_H_index) == 2
            
    nearest_atom_H_index = []
    nearest_atom_O_index = []
    for j in nearest_atom_index:
        if j in H_index and j not in water_H_index:
            nearest_atom_H_index.append(j)
        
        elif j in O_index:
            nearest_atom_O_index.append(j)
    
    
     

    
    return nearest_atom_H_index, nearest_atom_O_index, water_H_index

def get_h_bonds_signle_water(atoms, index, circle):
    
    """
    Return the number of H bonds near the surface
    
    """
    H_links = get_H_bonds_according_to_index(atoms, index, circle)[0]
    O_links = get_H_bonds_according_to_index(atoms, index, circle)[1]
    water_H_index = get_H_bonds_according_to_index(atoms, index, circle)[2]
    # oxygen_index = [key for key,v in O_links.items()]

    H_bond_pair = []
    for i  in  O_links:
        deg1 = []
        deg2 = []
 
            


        O_accepter_index = i
        O_doner_index = index
        
        H1_index = water_H_index[0]
        H2_index = water_H_index[1]

        
        angle1 = atoms.get_angle(i, index, H1_index, mic=True)
    
        dis1 = atoms.get_distance(i, index, mic=True)
        H_bond_length1 = atoms.get_distance(index, H1_index)
        H_bond_length2 = atoms.get_distance(index, H2_index)
        if angle1 < 30 and  dis1 < 3.5 and H_bond_length1 < 3.0:
            p = [index,H1_index,i]
            # pair1.append(p)
            deg1.append(angle1)
            H_bond_pair.append(p)

        angle2 = atoms.get_angle(i, index, H2_index, mic=True)
    
        dis2 = atoms.get_distance(i, index, mic=True)

        if angle2 < 30 and  dis2 < 3.5 and H_bond_length2 < 3.0:
            p = [index,H2_index,i]
            # pair1.append(p)
            deg1.append(angle2)
            H_bond_pair.append(p)
            # tpair = pair1 + pair2
            # H_bond_pair.append(tpair)
            # if pair1 != []:
            #     H_bond_pair.append(pair1)
            # if pair2 != []:
            #     H_bond_pair.append(pair2)
            # print(f'pair1: {pair1}, pair2: {pair2}')
        # for l in H_bond_pair:
        #     if l == []:
        #         H_bond_pair.remove(l)
    H_bond_pair1 = []
    
    for i in H_links:
        for j  in O_links:
            
            dis = atoms.get_distance(j, index, mic=True)
            angle = atoms.get_angle(index, j, i,mic=True)
            dis_OH = atoms.get_distance(j,i, mic=True)
            H_bond_length1 = atoms.get_distance(index, i)
            if angle < 30 and  dis < 3.5 and dis_OH <= 1.3 and H_bond_length1 < 3:
                p = [index,i,j]
                # pair1.append(p)

                H_bond_pair1.append(p)
                
    H =   H_bond_pair + H_bond_pair1 
    counts = len(H)
    print(f'H_bond_pair: {H}, counts: {counts}')
    print(len(O_links))
    return counts

def get_h_bonds(atoms):
    
    """
    Return the number of H bonds near the surface
    
    """
    
    O_links = get_nearest_Oxygen(atoms, 3.5, 'Ru', 0, 4.9)
    # oxygen_index = [key for key,v in O_links.items()]
    h, o = get_water_according_to_distance(atoms, 4.5, 'Ru')
    water_group, b, c, d = assert_distance(atoms, o, h, 1.2)
    H_bond_pair = []
    for k, v in O_links.items():


        for i, j  in water_group.items():
            deg1 = []
            deg2 = []
            # pair1 = []
            # pair2 = []
            if i != k:

                O_accepter = atoms.get_positions()[k]
                O_doner = atoms.get_positions()[i]
                
                H1_index = water_group[i][0]
                H2_index = water_group[i][1]
                H1 = atoms.get_positions()[H1_index]
                H2 = atoms.get_positions()[H2_index]
                
                angle1 = get_angle(O_accepter-O_doner, H1-O_doner)
            
                dis1 = np.linalg.norm(O_accepter - H1)
                
                if angle1 < 35 and  dis1 < 3.0:
                    p = [k,i,H1_index]
                    # pair1.append(p)
                    deg1.append(angle1)
                    H_bond_pair.append(p)

                angle2 = get_angle(O_accepter-O_doner, H2-O_doner)
                dis2 = np.linalg.norm(O_accepter - H2)
                if angle2 < 35 and  dis2 < 3.0:
                    p = [k,i,H2_index]
                    # pair2.append(p)
                    deg2.append(angle2)
                    H_bond_pair.append(p)

                # tpair = pair1 + pair2
                # H_bond_pair.append(tpair)
                # if pair1 != []:
                #     H_bond_pair.append(pair1)
                # if pair2 != []:
                #     H_bond_pair.append(pair2)
                # print(f'pair1: {pair1}, pair2: {pair2}')
            # for l in H_bond_pair:
            #     if l == []:
            #         H_bond_pair.remove(l)

    counts = len(H_bond_pair)
    print(f'H_bond_pair: {H_bond_pair}, counts: {counts}')
    print(len(O_links))
    return counts
def get_surface_water(atoms, cdl):
    surface_Z = 8.3
    O_index = [atom.index for atom in atoms if atom.symbol == 'O' and 0.1 <= atom.position[2] - surface_Z <= cdl]
    return O_index
    
def get_h_bonds_test(atoms, cdl, max_thikness):
    
    """
    Return the number of H bonds near the surface
    
    """
    
    O_links = get_surface_water(atoms, cdl)
    # oxygen_index = [key for key,v in O_links.items()]
    near_surface_water = len(O_links)
    counts = []
    for k in O_links:
        count = get_h_bonds_signle_water(atoms, k, 3.5)
        counts.append(count)
        
  
    return sum(counts), O_links
def get_neareaset_H_index(atoms, donor_index):
    
    H_index = [atom.index for atom in atoms if atom.symbol == 'H']
    index = []
    distance = atoms.get_distances(donor_index, H_index, mic=True)
    for i in range(len(distance)):
        if distance[i] <= 4.0:
            index.append(H_index[i])
            
    
    return index





# from ase.io import read
# atoms = read('CONTCAR')
# a = get_h_bonds_test(atoms)
# b = get_neareaset_H_index(atoms, 64)
# view(atoms)
# def ArrayBottomK(top_k,arr):
#     bottom_k_idx=arr.argsort()[0:top_k]
# #     return bottom_k_idx  
# print(a)

def H_bonds_plots(file, start_frame, frame=None):

    traj = Trajectory(file, 'r' )
    total_frame = len(traj)   
    if frame is None:
        frame = total_frame
    
    h_bonds = []
    water_molecule = []
    for i in range(start_frame, frame):
        # total_frame = (1000)


        img = traj[i]
        

        # img.set_pbc = [True]
        h_number, O = get_h_bonds_test(img,3.5,4.9)

        water_molecule.append(len(O))
        h_bonds.append(h_number)
        print(O)
    time = np.arange((frame-start_frame))
    h_bond = np.array(h_bonds)
    water_molecule = np.array(water_molecule)
    
    # print(h_bonds)
    # print(np.mean(h_bonds))
    per_H_bonds = h_bond / water_molecule
    fig = plt.figure( figsize = ( 15, 10 ) )
    ax1 = fig.add_subplot( 2, 1, 1 )
    ax1.plot( time/1000*0.5, h_bond, label = 'H distribution' )
    ax1.set_xlabel( ' time (ps)' )
    ax1.set_ylabel( 'Avarage H_bonds (a.u.)' )
    ax2 = fig.add_subplot( 2, 1, 2 )
    ax2.plot( time, per_H_bonds, label = 'H distribution' )
    ax2.axhline(np.mean(per_H_bonds),color='red', linestyle='--')
    print(np.mean(per_H_bonds))
    ax2.set_xlabel( ' time (ps)' )
    ax2.set_ylabel( 'Avarage H_bonds_per_water (a.u.)' )
    
# h_bonds = H_bonds_plots('RuNi-1k.traj', 1000, frame=11000)
# traj = Trajectory('Ru.traj', 'r' )
# atoms = traj[0]
# a = get_h_bonds(atoms)
# print(a)


def counts_density(lists, section):
    
    """
    Classify lists according to section intervals
    """
    if section is None:
        raise ValueError("section cannot be empty")
    count = []
    if isinstance(lists, list):
        
       cuts = pd.cut(np.array(lists), section)
       counts = pd.value_counts(cuts,sort=False)
       for k, v in dict(counts).items():
           count.append(v)
    else:
         cuts = pd.cut(lists, section)
         counts = pd.value_counts(cuts, sort=False)
         for k, v in dict(counts).items():
             count.append(v)
    return count

def get_z_distance(atoms, atom_symbol, surface_number,max_z,dpoints):
    
    """
    atom_symbol: slab symbol
    surface_number: the number of surface atoms
    max_z: the distance between surface and max_z
    dpoints: the number of points between surface and max_z
    
    """ 
       
    top_metal_index = atoms[[ at.index for at in atoms if at.symbol == atom_symbol ]]
    maxz_atom_symbol = np.mean(np.sort(top_metal_index.positions[ :, 2 ])[::-1][:surface_number])
    H_positionZ = [i.position[2]-maxz_atom_symbol for i in atoms if i.symbol == 'H']
    H_positionZ = [i for i in H_positionZ if i < max_z]
#    oxygen_index = [i.index for i in atoms if i.symbol == 'O' and i.position[2] - maxz_atom_symbol <= thickness]
    
    O_positionZ = [i.position[2]-maxz_atom_symbol for i in atoms if i.symbol == 'O']
    O_positionZ = [i for i in O_positionZ if i < max_z]
    section = np.linspace(0,max_z,dpoints+1).tolist()
    H_counts = counts_density(H_positionZ, section)
    O_counts = counts_density(O_positionZ, section)
    
    return np.array(H_counts), np.array(O_counts)
# traj = Trajectory('Ru.traj', 'r' )
# atoms = traj[0]
# a,b = get_z_distance(atoms, 'Ru', 9, 7, 61)
# time = np.linspace(0,7,60)
# fig = plt.figure( figsize = ( 3.5, 3.5 ) )
# ax1 = fig.add_subplot( 2, 1, 1 )
# ax1.plot( time, b, label = 'H distribution' )
# ax1.set_xlabel( 'Radial distance (A)' )
# ax1.set_ylabel( 'Density (a.u.)' )

def plot_oxygen_distance_to_nearest_surface(file,atom_symbol, surface_number,max_z,dpoints, start_frame, frame=None):
    traj = Trajectory(file, 'r')
    if frame is None:
        
        
       frame = len(traj)
    oxy = np.zeros(dpoints)
    hydro = np.zeros(dpoints)
    for i in range(start_frame,frame):
        img = traj[i]
        H_counts, oxy_counts = get_z_distance(img, atom_symbol, surface_number, max_z, dpoints)
        oxy += oxy_counts
        hydro += H_counts
    distance = np.linspace(0,max_z,dpoints)
    oxygen_density  = oxy/(frame-start_frame)
    hydrogen_density = hydro/(frame-start_frame)
    fig = plt.figure( figsize = ( 14.5, 8.5 ) )
    ax1 = fig.add_subplot( 2, 1, 1 )
    ax1.plot( distance, oxygen_density, label = 'H distribution' )
    ax1.set_xlabel( f'Distance From oxygen to nearest {atom_symbol} (A)' )
    ax1.set_ylabel( 'Oxygen Density (a.u.)' )
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot( distance, hydrogen_density, label = 'H distribution' )
    ax2.set_xlabel( f'Distance From hydrogen to nearest {atom_symbol} (A)' )
    ax2.set_ylabel( 'Hydrogen Density (a.u.)' )
    
    
# print(b)
# a = plot_oxygen_distance_to_nearest_surface('RuNi-1k.traj', 'Ru', 7, 7, 61,1000)
def water_dipole_orientation(atoms, min_pos, max_pos, threshold):
    
    h,o = get_water_according_to_space(atoms, min_pos, max_pos, 'Ru')
    cos_phi = []
    if len(o) == 0:

        cos_phi.append(0)
    else:
        water_group, distance_group, angle_list, position_group = assert_distance(atoms, o, h, threshold)
    
    
        for k, v in water_group.items():
            O_pos = atoms.get_positions()[k]
            H1_index = water_group[k][0]
            H2_index = water_group[k][1]
            H1_pos = atoms.get_positions()[H1_index]
            H2_pos = atoms.get_positions()[H2_index]
            vec1 = (H1_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            vec2 = (H2_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            water_bisetion_angle = vec1 + vec2
            # print(vec1)
            normal = np.array([0,0,1])
            cos_angle = np.dot(normal, water_bisetion_angle) / (np.linalg.norm(normal) * np.linalg.norm(water_bisetion_angle))
            cos_phi.append(cos_angle)

    
    return np.mean(cos_phi)
def water_dipole_near_atoms(atoms, position_z, symbol, threshold):
    
    h,o = get_water_according_to_atom(atoms, position_z, symbol, 1)
    cos_phi = []
    if len(o) == 0:

        cos_phi.append(0)
    else:
        water_group, distance_group, angle_list, position_group = assert_distance(atoms, o, h, threshold)
    
    
        for k, v in water_group.items():
            O_pos = atoms.get_positions()[k]
            H1_index = water_group[k][0]
            H2_index = water_group[k][1]
            H1_pos = atoms.get_positions()[H1_index]
            H2_pos = atoms.get_positions()[H2_index]
            vec1 = (H1_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            vec2 = (H2_pos - O_pos) / np.linalg.norm((H1_pos - O_pos))
            water_bisetion_angle = vec1 + vec2
            # print(vec1)
            normal = np.array([0,0,1])
            cos_angle = np.dot(normal, water_bisetion_angle) / (np.linalg.norm(normal) * np.linalg.norm(water_bisetion_angle))
            cos_phi.append(cos_angle)

    
    return cos_phi

def plot_water_dipole_based_atoms(file,position_z, symbol, start, frame=None):
    
    traj = Trajectory(file, 'r')
    if frame is None:
        
        
       frame = len(traj)
    angle = []
    for i in range(start,frame):
        atoms = traj[i]
        an = water_dipole_near_atoms(atoms, position_z, symbol, 1.2)
        angle.extend(an)
    
    plt.hist(angle,bins=1000)
    plt.xlabel( 'Cos_phi' )
    plt.ylabel( 'Frequency (a.u.)' )


#a = plot_water_dipole_based_atoms('test1.traj', 7, 'Ru',2000)

def plot_water_dipole_orientation(file, max_z, dpoints, start, frame=None):
    traj = Trajectory(file, 'r')
    if frame is None:
        
        
       frame = len(traj)
    Z_space = np.linspace(0, max_z, dpoints).tolist()
    # print(Z_space)
    total_angle = np.zeros(dpoints-1)
    for i in range(start,frame):
        angle = []
        for j in range(len(Z_space)-1):
        
            img = traj[i]
            cos_phi = water_dipole_orientation(img, Z_space[j], Z_space[j+1], 1.2)
            angle.append(cos_phi)
        
        total_angle += np.array(angle)
    
    distance = np.linspace(0,max_z,dpoints-1)
    angle_phi  = total_angle/frame
    fig = plt.figure( figsize = ( 3.5, 3.5 ) )
    ax1 = fig.add_subplot( 2, 1, 1 )
    ax1.plot( distance, angle_phi, label = 'H distribution' )
    ax1.set_xlabel( 'Z coordinate (A)' )
    ax1.set_ylabel( 'Cos_phi (a.u.)' )
    
def plot_water_dipole_orientation(file, max_z, dpoints, start, frame=None):
    traj = Trajectory(file, 'r')
    if frame is None:
        
        
       frame = len(traj)
    Z_space = np.linspace(0, max_z, dpoints).tolist()
    # print(Z_space)
    total_angle = np.zeros(dpoints-1)
    # total_degree = np.zeros(dpoints-1,180)
    results = np.zeros((dpoints-1, 180))  
    total_dict ={} # data type {0.1 : {0: 0, 1: 1, 2: 1, 3: 3, 4: 2, 5: 1, 6: 1, 7: 0, 8: 0}}
    for i in range(start,frame):
        angle = []
        # degrees = []
        result = np.zeros((dpoints-1, 180))
        for j in range(len(Z_space)-1):
        
            img = traj[i]
            cos_phi, degree = water_dipole_orientation(img, Z_space[j], Z_space[j+1], 1.2)
            angle.append(cos_phi)
            # degrees.append(degree)
            deg = degree
            
            # [deg.extend(i) for i in degree]
            deg.sort()

            for m in deg:
                if 0 <= m - 1 <= 180:
                   result[j,m-1] +=1
                # b[deg[m]] = degree.count(deg[m])
            
            # total_dict[i] = b
      
            # for k, v in total_dict.items():
            #     for ik, iv in v.items():
            #         result[j,ik] = iv
            
        
        total_angle += np.array(angle)
        results += result
    #     deg = []
    #     [deg.extend(i) for i in degrees]
    #     deg.sort()
    #     b = {}
    #     for m in range(len(deg)):
    #         b[deg[m]] = degree.count(deg[m])
        
    #     total_dict[i] = b
    # results = np.zeros((dpoints-1, 181, 1))        
    # for k, v in total_dict.items():
    #     for ik, iv in v.items():
    #         results[i,ik,0] = iv
        
        
            
    x = np.linspace(0, max_z, dpoints-1) 
    y = np.linspace(0, 180, 181)        
    X,Y =np.meshgrid(x,y)   
    print(X.shape)
    print(results[8:12,:])
    print(np.max(results[8:12,:]))
    
    distance = np.linspace(0,max_z,dpoints-1)
    angle_phi  = total_angle/frame
    fig = plt.figure( figsize = ( 3.5, 3.5 ) )
    ax1 = fig.add_subplot( 2, 1, 1 )
    ax1.plot( distance, angle_phi, label = 'H distribution' )
    
    ax1.set_xlabel( 'Z coordinate (Å)' )
    ax1.set_ylabel( 'Cos_phi (a.u.)' )
    with open("{}_distance-cos_phi.txt".format(file), "w") as file1:
        # 将时间和数据同时写入文件，每行一个时间和对应的数据
        for time, value in zip(distance, angle_phi):
            file1.write(f"{time} {value}\n")   
            
    if os.path.exists('{}.npy'.format(file)):
        print('Already has output file')
    else:
        
        np.save('{}.npy'.format(file),results.T)
        
        
    plt.show()
    import seaborn as sns
    y = range(0,180,30)
    f, ax = plt.subplots(figsize=(20, 8))
    sns.heatmap(results.T, cmap='rocket_r', cbar_kws={'label': 'Frequency'})
    plt.xlabel('',fontsize=20, color='k') #x轴label的文本和字体大小
    plt.ylabel('F(φ)',fontsize=20, color='k') #y轴label的文本和字体大小
    plt.xticks([0,7],fontsize=20) #x轴刻度的字体大小（文本包含在pd_data中了）
    plt.yticks([0,7],fontsize=20) #y轴刻度的字体大小（文本包含在pd_data中了）
    # plt.title('title',fontsize=20) #图片标题文本和字体大小
#    plt.imshow(np.log10(results),extent=[0.1,7, 1,181], cmap='viridis', aspect='auto')  # 使用 viridis 调色板，也可以选择其他的
    # plt.colorbar()  # 添加颜色条
    
    plt.savefig('RuNi-1k-1000-11000.png', dpi=600)
    plt.show()
    
    

# a = plot_water_dipole_orientation('RuNi-1k.traj',7,81,1000, frame=11000)    

# a = plot_water_dipole_orientation('RuNi-1k.traj',7,41,1)
# print(a)

# view(atoms)
# # print(atoms.get)
# print(water_group)


# atoms.write('POSCAR1')
# a = water_dipole_orientation(atoms, 0, 2.5, 1.2)
# print(a)


def radical_distance_plot(file,atom_symbol, surface_number,max_z,dpoints, start_frame, end_frame=None):
    
    traj = Trajectory(file, 'r' )
    if end_frame is None:
       total_frame = len(traj)   
    
    else:
        total_frame = end_frame

    H_array = np.zeros(dpoints) 
    O_array = np.zeros(dpoints)
    
    for i in range(start_frame,total_frame):
        # total_frame = (1000)
      
        img = traj[i]
        img.set_pbc = [True]
        h, o = get_z_distance(img, atom_symbol, surface_number, max_z,dpoints)
        H_array += h
        O_array += o
        
    
    time = np.linspace(0,max_z,dpoints)
    fig = plt.figure( figsize = ( 3.5, 3.5 ) )
    ax1 = fig.add_subplot( 2, 1, 1 )
    ax1.plot( time, H_array, label = 'H distribution' )
    ax1.set_xlabel( 'Radial distance (A)' )
    ax1.set_ylabel( 'Density (a.u.)' )
    ax1 = fig.add_subplot( 2, 1, 2 )
    ax1.plot( time, O_array, label = 'O distribution' )
    ax1.set_xlabel( 'Radial distance (A)' )
    ax1.set_ylabel( 'Density (a.u.)' )
    
# a = radical_distance_plot('RuNi-1k.traj', 'Ru', 9, 7, 41, 1)
def count_near_circle(atoms,pair_A, pair_B, min_r, max_r):
    
    A_link = get_global_nearest_Oxygen(atoms, pair_A, pair_B, min_r, max_r)
    number = []
    counts = []
    for k, v in A_link.items():
        number.append(k)
        if v is None:
            count = 0
            counts.append(count)
        else:
          count = len(v)
          counts.append(count)
    return counts

def count_surface_near_circle(atoms,pair_A, pair_B, min_r, max_r):
    
    A_link = get_global_nearest_Oxygen(atoms, pair_A, pair_B, min_r, max_r)
    number = []
    counts = []
    for k, v in A_link.items():
        number.append(k)
        if v is None:
            count = 0
            counts.append(count)
        else:
          count = len(v)
          counts.append(count)
    return counts

def pair_correlation_function(atoms, pair_A, pair_B ,  distance, interval=None):
    
    if interval is None:
        interval = 41
    
    section = np.linspace(0,distance,interval)  
    dr = distance/(interval-1)  

    mean = np.zeros(len(section))

    section = np.linspace(0,distance,interval).tolist()  
    mean_coordinate = []
    for j in section:
        a = count_near_circle(atoms,pair_A,pair_B,j,j+dr)
        mean_coordinate.append(np.sum(np.array(a)/(4*np.pi*dr*(j+0.001)**2)))
    mean += mean_coordinate
    
    
    return mean
 
    
    

# img1 = traj[9000]
# O_link = get_global_nearest_Oxygen(img1, 3.5, 'Ru', 3.0,3.9)
# lists = [i for i in range(10)]
# number = []
# counts = []
# for k, v in O_link.items():
#     number.append(k)
#     if v is None:
#         count = 0
#         counts.append(count)
#     else:
#       count = len(v)
#       counts.append(count)
# array = np.zeros(len(number))
# print(number)
# print(counts)
def pcf_plot(file,pair_A, pair_B,  distance, interval, start_frame, dframe, frame=None):
    
    section = np.linspace(0,distance,interval) 
 
    traj = Trajectory(file, 'r' )
    if frame is None:
       frame = len(traj)
    else:
        frame = frame
       
    mean2 = np.zeros(len(section))

        
    for i in range(start_frame,frame, dframe):
        atoms = traj[i]
        mean = pair_correlation_function(atoms, pair_A, pair_B,  distance, interval)
        mean2 += mean
      
    time = section
    fig = plt.figure( figsize = ( 3.5, 3.5 ) )
    ax1 = fig.add_subplot( 2, 1, 1 )
    ax1.plot( time, mean2/((frame-start_frame)/dframe), label = 'H distribution' )
    ax1.set_xlabel( 'Radial distance (A)' )
    ax1.set_ylabel( 'Density (a.u.)' )    
            
        
        
a = pcf_plot('RuNi-1k.traj', 'O', 'O', 6, 61, 9,1,19)


# img1.write('POSCAR')
# b = img1.get_positions()[88]
# cell = img1.get_cell()

# f = img1.get_positions()[64]
# d = f+cell[0]+cell[1]
# dis1 = np.linalg.norm(d - b)
# print(a)
# print(b)
# print(d)
# print(cell)
# a = radical_distance_plot('Ru.traj', 'Ru', 16, 10, 61)






