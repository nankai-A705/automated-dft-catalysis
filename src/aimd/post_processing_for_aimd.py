# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 20:02:46 2023

@author: win
"""

import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
import numpy as np
from ase.geometry import get_distances
from ase.visualize import view
from ase.geometry import wrap_positions

traj = Trajectory( 'test.traj', 'r' )
total_frame = len(traj)
tag_group = []
first_frame = traj[0]

#view(first_frame)
# print(total_frame)


def get_water_according_to_distance(atoms, thickness, atom_symbol):
    
    top_metal_index = atoms [[ at.index for at in atoms if at.symbol == atom_symbol ]]
    maxz_atom_symbol = np.max( top_metal_index.positions[ :, 2 ] )
    hydrogen_index = [i.index for i in atoms if i.symbol == 'H']
    oxygen_index = [i.index for i in atoms if i.symbol == 'O' and i.position[2] - maxz_atom_symbol <= thickness]
    
    if len(oxygen_index) == 0:
        
        raise ValueError('You should choose right frame for analysis or change the thickness of edl')
        
        
    else:
        
        water_group = {}
        for oxygen in oxygen_index:
            single_water_group = []
            single_hydrogen_group = []
            distances = atoms.get_distances(oxygen, hydrogen_index)
            bottom_k_index = ArrayBottomK(2, distances)
            for i in bottom_k_index:
                single_hydrogen_group.append(hydrogen_index[i])
                
       
        #single_hydrogen_group.append()
            water_group[oxygen] = single_hydrogen_group #{O index: H1 index, H2 index}
        return water_group 



def ArrayTopK(top_k,arr):
    top_k_idx=arr.argsort()[::-1][0:top_k]
    return top_k_idx

def ArrayBottomK(top_k,arr):
    bottom_k_idx=arr.argsort()[0:top_k]
    return bottom_k_idx

total_angle = []
total_distance = []
for i in range(total_frame):
    img = traj[i]
    img.set_pbc = [True]
    water_group = get_water_according_to_distance(img, 3.5, 'Ru')
    
    angle = []
    distance = []
    
    for k, v in water_group.items():
        
        a = img.get_angle(v[0], k, v[1])
        angle.append(a)
        dis = img.get_distance(k, v[0])
        O_cor = np.array(img.get_positions()[k])
        H_cor = np.array(img.get_positions()[v[0]])
        if dis > 2:
            O = np.array(O_cor)
            H = np.array(H_cor)
            index = np.argmax(abs(H-O))
            H_new = v[0] + img.get_cell()[index]
            dis = np.linalg.norm(O_cor - H_cor)
        
        distance.append(dis)
        
#     mean_angle = np.mean(np.array(angle))
#     mean_distance = np.mean(distance)
#     total_angle.append(mean_angle)
#     total_distance.append(mean_distance)
    

# angle = np.array(total_angle)
# distance = np.array(total_distance)
# time = np.arange( total_frame )
# fig = plt.figure( figsize = ( 3.5, 3.5 ) )
# ax1 = fig.add_subplot( 2, 1, 1 )
# ax1.plot( time, angle, label = 'water angle' )
# ax1.set_xlabel( 'Time/ps' )
# ax1.set_ylabel( 'angle [K]' )

# ax2 = fig.add_subplot( 2, 1, 2 )
# ax2.plot( time, distance, label = 'OH distance' )
# ax2.set_xlabel( 'Time/ps' )
# ax2.set_ylabel( 'distance [A]' )
# a = get_water_according_to_distance(first_frame, 10, 'Ru')  
# print(a)  
H = np.array([7.5,0.114,16.39])
O = np.array([12.8,9.3,15.4])
print(np.argmax(abs(H-O)))
# print(first_frame.get_positions(wrap=True)[140])      
cell = first_frame.get_cell()
a = wrap_positions([[14.265,6.093,13.896]], cell,pbc=[1,1,1])
print(a)        
        
    


# print(len(b))