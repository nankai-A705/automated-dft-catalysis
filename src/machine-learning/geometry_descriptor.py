# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 08:41:08 2022

@author: hongye qin
@email:290720931@qq.com
"""

from ase.build import fcc111, add_adsorbate, molecule
from ase.geometry import get_layers,get_distances
from ase.io import read, write
from ase.visualize import view
from collections import Counter
import numpy as np
from tools import get_polygon_area, distance_between_coordinate, get_volume_include_coordinate
Cu = fcc111('Cu',(4,4,5),periodic=True)
H2O = molecule('H2O')
add_adsorbate(Cu, H2O, 2.5, (7.696,2.82))
Cu.center(vacuum=10, axis=2)
Cu.symbols[68]='Ni'
ref = ['Fe','Co','Ni','Cu','Zn']

def adsorbates_center(atoms):
    
    """
    Return the center coordination of adsorbates
    
    """
    atom_index = [atom.index for atom in atoms if atom.tag == 0] # Return adsorbates atom index
    position = []
    for i in atom_index:
        position.append(atoms.get_positions()[i])
    
    positions_array = np.array(position)
    
    return np.mean(positions_array,axis=0)

def get_min_z_adsorbates(atoms):
    
    """
    Return minmum Z index of adsorbates atoms
    
    """
    atom_index = [atom.index for atom in atoms if atom.tag == 0] # Return adsorbates atom index
    position = []
    for i in atom_index:
        position.append(atoms.get_positions()[i][2])
    min_index = position.index(min(position[i] for i in range(len(position))))
    
    return atom_index[min_index]

def get_closest_atoms(atoms, ads_id,start_layer, end_layer,one_hot_code_for_layer,layer=None):
    
    """
       ads_id : adsorbates index 
       
       start_layer : from 0 to end_layer  
       
       end_layer : the layer index need to be considered
       
       one_hot_code_atom_for_layer: the number of atoms to be condiered for each of layer
       
       get_colsest_atoms(Cu, 20 ,1 , 4, 7)
       
       Should be mind that adsorbates need to place at center of slab
    """
    
    layer_index = [index for index in range(start_layer,end_layer)]
    if layer is not None:
       atom_index = [atom.index for atom in atoms if atom.tag != 0]
    #layer = get_layers(atoms, (0,0,1))
    else:
       atom_index = [atom.index for atom in atoms if atom.tag in layer_index]
    #atom_index = [j for j in range(len(atoms)) if layer[i] in range(end_layer)]
    #atoms_index.extend(atom_index)
    distances = atoms.get_distances(ads_id, atom_index)
   
    idx = [idx for _, idx in sorted(zip(distances, atom_index))][1:] # Return closest atom index except adsorbate   
    first_layer = []
    second_layer = []
    third_layer = []
    for i, id in enumerate(idx):
        if atoms[id].tag == layer_index[0]:
            first_layer.append(id)
        elif atoms[id].tag == layer_index[1]:
             second_layer.append(id)
        elif atoms[id].tag == layer_index[2]:
            third_layer.append(id)
    if len(first_layer) >= one_hot_code_for_layer:
       first_layer = first_layer[:one_hot_code_for_layer] 
    else:
        raise Exception('the atom number of first layer less than {} or far away from adsorbates, you should redefine one-hot code number[first, senceond, third]'.format(one_hot_code_for_layer))
    if len(second_layer) >= one_hot_code_for_layer:
       second_layer = second_layer[:one_hot_code_for_layer] 
    else:
        raise Exception('the atom number of second layer less than {} or far away from adsorbates, you should redefine one-hot code number[first, senceond, third].'.format(one_hot_code_for_layer))
    if len(third_layer) >= one_hot_code_for_layer:
       third_layer = third_layer[:one_hot_code_for_layer] 
    else:
        raise Exception('the atom number of third layer less than {} or far away from adsorbates, you should redefine one-hot code number[first, senceond, third]'.format(one_hot_code_for_layer))       
    
    one_hot_code_atom_index = [x for l in (first_layer, second_layer, third_layer) for x in l]
    
    tags = []
    for i in atom_index:
        if atoms[i].tag in layer_index:
           tags.append(atoms[i].tag) 
    #atoms_index_with_tag = zip(atom_index, tags) 
    return one_hot_code_atom_index, idx 


def get_closest_atoms_with_zone(atoms, ads_id,start_layer, end_layer,zone_number,layer=None):
    
    """
       ads_id : adsorbates index or adsorbates coordinate
       
       
       start_layer : from 0 to end_layer  
       
       end_layer : the layer index need to be considered
       
       zone_number: the number of atoms to be condiered for each of layer [6,3,3] first layer with 6 atoms
       
       get_colsest_atoms(Cu, 20 ,1 , 4, [6,3,3])
    """
    cell = atoms.get_cell()
    
    layer_index = [index for index in range(start_layer,end_layer)]
    if layer is not None:
       atom_index = [atom.index for atom in atoms if atom.tag != 0]
    #layer = get_layers(atoms, (0,0,1))
    else:
       atom_index = [atom.index for atom in atoms if atom.tag in layer_index]
    #atom_index = [j for j in range(len(atoms)) if layer[i] in range(end_layer)]
    #atoms_index.extend(atom_index)
    if isinstance(ads_id,int):
       distances = atoms.get_distances(ads_id, atom_index)
    else:
        positions = [atoms[i].position for i in atom_index]
        distances = get_distances(np.array(ads_id),positions,cell=cell,pbc=True)[1][0]
   
    idx = [idx for _, idx in sorted(zip(distances, atom_index))] # Return closest atom index except adsorbate   
    first_layer = []
    second_layer = []
    third_layer = []
    for i, id in enumerate(idx):
        if atoms[id].tag == layer_index[0]:
            first_layer.append(id)
        elif atoms[id].tag == layer_index[1]:
             second_layer.append(id)
        elif atoms[id].tag == layer_index[2]:
            third_layer.append(id)
    if len(first_layer) >= zone_number[0]:
       first_layer = first_layer[:zone_number[0]] 
    else:
        raise Exception('the atom number of first layer less than {} or far away from adsorbates, you should redefine one-hot code number[first, senceond, third]'.format(zone_number[0]))
    if len(second_layer) >= zone_number[1]:
       second_layer = second_layer[:zone_number[1]] 
    else:
        raise Exception('the atom number of second layer less than {} or far away from adsorbates, you should redefine one-hot code number[first, senceond, third].'.format(zone_number[1]))
    if len(third_layer) >= zone_number[2]:
       third_layer = third_layer[:zone_number[2]] 
    else:
        raise Exception('the atom number of third layer less than {} or far away from adsorbates, you should redefine one-hot code number[first, senceond, third]'.format(zone_number[2]))       
    
    one_hot_code_atom_index = [x for l in (first_layer, second_layer, third_layer) for x in l]
    
    tags = []
    for i in atom_index:
        if atoms[i].tag in layer_index:
           tags.append(atoms[i].tag) 
    #atoms_index_with_tag = zip(atom_index, tags) 
    return one_hot_code_atom_index, idx 

def parse_formula(atoms, layer):
    """
    atoms type : Atoms ase
    return : chemical species ['Fe','Co']
             species number   ['2','1']
    
    """
    slab_symbols = [atom.symbol for atom in atoms if atom.tag != 0]
    slab_atom_index = [atom.index for atom in atoms if atom.tag != 0]
    adsorbate_symbol = [atom.symbol for atom in atoms if atom.tag == 0]
    identical_symbol = [x for x in slab_symbols if x in adsorbate_symbol]
    if len(identical_symbol) == 0:
        print('The symbol at slab and adsorbate are not identical')
    else:
        print('You should be carfully for counting the number of slab atoms')
    symbol_dict = {}
    for sym in slab_symbols:
        symbol_dict[sym] = slab_symbols.count(sym)
    chemical_species = []
    species_number = []   
    for key in symbol_dict:
        chemical_species.append(key)
        species_number.append(symbol_dict[key])
    
    return chemical_species, species_number

def encode_with_neareast_atom_number(atoms,ads_id, start_layer,end_layer,zone_size):
    
    """
    zone_size = [first layer atoms,second layer atoms,third layer atoms]
    ref = ['Fe','Co','Ni','Cu','Zn']
    if the number of Fe in zone contain 3 atoms 
    return [3,0,0,0,0]
    """
    total_index = get_closest_atoms_with_zone(atoms, ads_id, start_layer, end_layer, zone_size)[0]
    first_layer_index = total_index[0:zone_size[0]]
    second_layer_index = total_index[int(zone_size[0]):zone_size[1]+zone_size[0]]
    third_layer_index = total_index[zone_size[1]+zone_size[0]:zone_size[2]+zone_size[1]+zone_size[0]]
    
    first_layer_sym = [atoms[index].symbol for index in first_layer_index]
    second_layer_sym = [atoms[index].symbol for index in second_layer_index]
    third_layer_sym = [atoms[index].symbol for index in third_layer_index]
    first_layer_count = Counter(first_layer_sym)
    second_layer_count = Counter(second_layer_sym)
    third_layer_count = Counter(third_layer_sym)
    
    natoms = sum(1 for atom in atoms if atom.tag != 0)
    first_layer_code = []
    for i in ref:
        if dict(first_layer_count).__contains__(i):
           first_layer_code.append(first_layer_count[i])  
        else:
           first_layer_code.append(0)
    second_layer_code = []
    for i in ref:
        if dict(second_layer_count).__contains__(i): # Determine if exists 'Ni', if exists return number of atoms, else return 0
           second_layer_code.append(second_layer_count[i])  
        else:
           second_layer_code.append(0)
    third_layer_code = []
    for i in ref:
        if dict(third_layer_count).__contains__(i):
           third_layer_code.append(third_layer_count[i])  
        else:
           third_layer_code.append(0)
    
    
    one_hot_code_atom_number = [x for l in (first_layer_code, second_layer_code, third_layer_code) for x in l] 

    return one_hot_code_atom_number


def get_surface_zone_area_near_ads(atoms, tolrence=0.8, adsorbate=None):
    
    """
               ---------
              |         \
             |     ads    \
              \          |
               \        |
                --------              
    return area and volume         
    """       
    
    surface_atom_index = [atom.index for atom in atoms if atom.tag == 1]
    
    surface_atom_position = [atom.position for atom in atoms if atom.tag == 1]
    ads_position = adsorbates_center(atoms)
    distance = distance_between_coordinate(ads_position, surface_atom_position)
    assert len(surface_atom_index) == len(distance)
    idx = [idx for _, idx in sorted(zip(distance, surface_atom_index))]# Return surface atom index 
    dist = [_ for _, idx in sorted(zip(distance, surface_atom_index))]
    
    # To make sure the area with same order we fix all in hexagon zone
    ads_position_xyz = adsorbates_center(atoms)
    
    if dist[1] - dist[0] > tolrence and dist[7] - dist[6] > tolrence: #  top
       idx = idx[1:7]
       area_atom_position_xy = [atoms[i].position[:2] for i in idx ]
       area_atom_position_xyz = [atoms[i].position for i in idx ]
       return get_polygon_area(area_atom_position_xy), get_volume_include_coordinate(ads_position_xyz, area_atom_position_xyz)
   
    elif dist[3] - dist[2] > tolrence and dist[2] -  dist[1] < 0.2 and dist[1] - dist[0] < 0.2:    # hollow
        idx = idx[:3]
        area_atom_position_xy = [atoms[i].position[:2] for i in idx ] 
        area_atom_position_xyz = [atoms[i].position for i in idx ] 
        return get_polygon_area(area_atom_position_xy) * 6 , get_volume_include_coordinate(ads_position_xyz, area_atom_position_xyz) * 6
    
    elif dist[4] - dist[2] > tolrence and dist[1] - dist[0] < 0.3 and dist[3]  - dist[2] < 0.3: # bridge
        idx = idx[:4]
        area_atom_position_xy = [atoms[i].position[:2] for i in idx ]
        area_atom_position_xyz = [atoms[i].position for i in idx ]
        return get_polygon_area(area_atom_position_xy) * 3 , get_volume_include_coordinate(ads_position_xyz, area_atom_position_xyz) * 3
    else:
        idx = idx = idx[:4]
        area_atom_position_xy = [atoms[i].position[:2] for i in idx ] 
        area_atom_position_xyz = [atoms[i].position for i in idx ]
        return get_polygon_area(area_atom_position_xy) * 3 , get_volume_include_coordinate(ads_position_xyz, area_atom_position_xyz) * 3

    
    
    
    
def encode_surface_area_close_adsorbates(atoms, ref):
   
    """
    ref : reference energy without doping or pure metal surface area
    strain caused by doping will be include in this module
    
    
    """
    
    return get_surface_zone_area_near_ads(atoms)[0] - ref

def encode_atoms_number_near_ads(atoms, cutoff=5):
    """
    Return the atoms number near the ads according to the cutoff ricial
    cutoff 
    """
    
    atoms_position_except_ads = [atom.position for atom in atoms if atom.tag != 0]
    atoms_index_except_ads = [atom.index for atom in atoms if atom.tag != 0]
    ads_position = adsorbates_center(atoms)
    distance = distance_between_coordinate(ads_position, atoms_position_except_ads)
    idx = [idx for _, idx in sorted(zip(distance, atoms_index_except_ads))]# Return surface atom index 
    dist = [_ for _, idx in sorted(zip(distance, atoms_index_except_ads))]
    atoms_number_near_ads = len([d for d in dist if d < cutoff])
    return atoms_number_near_ads

def encode_volume_include_ads(atoms):
    
    
    return get_surface_zone_area_near_ads(atoms)[1]

#view(Cu)
#j = encode_with_neareast_atom_number(Cu, 20, 1, 4, [6,3,3])
#i  = (Cu, 20,1, 4, [6,3,3])    
# print(encode_volume_include_ads(Cu))
# print(adsorbates_center(Cu))
# #print(Cu.tag)
# print(Cu.get_positions()[20])
# a = get_closest_atoms_with_zone(Cu, [7.696, 2.82, 20.837] ,1 , 4, [6,3,0])
# print(a[1])
# from ase.visualize import view
# view(Cu)








        