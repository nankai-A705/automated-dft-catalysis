# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 21:15:43 2021

@author: hongye qin
@Time:2021/08/24
"""

import itertools
import os
import numpy as np
from ase.atoms import Atoms
from ase.atom import Atom
from ase.build import fcc111
from ase.constraints import FixAtoms
from pymatgen import vis
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import VoronoiNN
import time
from ase.thermochemistry import HarmonicThermo
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from copy import deepcopy
from ase.geometry import get_layers
from pymatgen.io.ase import AseAtomsAdaptor
from ase.build import add_adsorbate
from FindSurfaceSite import FindSurfaceSite
from ase.db import connect
from ase.calculators.vasp import Vasp

class Generate_Adsorbates_Structure():
    
    def __init__(self,slab,adsorbates,site,tolerance=None,coverage=None,mode=None):
        self.slab = slab
        self.adsorbates = adsorbates
        self.coverage = coverage
        self.tolerance = tolerance
        self.mode = mode
        self.site = site

    def generate_structure_with_adsorbates_ase(self):
        """
        This functions generate slab with full coverage of adsrobates 
        """
        
        if isinstance(self.adsorbates, Atoms):
          ads = self.adsorbates
        elif isinstance(self.adsorbates, Atom):
          ads = Atoms([self.adsorbates])
        
        from ase.db import connect
        db = connect('surface.db')
        adsorbates = self.adsorbates
        
        if self.coverage == None:
           if self.mode == 'all':
            
               surface_sites = FindSurfaceSite(self.slab).get_surface_site_ase(self.slab,self.tolerance)[0]
               slabs = []
               for pos in surface_sites:
                   add_adsorbate(self.slab, adsorbates,2.0, pos[:2])
                   slabs.append(self.slab)
               return slabs
           if self.mode == 'top':
               top_site = FindSurfaceSite(self.slab).get_surface_site_ase_delaunay(self.tolerance)[0]
               slab = self.slab
               slabs = []
               for pos in top_site:
                   add_adsorbate(slab, adsorbates,2.0, pos[:2])
                   slabs.append(slab)
                   
               return slabs  
           if self.mode == 'hollow':
               hollow_site = FindSurfaceSite(self.slab).get_surface_site_ase_delaunay(self.tolerance)[1]
               slabs = []
               for pos in hollow_site:
                   add_adsorbate(self.slab, adsorbates,2.0, pos[:2])
                   slabs.append(self.slab)
               return slabs
                   
           if self.mode == 'bridge':
               bridge_site = FindSurfaceSite(self.slab).get_surface_site_ase_delaunay(self.tolerance)[2]
               slabs = []
               for pos in bridge_site:
                   add_adsorbate(self.slab, adsorbates,2.0, pos[:2])
                   slabs.append(self.slab)
               return slabs
            
           
        else:
            
        
            coverage = self.coverage
            surface_sites = FindSurfaceSite(self.slab).get_surface_site_ase(self.slab,self.tolerance)[0]
            number_of_sites = len(surface_sites)
            effective_sites = round(coverage* len(surface_sites))
            site_index = np.ceil(np.linspace(0,number_of_sites-1,round(number_of_sites/effective_sites)))
            add_position = []
            for index in site_index:
                add_position.append(surface_sites[int(index)-1])
            structure_with_adsorbates = []
            
            for pos in add_position:
                add_adsorbate(self.slab, ads,2,pos[:2],4)
                
    
            return  self.slab
        
    
    def generate_structure_with_adsorbates_pymatgen(self, top=False, hollow=False, bridge=False):
                                                    
        """
        Return All slab sites with adsorbates ["TOP","HOLLOW,"BRIDGE"]
        """
        if isinstance(self.slab, Atoms):
            slab = AseAtomsAdaptor().get_structure(self.slab)
        if isinstance(self.adsorbates, Atoms):
            adsorbate = AseAtomsAdaptor().get_molecule(self.adsorbates)
        
        ads_slab = AdsorbateSiteFinder(slab).generate_adsorption_structures(
             adsorbate)
        
        if top:
            top_sites = FindSurfaceSite(self.slab).get_surface_site_ontop()
            slabs= []
            for site in top_sites:
                top_site_with_adsorbate = AdsorbateSiteFinder(slab).add_adsorbate(adsorbate, site)
                slabs.append(top_site_with_adsorbate)
            return slabs 
            
            
        return top_site_with_adsorbate
        if hollow:
            hollow_sites = FindSurfaceSite(self.slab).get_surface_site_hollow()
            for site in hollow_sites:
                hollow_site_with_adsorbate = AdsorbateSiteFinder(slab).add_adsorbate(adsorbate, site)
                return hollow_site_with_adsorbate
                
            
        
        
    def generate_structure_with_single_adsorbates(self):
        
        add_adsorbate(self.slab, self.adsorbates, 2.0, self.site[:2])
        return self.slab
        
        


###############################################################################
# Below this line are tools function for simulation
###############################################################################
from pandas import Series       
        
def threshold_cluster(Data_set,threshold):
    
    """
    Automate classification list to different list according to threshold
    data = [1,2,7,9,2.5,7.8,7.9]
    threshold : 2
    return [[1,2,2.5],[7,7.8,7.9,9]
            index = [0,1,3],[2,4,5,6]
    """
    #Return List to Dataframe
    stand_array=np.asarray(Data_set).ravel('C')
    stand_Data=Series(stand_array)
    index_list,class_k=[],[]
    while stand_Data.any():
        if len(stand_Data)==1:
            index_list.append(stand_Data.index.to_list())
            class_k.append(stand_Data.to_list())
            stand_Data=stand_Data.drop(stand_Data.index)
        else:
            class_data_index=stand_Data.index[0]
            class_data=stand_Data[class_data_index]
            stand_Data=stand_Data.drop(class_data_index)
            if (abs(stand_Data-class_data)<=threshold).any():
                args_data=stand_Data[abs(stand_Data-class_data)<=threshold]
                stand_Data=stand_Data.drop(args_data.index)
                new_list = [class_data_index]+args_data.index.to_list()
                index_list.append(new_list)
                class_k.append([class_data]+args_data.to_list())
            else:
                index_list.append([class_data_index])
                class_k.append([class_data])
                
    return index_list,class_k    
                
def get_new_tags(raw_data,threshold):
    
    """
    This function could return tags list according to atom position z value
    examples: raw_data : [1.2,2.6,7.8,1.2]
    usage: get_new_tags(raw_data,1)
    return [0,1,2.0]
    """
    index, data =  threshold_cluster(raw_data,threshold)
    tag = [i[0] for i in data] 
    sort = sorted(tag)
    tags = []
    for i in tag:
        tags.append(sort.index(i))

    new_tags = []
    for i in range(len(index)):
        l = [tags[i]] * int(len(index[i]))
        
        new_tags.extend(l)
    
    atom_index = []
    for i in range(len(index)):
        atom_index.extend(index[i])
    zip_tag = zip(atom_index, new_tags)
    zip_tag_sort = sorted(zip_tag)
    idx, tag_index = zip(*zip_tag_sort)
    return tag_index


def fix_layers(atoms, miller = (0, 0, 1), tol = 1.0, n = [0, 2]):
    '''
    Return Fixed Atoms according layer index
    '''
    layers = get_layers(atoms, miller, tol)[0]
    index = [j for j in range(len(atoms)) if layers[j] in range(n[0], n[1])]
    constraint = FixAtoms(indices=index)
    atoms.set_constraint(constraint)
    return  atoms   

def fix_layer_except_adsorbate(atoms):
    """
    Fixed all of atoms except for adsorbates
    """
    
    atom_index = [atom.index for atom in atoms ]
    positions_z = [atoms.positions[i][2] for i in atom_index]
    atoms_tag = get_new_tags(positions_z, 0.5) #tolerance need to be set with right value especially for hollow and bridge sites
    atoms.set_tags(atoms_tag)

    tags = atoms.get_tags().tolist()
    min_tag = min(tags, key=tags.count)
    fix_tags = [i for i in tags if i != min_tag]
    constraint = FixAtoms(indices=[atom.index for atom in atoms if atom.tag in fix_tags])
    atoms.set_constraint(constraint)  
    return atoms

def get_vib_zpe(atoms,temperature):
    """
    Get zpe, entropy, internal_energy, helmholtz of adsorbates on surface
    """
    
#   atoms = fix_layer_except_adsorbate(atoms)
    
    vib_calculator = Vasp(ibrion=5,
                          nfree=2,
                          potim=0.01,
                          nsw=1,
                          ismear=1,
                          sigma=0.2,
                          ediff=0.0000001,
                          encut=500,
                          lwave=False,
                          algo='Fast',
                          lasph=True,
                          lreal='Auto',
                          kpts=(1, 1, 1),
                          nelm=60,
                          ivdw=11,
                          lcharg=False,
                          ldipol=True,
                          idipol=3,
                          xc='pbe',
                          isym=-1,
                          lmaxmix=4)
    atoms.calc = vib_calculator
    atoms.get_potential_energy()
    vibrations = vib_calculator.get_vibrations() #
    harmonic_mode_energies=vibrations.get_energies()
    Thermo = HarmonicThermo(harmonic_mode_energies)
    entropy = Thermo.get_entropy(temperature)
    zpe = Thermo.get_ZPE_correction()
    internal_energy = Thermo.get_internal_energy(temperature)
    helmholtz = Thermo.get_helmholtz_energy(temperature)
    time.sleep(2)
    return zpe, entropy, internal_energy, helmholtz

def get_unique_id(atoms):
    pass

    
def get_atom(name, idx):
    
    """
    Get Atom object according to id 
    
    
    """
    db = connect(name)
    atoms = db.get_atoms(id=idx)
    return atoms
    
def get_id(file):
    
    bulk = connect(file)
    id_list = []
    for row in bulk.select():
        id_list.append(row.id)
    return id_list, len(id_list)

file = 'RuNi-relaxed.db'
db = connect(file)
adsorbates_db = connect('top-adsorbates-relaxed.db')
vasp = Vasp(ibrion=2,
            nsw=500,
            isif=2,
            ispin=2,
            isym=-1,
            ediff=1e-5,
            ediffg=-0.02,
            algo='Fast',
            lasph=True,
            ismear=1,
            lreal='Auto',
            sigma=0.2,
            kpts=(5, 5, 1),
            encut=600,
            nelm=60,
            ivdw=11,
            lorbit=11,
            lwave=False,
            lcharg=False,
            ldipol=True,
            idipol=3,
            xc='pbe',
           # ncore=16,
            lmaxmix=4)

number_of_atoms = get_id(file)[1] 
count = 0
for i in range(3,number_of_atoms):
    atoms = get_atom(file, i+1)
    top_site = FindSurfaceSite(atoms).get_surface_site_ase_delaunay(1)[0]
    hollow_site = FindSurfaceSite(atoms).get_surface_site_ase_delaunay(1)[1] 
    bridge_site = FindSurfaceSite(atoms).get_surface_site_ase_delaunay(1)[2]
    site_dict = {'top_site' : top_site}
#   site_dict = { 'hollow_site' : hollow_site, 'bridge_site' : bridge_site}
    for site_name, sites in site_dict.items():

        
        for j, pos in enumerate(sites):
            atom = get_atom(file, i+1)
            
            ads_slab = Generate_Adsorbates_Structure(atom,'H',pos).generate_structure_with_single_adsorbates()
            ads_slab = fix_layers(ads_slab)
            count += 1
            
            id_list = get_id(file)[0]
            if os.path.exists('vasprun.xml'):
                os.remove('vasprun.xml')
                os.remove('WAVECAR')
                os.remove('CHGCAR')
                os.remove('CHG')
                ads_slab.calc = vasp #relax slab with adsorbates
                energy = ads_slab.get_potential_energy()
                adsorbates_db.write(ads_slab,label_id=id_list[i],site_label=site_name)
                if os.path.exists('vasprun.xml'):
                    os.remove('vasprun.xml')
                    os.remove('WAVECAR')
                    os.remove('CHGCAR')
                    os.remove('CHG')
                    fix_slab = fix_layer_except_adsorbate(ads_slab)
                    zpe = get_vib_zpe(fix_slab,298.15)[0]
                    entropy = get_vib_zpe(fix_slab,298.15)[1]
                    internal_energy= get_vib_zpe(fix_slab,298.15)[2]
                    helmholtz = get_vib_zpe(fix_slab,298.15)[3]
                    adsorbates_db.update(count,zpe=zpe,entropy=entropy,internal_energy=internal_energy, helmholtz=helmholtz)
                else:
                    fix_slab = fix_layer_except_adsorbate(ads_slab)
                    zpe = get_vib_zpe(fix_slab,298.15)[0]
                    entropy = get_vib_zpe(fix_slab,298.15)[1]
                    internal_energy= get_vib_zpe(fix_slab,298.15)[2]
                    helmholtz = get_vib_zpe(fix_slab,298.15)[3]
                    adsorbates_db.update(count,zpe=zpe,entropy=entropy,internal_energy=internal_energy, helmholtz=helmholtz)
                    
            
                
            else:
                ads_slab.calc = vasp #relax slab with adsorbates
                energy = ads_slab.get_potential_energy() 
                adsorbates_db.write(ads_slab,label_id=id_list[i],site_label=site_name)
                if os.path.exists('vasprun.xml'):
                    os.remove('vasprun.xml')
                    os.remove('WAVECAR')
                    os.remove('CHGCAR')
                    os.remove('CHG')
                    fix_slab = fix_layer_except_adsorbate(ads_slab)
                    zpe = get_vib_zpe(fix_slab,298.15)[0]
                    entropy = get_vib_zpe(fix_slab,298.15)[1]
                    internal_energy= get_vib_zpe(fix_slab,298.15)[2]
                    helmholtz = get_vib_zpe(fix_slab,298.15)[3]
                    adsorbates_db.update(count,zpe=zpe,entropy=entropy,internal_energy=internal_energy, helmholtz=helmholtz)
                else:
                    fix_slab = fix_layer_except_adsorbate(ads_slab)
                    zpe = get_vib_zpe(fix_slab,298.15)[0]
                    entropy = get_vib_zpe(fix_slab,298.15)[1]
                    internal_energy= get_vib_zpe(fix_slab,298.15)[2]
                    helmholtz = get_vib_zpe(fix_slab,298.15)[3]
                    adsorbates_db.update(count,zpe=zpe,entropy=entropy,internal_energy=internal_energy, helmholtz=helmholtz)
                
        
       
        
        
        
        
        
        
######################################################################
#Below this line is test file
######################################################################
# number_of_atoms = get_id(file)[1] 
# def get_update_id(db):
#     number_of_atoms = get_id(db)[1] 
#     idx = 0
#     for i in range(number_of_atoms):
#         atoms = get_atom(db, i+1)
#         sites = FindSurfaceSite(atoms).get_surface_site_ase_delaunay(1)[0]
#         for j, pos in enumerate(sites):
            
#             idx +=1
    
    
#             return idx

# a = get_unique_id(file)
# print(a)
# for j, pos in enumerate(top_site):
#     atoms = atoms = get_atom(file, 1)
#     fix_atom = fix_layers(atoms)
#     ads_slab = Generate_Adsorbates_Structure(fix_atom,'H',pos).generate_structure_with_single_adsorbates()
#     fix_atom = fix_layers(ads_slab)
#     new_atom = fix_layer_except_adsorbate(fix_atom)
#     print(new_atom)
# zpe = get_vib_zpe(atoms,298.15)[0]
# entropy = get_vib_zpe(atoms,298.15)[1]
# s = get_vib_zpe(atoms,298.15)[2]
# i = get_vib_zpe(atoms,298.15)[3]
# print(zpe)
# print(s)
# print(entropy)
# print(i)

# slabs = []
# for pos in top_site:
#     Cu = fcc111('Cu',(3,3,6),6)
#     surface_sites,surface_symbols, surface_sites_with_strip, surface_symbols_with_strip = FindSurfaceSite(Cu).get_surface_site_ase(Cu,tolerance=0.1)       
# #        def fix_atom_according_layer
#     adsorbate = Atoms('OH', positions = [[0, 0, 0], [0.6, 0.7, 0.15]])       
#     Cu.center(vacuum=10.0, axis=2)
#     Cu = fix_layers(Cu)
#     ads_slab = Generate_Adsorbates_Structure(Cu, adsorbate,pos[:2]).generate_structure_with_single_adsorbates()
#     print(ads_slab)
# print(slabs)
