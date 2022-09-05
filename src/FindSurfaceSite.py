# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 21:02:35 2021

@author: hongye qin
@email： hongyechin@qq.com
"""

import itertools
import os
import numpy as np
from ase.build import fcc111
from ase.io import read
from pymatgen import vis
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs, Slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import in_coord_list_pbc
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from ase.atoms import Atoms
from ase.atom import Atom
import random
from ase.build import add_adsorbate
from ase.constraints import FixAtoms





def add_adsorbates_with_fixed_layer(slab, ads,heights,position,layers):
    add_adsorbate(slab, ads, height=heights, position=position)
    mask = [atom.tag < layers for atom in slab]
    c = FixAtoms(mask=mask)
    slab.set_constraint(c)
    
class FindSurfaceSite:
    
    """
    This Class
    """
    
    def __init__(self,slab):
        self.slab = slab
        self.mols = {
                    'O' : Atoms('O'),
                    'OH' : Atoms('OH', positions = [[0, 0, 0], [0.6, 0.7, 0.15]]),
                    'OOH' : Atoms('OOH', positions = [[0, 0, 0], [0.8, 0.9, 0.0], [1.4, 1.6, 0.45]]),
                    }
                    
            
         
    
    def get_c_vec(self,slab):
        """
         Convenience function which returns the unit vector aligned
         with the miller index
         Notes: This function only works for ASE Module.
         For example: a_vector: [1,0,0], b_vector: [0,1,0]
         return np.cross(a_vector,b_vector) c_vector:[0,0,1]
    
       """
        vec = np.cross(slab.get_cell()[0],slab.get_cell()[1])
        return vec / np.linalg.norm(vec)

    def get_surface_site_ase(self,slab, tolerance,**excludes):
        
        """
          This function will return a surface site coorodinate list which contain excludes['symbol'],
          tolerance: surface site Z value mines subsurface site Z value
          for examples: if someone wants to get the SnO2(110) surface coordinate without O coordinate
          Input file: get_surface_site(slab,0.1,'O')
          
        """
        if isinstance(self.slab, Atoms):
           self.slab = self.slab
        elif isinstance(self.slab, Structure):
            self.slab = AseAtomsAdaptor().get_atoms(self.slab)
        # Sometimes The slab with vacuum created by pymatgen may want to find surface site with this module
            
        Z_coordinate = np.dot(slab.get_positions(),self.get_c_vec(slab)) # Get all atoms Z coorodinate
        mask = (Z_coordinate - np.amax(Z_coordinate)) >= -tolerance# Return surface atoms index
        mask_list = np.where(mask)[0]
        surface_sites = []
        surface_symbols = [] # Total surface sites
        surface_symbols_with_strip = []
        surface_sites_with_strip = [] # Surface sites without excludes
        
        for i in mask_list:
            surface_sites.append(slab.get_positions()[i])
            surface_symbols.append(slab.symbols[i])
            
            if slab[i].symbol in excludes: continue
            surface_symbols_with_strip.append(slab.symbols[i])
            surface_sites_with_strip.append(slab.get_positions()[i])
        if excludes == None:
            return surface_sites,surface_symbols
        else:
            return surface_sites,surface_symbols, surface_sites_with_strip, surface_symbols_with_strip
        
    def get_surface_site_pymatgen(self):
        
        if isinstance(self.slab, Slab):
            self.slab = self.slab
        else:
            self.slab = AseAtomsAdaptor().get_structure(self.slab)
        
        asf_slabs = AdsorbateSiteFinder(self.slab)
        surface_sites = asf_slabs.find_adsorption_sites(distance=2.0) #Type : Dict
        return surface_sites 
    
    def get_surface_site_ontop(self):
        
        surface_sites = self.get_surface_site_pymatgen()
        ontop_sites = []
        for pos in surface_sites['ontop']:
            ontop_sites.append(pos)
            
        return ontop_sites #Type : List
    def get_surface_site_hollow(self):
        
        surface_sites = self.get_surface_site_pymatgen()
        hollow_sites = []
        for pos in surface_sites['hollow']:
            hollow_sites.append(pos)
            
        return hollow_sites #Type : List
    
    def get_surface_site_bridge(self):
        
        surface_sites = self.get_surface_site_pymatgen()
        bridge_sites = []
        for pos in surface_sites['bridge']:
            bridge_sites.append(pos)
            
        return bridge_sites #Type : List
       
    def get_surface_site_ase_delaunay(self, tolerance):
        
        """
        This function returns top, hollow, bridge sites of slab(type:Atoms) which 
        generated with ASE module
        
        """
        
        surface_site = self.get_surface_site_ase(self.slab, tolerance)[0]
        surface_sites = np.array(surface_site) 
        from scipy.spatial import Delaunay

        top = surface_sites
        tri = Delaunay(surface_sites[:, :2])     # points set
        pos_nodes = surface_sites[tri.simplices] # triangular nodes coordination

        hollow = np.array([t.sum(axis=0)/3 for t in pos_nodes])

        bridge = []
        for i in pos_nodes:
            bridge.append((i[0]+i[1])/2)
            bridge.append((i[0]+i[2])/2)
            bridge.append((i[1]+i[2])/2)
        bridge = np.array(bridge)   
        return top, hollow, bridge
        
        
        
       
# if __name__ == '__main__':        
        
#    Cu = Structure.from_file("Cu.cif")
#    slabgen = SlabGenerator(Cu, (1,1,1), 10, 10)
#    slabs = slabgen.get_slabs()
#    #print(len(slabs))
#    Cu_110 = fcc111('Cu',(3,3,6),6)
#    Cu_110.center(vacuum=10, axis=2)
#    Cu_111 = slabs[0]
#    if isinstance(Cu_110, Slab):
#        Cu_110 = Cu_110
#    else: 
#        Cu_110 = AseAtomsAdaptor().get_structure(Cu_110)
#    print(Cu_110)
#    surface = FindSurfaceSite(Cu_111).get_surface_site_ontop
   
#    #sites = surface.get_surface_site_pymatgen()
#    print(surface)
   
   
#    """
#    Test Atoms type
#    """
#    Cu_ase = fcc111('Cu',(3,3,6),6)
#    Cu_ase.center(vacuum=10, axis=2)
#    ase_site = FindSurfaceSite(Cu_ase).get_surface_site_ase_delaunay(0.02)
#    top = ase_site[0]
#    hollow = ase_site[1]
#    bridge = ase_site[2]
   
   
#    fig = plt.figure()
#    ax = fig.add_subplot(1, 1, 1)
#    ax.triplot(top[:,0], top[:,1], triangles = Delaunay(top[:, :2]).simplices, color='grey',)
#    ax.plot(hollow[:,0],  hollow[:,1],  'ok', label='Hollow')
#    ax.plot(bridge[:,0],  bridge[:,1],  'or', label='Bridge')
#    ax.plot(top[:,0], top[:,1], 'ob', label='Atop')
#    plt.legend(loc='lower left')
#    plt.axis('scaled')
        
        
        
        
        
        
        
        
        
        
        
        
    