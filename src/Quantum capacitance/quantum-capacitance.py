# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 20:00:53 2023

@author: win
"""
import numpy as np
import scipy.constants as cn
import matplotlib.pyplot as plt
from ase.visualize import view
import pandas as pd
from DosAnalysis import Dos
from ase.io import read

# 1 F = C*2 / J  https://en.wikipedia.org/wiki/Farad

Boltzmann = 8.617333262e-5
uF_per_cm2_constant = 6.241418e40 # Constant to convert default units to mF/cm^-2
F_per_gram_constant = 3.75867e42 # Constant to convert default units to F/g
F_per_gram_constant = 3.75867e42 # Constant to convert default units to F/g




# dos = Dos(xml=xml)
class Quantum_Capacitance(Dos):
    
    def __init__(self, xml):
        
        self.E_min = -7.0
        self.E_max = 7.0
        self.parse = Dos(xml)
        self.dos = self.get_range_dos()
        self.mass = self.parse.get_system_mass()
        self.E = self.get_range_energy()
        self.T = 300
        self.area = self.parse.get_cross_section_area()
        self.mass = self.parse.get_system_mass()                          
        self.phi_min=-1.0
        self.phi_max = 1.0

        

    
    def get_dos(self):
        
        
        tdos, up_dos, down_dos  = self.parse.get_tdos(spin=2)
        up_dos = pd.DataFrame(up_dos).loc[:,'total']
        down_dos = pd.DataFrame(down_dos).loc[:,'total']
        
        
        return np.abs(up_dos) + np.abs(down_dos)
    
    def get_energy(self):
        
        dos = self.parse.get_tdos(spin=2)
        
        return  pd.DataFrame(dos[1]).loc[:,'energy']
    

          
        
    
    def sech(self,x):
        
        return 2 / (np.exp(x) + np.exp(-x))



    def get_Cq(self, phi):
        
        E = np.array(self.E)
        dos = np.array(self.dos)
        prefix = cn.e ** 2 / (4 * Boltzmann * self.T)
        sech_new = []
        for e in E:
            sech2 = self.sech((e - phi)/(2*Boltzmann * self.T)) ** 2
            sech_new.append(sech2)
        sech1 = np.array(sech_new)
        inner = dos * sech1
        integrate = np.trapz(inner, E)
        mF_per_cm2 = integrate * prefix * uF_per_cm2_constant / self.area
        F_per_gram = integrate * prefix * F_per_gram_constant / self.mass
        
        
        return mF_per_cm2, F_per_gram
    
    def get_Cq_phi(self):
        
        phi = np.linspace(self.phi_min, self.phi_max, 200).tolist()
        Cq = []
        for i in phi:
            cq = self.get_Cq(i)[0]
            Cq.append(cq)
        
        return phi, Cq
    
    def get_gram_Cq_phi(self):
        
        phi = np.linspace(self.phi_min, self.phi_max, 200).tolist()
        Cq = []
        for i in phi:
            cq = self.get_Cq(i)[1]
            Cq.append(cq)
        
        return phi, Cq
    
    def get_range_energy(self):
        
        tdos, up_dos, down_dos = self.parse.get_tdos(spin=2)
        
        energy = up_dos[(up_dos['energy'] >= self.E_min) & (up_dos['energy'] <= self.E_max)].loc[:,'energy']
        

        
        return energy
    

        
       
        
    def get_range_dos(self):
        
        tdos, up_dos, down_dos = self.parse.get_tdos(spin=2)
        # up_dos = pd.DataFrame(up_dos).loc[:,'total']
        # down_dos = pd.DataFrame(down_dos).loc[:,'total']
        
        dos_up = up_dos[(up_dos['energy'] >= self.E_min) & (up_dos['energy'] <= self.E_max)].loc[:,'total']
        dos_down = down_dos[(down_dos['energy'] >= self.E_min) & (down_dos['energy'] <= self.E_max)].loc[:,'total']
        
        return abs(dos_up) + abs(dos_down)
    
    def Fermi_Dirc_distribution(self, phi):
        
        E = np.array(self.E)
        

        
        f = 1 / ((np.exp(E-phi)/ Boltzmann * self.T) + 1)
        
        return f
    
    def get_delt_Q(self):
        
        E = np.array(self.E)
        
        phi = np.linspace(self.phi_min, self.phi_max, 200).tolist()
        
        dos = np.array(self.dos)
        fermi_dirc = []
        for i in phi:
            fermi_dirc_difference = self.Fermi_Dirc_distribution(i) - self.Fermi_Dirc_distribution(0)
            fermi_dirc.append(fermi_dirc_difference)
        
        inner = np.array(fermi_dirc) * dos
        
        integrate = np.trapz(inner, E)
        Q  = 1e24 * integrate * cn.e / self.area # 1e24  C transfrom uC 10e6 A2 to cm2 10e16
        
        
        return phi, Q
    

    def plot_phi_Cq(self):
        
        
        phi, Cq = self.get_Cq_phi()
        
        font = {'family' : 'arial',
                'weight' : 'bold',
                'size'   : 16}

        plt.rc('font', **font)
        
        plt.scatter(phi, Cq)

        plt.xlabel('Potential (V)')
        plt.ylabel('Quantum Capacitance $\mathregular{(µF/cm2)}$') 
        plt.show()
        # plt.save('phi_Cq.png',dpi=400)
        
    def plot_gram_phi_Cq(self):
        
        
        phi, Cq = self.get_gram_Cq_phi()
        
        font = {'family' : 'arial',
                'weight' : 'bold',
                'size'   : 16}

        plt.rc('font', **font)
        
        plt.scatter(phi, Cq)

        plt.xlabel('Potential (V)')
        plt.ylabel('Quantum Capacitance $\mathregular{(µF/g)}$') 
        plt.show()
        # plt.save('phi_Cq.png',dpi=400)
    def plot_E_dos(self):
        
        
        
        font = {'family' : 'arial',
                'weight' : 'bold',
                'size'   : 16}

        plt.rc('font', **font)
        
       
        plt.plot(self.E,self.dos)
        plt.xlabel('Energy (V)')
        plt.ylabel('Density of state') 
        plt.show()
    def plot_phi_Q(self):
        
        
        
        font = {'family' : 'arial',
                'weight' : 'bold',
                'size'   : 16}

        plt.rc('font', **font)
        
        phi, Q = self.get_delt_Q()
        plt.plot(phi,Q)
        plt.xlabel('Energy (V)')
        plt.ylabel('uC/cm2') 
        plt.show()    
    
        
        




Q = Quantum_Capacitance('vasprun.xml')
# a = Q.get_Cq(0.2)
b = Q.plot_gram_phi_Cq()
c = Q.plot_phi_Cq()
# d = Q.plot_phi_Q()
# a = Q.mass
# print(a)
    
# E = np.linspace(-5, 5, 101)
# np.random.seed(0)
# dos = np.random.randint(5, size=(101))
# # plt.scatter(E, dos)
# plt.show()
# Ef = np.linspace(0, 0.8,100).tolist()
# c= []
# print(Ef)
# for i in Ef:
#     C = Cq_2D(E,dos,i,300,20)
#     c.append(C)
# # c = Cq(E,dos,0.6,300,20)
# plt.scatter(Ef,c)
# # print(c)
# plt.show()
# def get_cross_section_area(v1,v2):
    
#     return  abs(v1[0]*v2[1] - v1[1]*v2[0])

# a = get_cross_section_area(cell[0], cell[1]) 

# a = sech(E)
# print(a)
# y1 = np.array([0,1,2,3])
# y2 = np.array([1,2,4,5])
# print(y1*y2)