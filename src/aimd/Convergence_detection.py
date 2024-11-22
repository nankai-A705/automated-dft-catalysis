# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
author: Hongye Qin
email: 290720931@qq.com
"""

import os
import numpy as np
import matplotlib.pyplot as plt

class aimd_outcar:
    
    def __init__(self, file=None):
        if file is None:
            self.outcar = 'OUTCAR'
        else:
            self.outcar = file
 
            
    
    
    def get_free_energy(self):
        
        outcar = [line.strip() for line in open(self.outcar)]
        
        text = os.popen("grep 'free  energy' OUTCAR|awk ' {print $5}'").read()
        text.replace('\n'," ")
        data = text.split()
        energy = [float(value) for value in data]
       

        potim = []
        for line in open(self.outcar):
            # if 'free energy' in line:
            #     line = line.replace('\n'," ")
            #     #print(line)
            #     lines = line.split()
            #     print(float(lines[4]))
                
                
                # energy.append(float(lines[4]))
            if 'POTIM  =' in line:
                lines = line.split()
                potim.append(lines[2])
                
            
                
     
        self.energy = energy
        self.potim = float(potim[0])
        
        return self.energy, self.potim
        
    
    
    def plot_time_energy(self, start, end=None):
        
        self.energy, self.potim = self.get_free_energy()
        
        if end is None:
            
           end = (len(self.energy) + 1) * self.potim 
           
        else:
            end = end
        
        time = np.arange(self.potim*start, end, self.potim)
        energy = self.energy[start-1:]
        # print(time)
        # print(end)
        # print(len(energy))
        max_y = max(energy) + 5
        min_y = min(energy) - 5
        plt.plot(time, energy)
        plt.xlabel( ' Time (ps)' )
        plt.ylabel( 'Energy (eV)' )
        plt.ylim(min_y, max_y)
        
        
        
    

a = aimd_outcar()
# b=a.get_free_energy()
c = a.plot_time_energy(10)
# print(len(b))
# print(len(a.energy)*0.5)
# print(np.arange(0, 1.5,0.5))           