#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 22:37:55 2024

@author: hongye qin
"""

import os
import matplotlib.pyplot as plt
import numpy as np

if os.path.exists('grad.dat'):
    os.system('rm grad.dat')
os.system("grep cc REPORT |awk '{print $3}' > xxx")
os.system("grep b_m REPORT |awk '{print $2}' > fff")
os.system("paste xxx fff >> grad.dat")
os.system("sort -n grad.dat > grad_.dat")
os.system("mv grad_.dat grad.dat")
os.system("rm xxx")
os.system("rm fff")


def get_data(file):
    """
    get x , y coordinate, total points, wave numbers, intensity
    
    """
    x = []
    y = []

    f = open(file)
    for line in f:
        line = line.strip('\n')
        line = line.split('\t')
        num = len(line)
        if num == 2:
           x.append(float(line[0]))
           y.append(float(line[1])) 
    f.close()
    
    ttg = []
    tg = 0.0
    for i in range(1, len(x)):
        gg = 0.5*(x[i]-x[i-1])*(y[i]+y[i-1])
        tg += gg
        ttg.append(tg)
    transition_state_index = ttg.index(max(ttg))
    transition_state_energy = max(ttg)
    print("transition_state_index:", transition_state_index)
    print("transition_state energy:", transition_state_energy)
    fig, ax = plt.subplots()
    x = x[1:]
    ax.plot(x, ttg, linewidth=2.0)
    ax.set_xlabel('Collective variable (Ang)')
    ax.set_ylabel('Free energy (eV)')
    with open("output.txt", "w") as file:
    # 将时间和数据同时写入文件，每行一个时间和对应的数据
         for time, value in zip(x, ttg):
             file.write(f"{time} {value}\n")   
    plt.show()

a = get_data('grad.dat')

def get_bluemoon():
    # read from a single-step REPORT file with IBLUEOUT=T
    # returns: cv, lambda, |z|^(-1/2), GkT, |z|^(-1/2)*(lambda+GkT), T_inst
    data_raw = [line.strip() for line in open('REPORT')]
    cv = []
    lamda = []
    mass = []
    gkt = []
    bm = []
    T = []
    
    for line in data_raw:
        if 'cc>' in line:
            line = line.strip('\n')
            line = line.split('\n')
            s = line[0]
            part = s.split()
            
            cv.append(float(part[2]))
            
        
        if 'b_m>' in line:
           line = line.strip('\n')
           line = line.split('\t')
           s = line[0]
           part = s.split()
           a = float(part[1])
           b = float(part[2])
           c = float(part[3])
           d = float(part[4])
           lamda.append(a)
           mass.append(b)
           gkt.append(c)
           bm.append(d)
        
        if 'tmprt>' in line:
            line = line.strip('\n')
            line = line.split('\t')
            s = line[0]
            part = s.split()
            T.append(float(part[2]))
           


    
    return cv, lamda, mass, gkt, bm, T
        
cv, bm1, bm2, bm3, bm4, temp = get_bluemoon()      

fe_grad = np.array(bm4)/np.array(bm2)



fig, ax = plt.subplots()

ax.plot(np.array(cv), fe_grad, linewidth=2.0)

ax.set_xlabel('Collective variable (Ang)')
ax.set_ylabel('Force (eV)')
# print(a)
# atoms = Trajectory('test.traj')
# atom1  = atoms[1].write('POSCAR0')
# atom2 = atoms[3177].write('POSCAR1')
# atom3 = atoms[4046].write('POSCAR2')
 