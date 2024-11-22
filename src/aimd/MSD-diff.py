#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 19:26:28 2024

@author: hongye Qin
"""

import sys
import re
import numpy as np
import matplotlib.pyplot as plt

#setting grid for histogram

potim = 0.5                               #timestep from INCAR file
readfile = open("XDATCAR-test","r")        #input XDATCAR file in format XDATCAR.TEMP
#temp=re.sub("XDATCAR.",,sys.argv[1])  #extracts temperature from input file name
z=0                                     #counter
natoms=0                                #number of atoms in XDATCAR file
posion = []                             #atom positions in Cartesian coordinates
confcount = 0                           #number of structures in XDATCAR file
direct=[]                               #number of time steps for each structure in XDATCAR file
a=[]                                    #lattice parameter in 1st dimension
b=[]                                    #lattice parameter in 2nd dimension
c=[]                                    #lattice parameter in 3rd dimension
#read in XDATCAR file
line=readfile.readline()

while (line):
  z=z+1
  line.strip()
  line=re.sub('^',' ',line)
  y=line.split()
  if (z==2):
     scale=float(y[0])
  if (z==3):
     a.append(float(y[0]))
     a.append(float(y[1]))
     a.append(float(y[2]))
     a_len=(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])**0.5
  if (z==4):
     b.append(float(y[0]))
     b.append(float(y[1]))
     b.append(float(y[2]))
     b_len=(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])**0.5
  if (z==5):
     c.append(float(y[0]))
     c.append(float(y[1]))
     c.append(float(y[2]))
     c_len=(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])**0.5
  if (z==7):

      natoms = np.sum(np.array([int(i) for i in y]))

  if (y[0]=="Direct"):
     direct.append(int(y[2]))
     posion.append([])
     for i in range(0,natoms):
        line=readfile.readline()
        line.strip("\n")
        line=re.sub('^',' ',line)
        f=line.split()
        cartpos_x=a[0]*float(f[0])+a[1]*float(f[1])+a[2]*float(f[2])
        cartpos_y=b[0]*float(f[0])+b[1]*float(f[1])+b[2]*float(f[2])
        cartpos_z=c[0]*float(f[0])+c[1]*float(f[1])+c[2]*float(f[2])
        #positions of ions for each structure are obtained here
        posion[confcount].append([cartpos_x,cartpos_y,cartpos_z])
     confcount=confcount+1
  line=readfile.readline()
readfile.close

#calculate diffusion coefficient
#skip first 10 configurations corresponding to 300 fs
atom_index = [i for i in range(63,142)]
print(natoms)

MSD =  []
for i in range(10,confcount):
    d=0.0
    for j in atom_index:
      x_diff=posion[i][j][0]-posion[0][j][0]
      #if length is larger than 0.5 (in crystallographic coordinates) then we have to shift atom
      #due to periodic image to obtain the shortest distance.
      if (abs(x_diff)>(0.5*a_len)):
         if (x_diff<0):
            x_diff=x_diff+a_len
         elif (x_diff>0):
            x_diff=x_diff-a_len
      y_diff=posion[i][j][1]-posion[0][j][1]
      if (abs(y_diff)>(0.5*b_len)):
         if (y_diff<0):
            y_diff=y_diff+b_len
         elif (y_diff>0):
            y_diff=y_diff-b_len
      z_diff=posion[i][j][2]-posion[0][j][2]
      if (abs(z_diff)>(0.5*c_len)):
         if (z_diff<0):
            z_diff=z_diff+c_len
         elif (x_diff>0):
            z_diff=z_diff-c_len
      d=d+x_diff**2.0+y_diff**2.0+z_diff**2.0  
   # print(d)
    MSD.append(d/len(atom_index))



#print diffusion coefficient (in Ang^2/ps) vs temperature (in K)
dtotal = d
d=dtotal/(confcount-1-10)/len(atom_index)/6.0
time=(direct[confcount-1]-direct[10])*potim/10**3.0 #conversion to ps
print(d/time)
print(len(MSD))
import numpy as np
fig = plt.figure( figsize = ( 15, 10 ) )
ax1 = fig.add_subplot( 2, 1, 1 )
ax1.plot( np.arange(10,confcount), MSD, label = 'H distribution' )
print(MSD[0])
ax1.set_xlabel( ' time (ps)' )
ax1.set_ylabel( 'Avarage H_bonds (a.u.)' )
