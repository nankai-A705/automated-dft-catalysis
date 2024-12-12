# -*- coding: utf-8 -*-
"""
Created on Thu May 16 11:34:59 2024

@author: Win
"""
import numpy as np
def cubic(h1,k1,l1,h2,k2,l2):
    
    cos_phi =  (h1*h2 + k1*k2 + l1*l2) / np.sqrt((h1**2 + k1**2 + l1**2)*(h2**2 + k2**2 + l2**2))
    return np.degrees(np.arccos(cos_phi))


# A = cubic(1, 0, 0, 0, 1, 0)
# print(A)

def orthorhombic(lattic_a, lattic_c, f1, f2):
    
    Numerator = (f1[0]*f2[0]+f1[1]*f2[1]) + 0.5*(f1[0]*f2[2] + f2[0]*f1[1]) + (f1[2]*f2[2]*0.75*lattic_a**2)/lattic_c**2
    Denominator1 = f1[0]**2+f1[1]**2 + f1[0]*f1[1] + 0.75 *f1[2]**2 * lattic_a**2 / lattic_c**2 
    Denominator2 = f2[0]**2+f2[1]**2 + f2[0]*f2[1] + 0.75 *f2[2]**2 * lattic_a**2 / lattic_c**2 
    cos_phi = Numerator/ np.sqrt(Denominator1*Denominator2)
    
    return np.degrees(np.arccos(cos_phi))



def Hexagonal(lattic_a, lattic_c, f1, f2):
    
    Numerator = (f1[0]*f2[0]+f1[1]*f2[1])+0.5*(f1[0]*f2[1] + f1[1]*f2[0]) + 0.75 * (lattic_a**2/lattic_c**2)*f1[2]*f2[2]
    
    print(Numerator)
    
    Denominator1 = (f1[0]**2+f1[1]**2) + f1[0]*f1[1]  + 0.75 * (lattic_a**2/lattic_c**2)*f1[2]**2
    Denominator2 = (f2[0]**2+f2[1]**2) + f2[0]*f2[1] + 0.75 * (lattic_a**2/lattic_c**2)*f2[2]**2 
    cos_phi = Numerator/ np.sqrt(Denominator1*Denominator2)
    print(cos_phi)
    
    return np.degrees(np.arccos(cos_phi))



b = Hexagonal(4.603, 4.304, [0,0,2], [1,1,1])
print(b)