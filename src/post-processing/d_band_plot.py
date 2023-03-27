



from ase.io import read

from xml_parser import get_pdos
import numpy as np
"""
                                                    
"""
def smooth(y, box_pts):
    """Smooth the noisy density of state distribution using convolution operater.
    Args:
       y (array): one-dimensional input array
       box_pts (int): average box size parameter
    Returns:
       array: one-dimentional array of smoothed density of state values
    """
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def moment(x, y, n):
    """moments function to calculate the porbability distribution characteristics of density of states.
    Args:
       x (array): energy values
       y (array): density of states
       n (int): order parameter of moments function
    Returns:
       float: moment descriptor 
    """
    p = x**n * y
    return np.trapz(p, x)/np.trapz(y, x)

def density_moments(energies, dos, mode=None):
    """Calculate the moment descriptors for the density of states distributions.
    Args:
       energies (array): energy values with respect to fermi energy
       dos (array): density of states
    returns:
       list: moment characteristics including filling, center, sigma, skewness, kurtosis
    """
    # smooth the noisy density state 
    dos_rev = smooth(dos[:], 15)# smoothing function
    # determine the index of first non-positive value  
    Ind = np.argmax(energies>0)
    # calculate the moment descriptors for the density of states 
    if mode is None:
       dos_rev = smooth(dos[:], 15)# smoothing function
       filling = np.trapz(dos_rev[0:Ind], energies[0:Ind])/np.trapz(dos_rev, energies)
       center = moment(energies[:], dos_rev[:], 1)             
       sigma_c = np.sqrt(moment(energies[:]-center, dos_rev[:], 2))
       skewness = moment(energies[:]-center, dos_rev[:], 3)/sigma_c**3 
       kurtosis = moment(energies[:]-center, dos_rev[:], 4)/sigma_c**4 
       return [filling, center, sigma_c, skewness, kurtosis]
    if mode == 'center':
       return moment(energies, dos_rev, 1)
    if mode == 'c_center':
        return moment(energies, dos, 1)

atoms = read('CONTCAR')
atom_index = [atom.index for atom in atoms]
atom_symbol = [atom.symbol for atom in atoms]
metal_index = [atom.index for atom in atoms if atom.symbol == 'Ni']
d_band = []

for j in atom_index:
    if atoms[j].symbol != 'N':       
        energy = get_pdos('./','vasprun.xml',j+1,['dxy','dyz','dxz','x2-y2','dz2'])[0]
        dos = get_pdos('./','vasprun.xml',j+1,['dxy','dyz','dxz','x2-y2','dz2'])[1]
        d_band.append(density_moments(energy, dos, 'center'))
    else:
        # energy = get_pdos('./','vasprun.xml',j,['s'])[0]
        # dos = [0] * len(energy)      
        # d_band.append(density_moments(energy, dos, 'center'))
        d_band.append(0)

f = open('dbands.txt','w')
f.write(str(d_band))
f.close()

