# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 21:33:22 2023

@author: Hongye Qin
@email: 290720931@qq.com
"""





import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq, fftshift
from scipy import signal
import math

# Boltzmann Constant in [eV/K]
kb = 8.617332478E-5
# electron volt in [Joule]
ev = 1.60217733E-19
# Avogadro's Constant
Navogadro = 6.0221412927E23

c = 2.9979245899e10
def zero_padding(sample_data):
    """ A series of Zeros will be padded to the end of the dipole moment
    array (before FFT performed), in order to obtain a array with the
    length which is the "next power of two" of numbers.
    
    #### Next power of two is calculated as: 2**math.ceil(math.log(x,2))
    #### or Nfft = 2**int(math.log(len(data_array)*2-1, 2))
    """
    return int(2 ** math.ceil(math.log(len(sample_data), 2)))
################################################################################      
class xdatcar:
    """ Python Class for VASP XDATCAR """

    def __init__(self, File=None):
        if File is None:
            self.xdatcar = 'XDATCAR'
        else:
            self.xdatcar = File

        # time step of MD
        self.atom_index = None
        self.potim = None
        # mass per type
        self.mtype = None
        self.readoutcar()

        self.TypeName = None
        self.ChemSymb = None
        self.Ntype = None
        self.Nions = None
        self.Nelem = None
        self.Niter = None

        # position in Direct Coordinate
        self.position = None
        # position in Cartesian Coordinate
        self.positionC = None
        # Velocity in Angstrom per Femtosecond
        self.velocity = None
        self.std = 1000
        self.readxdat()

        self.mass_and_name_per_ion()
        # Temperature
        self.Temp = np.zeros(self.Niter-1)
        # Kinetic Energy
        self.Ken = np.zeros(self.Niter-1)
        # Time in femtosecond
        self.Time = np.arange(self.Niter-1) * self.potim
        self.getTemp()

        # Velocity Autocorrelation Function 
        self.VAF = None
        self.VAF2= None
        # Pair Correlation Function
        # self.PCF = None

    def mass_and_name_per_ion(self):
        # mass per ion
        self.mions = []
        self.ChemSymb = []

        if self.TypeName is None:
            self.TypeName = [chr(i) for i in range(65,91)][:self.Ntype]

        for i in range(self.Ntype):
            self.mions += [np.tile(self.mtype[i], self.Nelem[i])]
            self.ChemSymb += [np.tile(self.TypeName[i], self.Nelem[i])]

        self.mions = np.concatenate(self.mions)
        self.ChemSymb = np.concatenate(self.ChemSymb)

    def readxdat(self):
        """ Read VASP XDATCAR """

        inp = [line for line in open(self.xdatcar) if line.strip()]
        scale = float(inp[1])
        
        self.cell = np.array([line.split() for line in inp[2:5]],
                              dtype=float)
        self.cell *= scale

        ta = inp[5].split()
        tb = inp[6].split()
        # print ta, tb
        if ta[0].isalpha():
            self.TypeName = ta
            self.Ntype = len(ta)
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
        else:

            self.Nelem = np.array(tb, type=int)
            self.Nions = self.Nelem.sum()
            self.Ntype = len(tb)
            self.TypeName = None

        pos = np.array([line.split() for line in inp[7:]
                        if not line.split()[0].isalpha()],
                        dtype=float)
        self.position = pos.ravel().reshape((-1,self.Nions,3))
        self.Niter = self.position.shape[0]

        dpos = np.diff(self.position, axis=0)
        self.positionC = np.zeros_like(self.position)
        # apply periodic boundary condition
        dpos[dpos > 0.5] -= 1.0
        dpos[dpos <-0.5] += 1.0
        # Velocity in Angstrom per femtosecond
        for i in range(self.Niter-1):
            self.positionC[i,:,:] = np.dot(self.position[i,:,:], self.cell)
            dpos[i,:,:] = np.dot(dpos[i,:,:], self.cell) / self.potim

        self.positionC[-1,:,:] = np.dot(self.position[-1,:,:], self.cell)
        self.velocity = dpos
        
    def read_partial_xdat(self, atom_index):
        atom_index = [i-1 for i in atom_index]
        
        inp = [line for line in open(self.xdatcar) if line.strip()]
        scale = float(inp[1])
        
        self.cell = np.array([line.split() for line in inp[2:5]],
                              dtype=float)
        self.cell *= scale

        ta = inp[5].split()
        tb = inp[6].split()
        # print ta, tb
        if ta[0].isalpha():
            self.TypeName = ta
            self.Ntype = len(ta)
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
        else:

            self.Nelem = np.array(tb, type=int)
            self.Nions = self.Nelem.sum()
            self.Ntype = len(tb)
            self.TypeName = None
        
        
            

        pos = np.array([line.split() for line in inp[7:]
                        if not line.split()[0].isalpha()],
                        dtype=float)
        self.position = pos.ravel().reshape((-1,self.Nions,3))
        self.Niter = self.position.shape[0]
        
        self.partial_pos = self.position[:,atom_index,:]
        


        dpos = np.diff(self.partial_pos, axis=0)
        self.positionC = np.zeros_like(self.partial_pos)
        # apply periodic boundary condition
        dpos[dpos > 0.5] -= 1.0
        dpos[dpos <-0.5] += 1.0
        # Velocity in Angstrom per femtosecond
        for i in range(self.Niter-1):
            self.positionC[i,:,:] = np.dot(self.partial_pos[i,:,:], self.cell)
            dpos[i,:,:] = np.dot(dpos[i,:,:], self.cell) / self.potim

        self.positionC[-1,:,:] = np.dot(self.partial_pos[-1,:,:], self.cell)
        self.partial_velocity = dpos
        self.partial_atoms_number = len(atom_index)
        print(self.positionC.shape)
    

            


       

    def readoutcar(self):
        """ read POTIM and POMASS from OUTCAR """

        if os.path.isfile("OUTCAR"):
            # print "OUTCAR found!"
            # print "Reading POTIM & POMASS from OUTCAR..."

            outcar = [line.strip() for line in open('OUTCAR')]
            lp = 0; lm = 0; 
            for ll, line in enumerate(outcar):
                if 'POTIM' in line:
                    lp = ll
                if 'Mass of Ions in am' in line:
                    lm = ll + 1
                if lp and lm:
                    break

            # print outcar[lp].split(), lp, lm
            self.potim = float(outcar[lp].split()[2])
            self.mtype = np.array(outcar[lm].split()[2:], dtype=float)

    def getTemp(self, Nfree=None):
        """ Temp vs Time """

        for i in range(self.Niter-1):
            ke = np.sum(np.sum(self.velocity[i,:,:]**2, axis=1) * self.mions / 2.)
            self.Ken[i] = ke * 1E7 / Navogadro / ev
            if Nfree is None:
                Nfree = 3 * (self.Nions - 1)
            self.Temp[i] = 2 * self.Ken[i] / (kb * Nfree)
    
    def get_partial_diffusion_coefficient(self):
        

        
        self.partial_pos[self.partial_pos > 0.5] -= 1.0
        self.partial_pos[self.partial_pos <-0.5] += 1.0
        self.partial_diffusion = np.zeros_like(self.partial_pos)
        for i in range(self.Niter):
               self.partial_diffusion[i,:,:] = np.dot((self.partial_pos[i,:,:] - self.partial_pos[0,:,:]), self.cell)
        print(self.partial_diffusion)
        total_diffusion = np.sum(np.square(self.partial_diffusion),axis=1)
        print('partial diffusion shape:', self.partial_diffusion.shape)
        print('square diffusion shape:', total_diffusion.shape)
        
        diffusion_coefficient = np.sum(np.square(self.partial_diffusion))/self.Niter/ self.partial_atoms_number/6
        print('Niter:', self.Niter)
        time = self.Niter * self.potim / 10**3
        
        return diffusion_coefficient/time , total_diffusion

    def getVAF(self):
        """ Velocity Autocorrelation Function """

        # VAF definitions
        # VAF(t) = Natoms^-1 * \sum_i <V_i(0) V_i(t)>
        ############################################################
        # Fast Fourier Transform Method to calculate VAF
        ############################################################
        # The cross-correlation theorem for the two-sided correlation:
        # corr(a,b) = ifft(fft(a)*fft(b).conj()

        # If a == b, then this reduces to the special case of the
        # Wiener-Khinchin theorem (autocorrelation of a):

        # corr(a,a) = ifft(abs(fft(a))**2)
        # where the power spectrum of a is simply:
        # fft(corr(a,a)) == abs(fft(a))**2
        ############################################################
        # in this function, numpy.correlate is used to calculate the VAF

        self.VAF2 = np.zeros((self.Niter-1)*2 - 1)
        for i in range(self.Nions):
            for j in range(3):
                self.VAF2 += np.correlate(self.velocity[:,i,j],
                                          self.velocity[:,i,j], 
                                          'full')
        # two-sided VAF
        self.VAF2 /=  np.sum(self.velocity**2)
        self.VAF = self.VAF2[self.Niter-2:]
        
    def get_partial_VAF(self,atom_index):
        """ Velocity Autocorrelation Function """

        # VAF definitions
        # VAF(t) = Natoms^-1 * \sum_i <V_i(0) V_i(t)>
        ############################################################
        # Fast Fourier Transform Method to calculate VAF
        ############################################################
        # The cross-correlation theorem for the two-sided correlation:
        # corr(a,b) = ifft(fft(a)*fft(b).conj()

        # If a == b, then this reduces to the special case of the
        # Wiener-Khinchin theorem (autocorrelation of a):

        # corr(a,a) = ifft(abs(fft(a))**2)
        # where the power spectrum of a is simply:
        # fft(corr(a,a)) == abs(fft(a))**2
        ############################################################
        # in this function, numpy.correlate is used to calculate the VAF

        self.partial_VAF2 = np.zeros((self.Niter-1)*2 - 1)
        for i in range(len(atom_index)):
            for j in range(3):
                self.partial_VAF2 += np.correlate(self.partial_velocity[:,i,j],
                                          self.partial_velocity[:,i,j], 
                                          'full')
        # two-sided VAF
        self.partial_VAF2 /=  np.sum(self.partial_velocity**2)
        self.partial_VAF = self.partial_VAF2[self.Niter-2:]
        print(self.partial_VAF.shape)
    
    def choose_window(self,data):
        sigma = 2 * math.sqrt(2 * math.log(2))

        window_function = signal.gaussian(len(data), self.std/sigma, sym=False)
        
        return window_function
    def calc_ACF(self,array_1D):
    # Normalization
        yunbiased = array_1D - np.mean(array_1D, axis=0)
        ynorm = np.sum(np.power(yunbiased,2), axis=0)
    #    print("the average value of input data array", ynorm)
        
        autocor = signal.fftconvolve(array_1D,
                                     array_1D[::-1],
                                     mode='full')[len(array_1D)-1:] / ynorm
        return autocor

    
    def calc_FFT(self,array_1D, window):
        """
        This function is for calculating the "intensity" of the ACF at each frequency
        by using the discrete fast Fourier transform.
        """
        ####
        #### http://stackoverflow.com/questions/20165193/fft-normalization
        ####
        # window = self.choose_window(array_1D)
        WE = sum(window) / len(array_1D)
        wf = window / WE
        # convolve the blackman-harris window function.
        
        sig = array_1D * wf
        	
        # A series of number of zeros will be padded to the end of the \
        # VACF array before FFT.
        N = zero_padding(sig)
        
        yfft = np.fft.fft(sig, N, axis=0) / len(sig)
        #    yfft = np.fft.fft(data, n=int(N_fft), axis=0)/len(data) # no window func.
        print("shape of yfft {:}".format(np.shape(yfft)))
        return np.square(np.absolute(yfft))

    def phononDos(self, unit='THz', sigma=5):
        """ Phonon DOS from VAF """

        N = self.Niter - 1
        # Frequency in THz
        omega = fftfreq(2*N-1, self.potim) * 1E3
        # Frequency in cm^-1
        if unit.lower() == 'cm-1':
            omega *= 33.35640951981521
        if unit.lower() == 'mev':
            omega *= 4.13567
        # from scipy.ndimage.filters import  gaussian_filter1d as gaussian
        # smVAF = gaussian(self.VAF2, sigma=sigma)
        # pdos = np.abs(fft(smVAF))**2
        if self.VAF2 is None:
            self.getVAF()
        pdos = np.abs(fft(self.VAF2 - np.average(self.VAF2)))**2

        return omega[:N], pdos[:N]
    def partial_phononDos(self, unit='THz', sigma=5):
        """ Phonon DOS from VAF """

        N = self.Niter - 1
        # Frequency in THz
        omega = fftfreq(2*N-1, self.potim) * 1E3
        # Frequency in cm^-1
        if unit.lower() == 'cm-1':
            omega *= 33.35640951981521
        if unit.lower() == 'mev':
            omega *= 4.13567
        # from scipy.ndimage.filters import  gaussian_filter1d as gaussian
        # smVAF = gaussian(self.VAF2, sigma=sigma)
        # pdos = np.abs(fft(smVAF))**2
        if self.partial_VAF2 is None:
            self.get_partial_VAF()
        pdos = np.abs(fft(self.partial_VAF2 - np.average(self.partial_VAF2)))**2
        print(pdos.shape)
        # from scipy.ndimage.filters import  gaussian_filter1d as gaussian
        # smVAF = gaussian(self.partial_VAF2, sigma=sigma)
        # pdos = np.abs(fft(smVAF))**2

        return omega[:N], pdos[:N]

    def PCF(self, bins=50, Niter=10, A='', B=''):
        """ Pair Correlation Function """

        if not A:
            A = self.TypeName[0]
        if not B:
            B = A

        whichA = self.ChemSymb == A
        whichB = self.ChemSymb == B
        indexA = np.arange(self.Nions)[whichA]
        indexB = np.arange(self.Nions)[whichB]
        posA = self.position[:,whichA,:]
        posB = self.position[:,whichB,:]

        steps = range(0, self.Niter, Niter)
        rABs = np.array([posA[i,k,:]-posB[i,j,:]
                         for k in range(indexA.size)
                         for j in range(indexB.size)
                         for i in steps
                         if indexA[k] != indexB[j]])
        # periodic boundary condition
        rABs[rABs > 0.5] -= 1.0
        rABs[rABs <-0.5] += 1.0
        # from direct to cartesian coordinate
        rABs = np.linalg.norm(np.dot(self.cell, rABs.T), axis=0)
        # histogram of pair distances
        val, b = np.histogram(rABs, bins=bins)
        # density of the system
        rho = self.Nions / np.linalg.det(self.cell)
        # Number of A type atom
        Na = self.Nelem[self.TypeName.index(A)]
        # Number of B type atom
        Nb = self.Nelem[self.TypeName.index(B)]
        dr = b[1] - b[0]
        val = val * self.Nions / (4*np.pi*b[1:]**2 * dr) / (Na * Nb * rho) / len(steps)

        return val, b[1:]

################################################################################      
delta_t = 0.5 * 1e-15
def visualization(derivative, ACF, wavenumber, intensity):
    derivative = derivative * delta_t
    plt.subplot(3,1,1)
    L1 = np.arange(len(derivative))
    plt.plot(L1, derivative, color="red", linewidth=1.5)
    plt.axis([0, len(derivative), 
              1.1 * np.min(derivative),
              1.1 * np.max(derivative)],
             )
    plt.xlabel("Data Points", fontsize=15)
    plt.ylabel("Derivative of Dipole (a.u.)", fontsize=15)

    plt.subplot(3,1,2)
    L2 = np.arange(len(ACF))
    plt.plot(L2, ACF, color='red', linewidth=1.5)
    plt.axis([0, len(ACF),
              1.1 * np.min(ACF),
              1.1 * np.max(ACF)],
             )
    plt.xlabel("Data Points", fontsize=15)
    plt.ylabel("VACF (a.u.)", fontsize=15)

    plt.subplot(3,1,3)
    plt.plot(wavenumber, intensity, color="black", linewidth=1.5)
    plt.axis([0, 4000, 0, 0.001]
             )
    plt.xlabel("Wavenumber (cm$^{-1}$)", fontsize=15)
    plt.ylabel("Intensity (a.u.)", fontsize=15)
    plt.subplots_adjust(hspace = 0.5)
    plt.show()

def calc_derivative(array_1D, delta_t):
    ''' The derivatives of the angle_array were obtained by using the
    finite differences method.
    '''
    dy = np.gradient(array_1D)
    return np.divide(dy, delta_t)
# test code of the above class
if __name__ == '__main__':
    inp = xdatcar()
    inp.read_partial_xdat([i for i in range(144)])
    a,b = inp.get_partial_diffusion_coefficient()
    print(inp.partial_atoms_number)
    print(a)
    print(np.sum(b))
    # inp.get_partial_VAF([65,92,93])
#     coord = inp.positionC
#     window = inp.choose_window(coord)
#     for i in range(len(coord[0,:,:])):
#         normal_vectors = np.linalg.norm(coord, axis=-1)
#         derivative = calc_derivative(normal_vectors[:,i], delta_t)
#         ACF = inp.calc_ACF(derivative)
        
#         yfft_i = inp.calc_FFT(ACF, window)
        
#         if i == 0:
#              yfft = yfft_i
#         else:
#              yfft += yfft_i
#     # wavenumber = np.fft.fftfreq(len(yfft), delta_t * c)[0:int(len(yfft) / 2)]
# #    wavenumber = wavenumber * scaling_factor
#     wavenumber = np.fft.fftfreq(len(yfft), delta_t * c)[0:int(len(yfft) / 2)]
# #    wavenumber = wavenumber * scaling_factor

#     intensity = yfft[0:int(len(yfft)/2)]
#     visualization(derivative, ACF, wavenumber, intensity)

    # plt.plot((np.abs(fft(inp.VAF[inp.Niter-2:]))**2))
    # print inp.VAF.shape
    # plt.plot(inp.Time, inp.VAF, 'ko-', lw=1.0, ms=2,
    #         markeredgecolor='r', markerfacecolor='red')
    # 
    # plt.xlabel('Time [fs]')
    # plt.ylabel('Velocity Autocorrelation Function')

    x, y = inp.partial_phononDos('cm-1')
    plt.plot(x, y, 'ko-')
    plt.xlim(0, 5000)
    plt.ylim(-0.5, 500.0)

    # val, b = inp.PCF(100, 1)
    # plt.plot(b, val)
    # plt.axhline(y=1, color='r')
    # plt.xlim(0, 5)
    # plt.show()





