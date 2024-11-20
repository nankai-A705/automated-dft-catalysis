import numpy as np
import pandas as pd
from lxml import etree
from time import time
import os
import xml.etree.ElementTree as ET

class Dos:
    def __init__(self,xml):
        '''

        :param xml:
        '''
        if xml == None:
            xml = 'vasprun.xml'
        
        self.file = 'vasprun.xml'
        
        self.xml = etree.parse(xml)
            
        self.efermi = self.get_efermi()
        self.spin = self.get_spin()
        self.symbols, self.atomtype = self.get_symbols_and_atomtype()

    def get_efermi(self):
        '''
        get efermi energy
        :return: float
        '''
        fermi_path = '//calculation/dos/i[@name="efermi"]/text()'
        efermi = self.xml.xpath(fermi_path)
        efermi = float(efermi[0])
        return efermi

    def get_spin(self):
        '''
        to get spin from the xml file
        :return: int
        '''
        spin_path = '///incar/i[@name="ISPIN"]/text()'
        spin = int(self.xml.xpath(spin_path)[0].strip())
        return spin

    def get_symbols_and_atomtype(self):
        '''
        e.g. symbols = ['Na','Na','Mn','Mn','Mn','O',O']    : types of atoms in order
            atomtype = ['Na','Mn','O']
        :return:
        '''
        symbol_path = '//atominfo/array[@name="atoms"]/set/rc/c/text()'
        symbols = self.xml.xpath(symbol_path)
        symbols = [s.strip() for s in symbols if s.strip().isalpha()]
        atomtype = []
        for symbol in symbols:
            if symbol not in atomtype:
                atomtype.append(symbol)
        return symbols, atomtype
    def get_index_range(self, atom=''):
        '''
        get the index range of the target atom in self.symbols
        :param atom:
        :return: index range of a specific atom like (1,16)
        '''
        first_index = self.symbols.index(atom) + 1
        last_index = len(self.symbols) - self.symbols[::-1].index(atom)
        index_range = (first_index, last_index)
        return index_range

    def get_tdos(self, spin=None):
        '''
        to calculate the total DOS
        :return: a dataframe with two colunms ['energy', 'total']
        '''
        if spin == None:
            spin = self.spin
        total_up_path = '//calculation/dos/total/array/set/set[@comment="spin 1"]/r/text()'
        total_down_path = '//calculation/dos/total/array/set/set[@comment="spin 2"]/r/text()'
        # spin up
        tdos_up_datas = self.xml.xpath(total_up_path)
        tdos_up_array = np.array([list(map(float, s.split())) for s in tdos_up_datas])
        tdos_up_df = pd.DataFrame(tdos_up_array[:,:2],columns=['energy','total'])
        tdos = tdos_up_df
        # spin down
        if spin == 2:
            tdos_down_datas = self.xml.xpath(total_down_path)
            tdos_down_array = np.array([list(map(float, s.split())) for s in tdos_down_datas])
            tdos_down_df = pd.DataFrame(tdos_down_array[:, :2], columns=['energy', 'total'])
            tdos_down_df['total'] *= (-1)

            # spin up + spin down
            tdos = pd.concat([tdos_up_df, tdos_down_df], ignore_index=True)

        tdos['energy'] -= self.efermi
        tdos_up_df['energy'] -= self.efermi
        tdos_down_df['energy'] -= self.efermi
        
        

        
        return tdos, tdos_up_df, tdos_down_df

    def get_one_ion_pdos(self, ion_index, spin, orbitals=[]):
        '''
        to obtain the density of states for a specific atomic orbital and spin
        :param ion_index: index of a specific ion
        :param spin:
        :param orbitals: which orbitals are needed
        :return: a dataframe with column names ['energy']+orbital
        '''
        one_ion_path = f'//calculation/dos/partial/array/set/set[@comment="ion {ion_index}"]/set[@comment="spin {spin}"]/r/text()'
        ion_dos_datas = self.xml.xpath(one_ion_path)
        ion_dos_array = np.array([list(map(float, s.split())) for s in ion_dos_datas])
        column_names = self.get_column_names()
        ion_dos_df = pd.DataFrame(ion_dos_array,columns=column_names)
        if 'p' in orbitals:
            ion_dos_df['p'] = ion_dos_df['py'] + ion_dos_df['pz'] + ion_dos_df['px']
        if 'd' in orbitals:
            ion_dos_df['d'] = ion_dos_df['dxy'] + ion_dos_df['dyz'] + ion_dos_df['dz2'] + ion_dos_df['dxz'] + ion_dos_df['x2-y2']
        if 'f' in orbitals:
            ion_dos_df['f'] = ion_dos_df['f1'] + ion_dos_df['f2'] + ion_dos_df['f3'] + ion_dos_df['f4'] + ion_dos_df['f5'] + ion_dos_df['f6'] + ion_dos_df['f7']

        # only choose specific orbitals we want
        ion_dos_df = ion_dos_df[['energy']+orbitals]
        ion_dos_df['energy'] -= self.efermi
        if spin == 2:
            for orbital in orbitals:
                ion_dos_df[orbital] *= (-1)
        return ion_dos_df

    def get_column_names(self):
        '''
        to obtain the calculated orbital names from the XML file
        :return: a list
        '''
        column_name_path = '//calculation/dos/partial/array/field/text()'
        column_names = self.xml.xpath(column_name_path)
        column_names = [name.strip() for name in column_names]
        return column_names
    def get_system_mass(self):
        atom_info = '//atominfo/array[@name="atomtypes"]/set/rc/c/text()'
        atoms = self.xml.xpath(atom_info)
        total_type = len(set(self.atomtype))
        c_number = int((len(atoms)/total_type))
        mass = 0
        
        for i in range(total_type):
            mass += float(atoms[c_number*i]) * float(atoms[c_number*i+2])
            
        
        return  mass
    
        

    def get_pdos(self, atom='', spin=None, orbitals=[]):
        '''
        to calculate the Projected DOS
        :param atom: e.g. 'Na'
        :param spin:
        :param orbitals: which orbitals are needed
        :return:
        '''
        if spin == None:
            spin = self.spin
        index = self.get_index_range(atom)
        pdos = pd.DataFrame()
        up_pdos = pd.DataFrame()
        for i in range(index[0],index[1]+1):
            one_up_pdos = self.get_one_ion_pdos(ion_index=i, spin=1, orbitals=orbitals)
            if i == index[0]:
                up_pdos = one_up_pdos
                continue
            for orbital in orbitals:
                up_pdos[orbital] += one_up_pdos[orbital]
        if spin == 2:
            down_pdos = pd.DataFrame()
            for i in range(index[0],index[1]+1):
                one_down_pdos = self.get_one_ion_pdos(ion_index=i, spin=2, orbitals=orbitals)
                if i == index[0]:
                    down_pdos = one_down_pdos
                    continue
                for orbital in orbitals:
                    down_pdos[orbital] += one_down_pdos[orbital]
            pdos = pd.concat([up_pdos, down_pdos], ignore_index=True)
            return pdos,up_pdos,down_pdos
        elif spin == 1:
            pdos = up_pdos
            return [pdos]
    def get_d_pdos(self, atom='', spin=None):
        '''
        to get d-projected dos of a specific atom
        :param atom: chemical symbol like 'Mn','Na'
        :param spin:
        :return:
        '''
        if spin == None:
            spin = self.spin
        orbital = ['d']
        if spin == 2:
            d_pdos, up_pdos, down_pdos = self.get_pdos(atom=atom, spin=spin, orbitals=orbital)
            return d_pdos, up_pdos, down_pdos
        elif spin == 1:
            d_pdos = self.get_pdos(atom=atom, spin=spin, orbitals=orbital)
            return d_pdos
        
    def get_basis_vectors(self):
        with open(os.path.join('.',
                           self.file)) as f:
            tree = ET.parse(f)
            
        
        
            basis_set_path = './/structure[@name="finalpos"]/crystal/varray[@name="basis"]/v'
            basis = tree.findall(basis_set_path)
            print(basis)
            basis_vectors = []
            for b in basis:
            
                data = b.text.strip()
                
                basis_vector = [float(val) for val in data.split()]
                basis_vectors.append(basis_vector)
            
        return basis_vectors
        
        
    def get_cross_section_area(self):
        
        basis_sets = self.get_basis_vectors()
        area = abs(basis_sets[0][0]*basis_sets[1][1] - basis_sets[0][1]*basis_sets[1][0])
        
        return area
        
        
    def write_to_excel(self, filename, dos_types={}):
        '''
        to write the output results to an Excel file
        :param filename:
        :param dos_types: key: sheet_name; value:a dataframe of dos results
        :return:
        '''
        excel_writer = pd.ExcelWriter(filename, engine='openpyxl')
        for atom, dos_type in dos_types.items():
            dos_type.to_excel(excel_writer, sheet_name=f'{atom}', index=False)
        excel_writer.close()

    def smooth(self, y, box_pts):
        """
        Smooth the noisy density of state distribution using convolution operater.
        Args:
           y (array): one-dimensional input array
           box_pts (int): average box size parameter
        Returns:
           array: one-dimentional array of smoothed density of state values
        """
        box = np.ones(box_pts) / box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

    def moment(self,x, y, n):
        """
        moments function to calculate the porbability distribution characteristics of density of states.
        Args:
           x (array): energy values
           y (array): density of states
           n (int): order parameter of moments function
        Returns:
           float: moment descriptor
        """
        p = x ** n * y
        return np.trapz(p, x) / np.trapz(y, x)

    def density_moments(self, dos_df, mode=None):
        '''

        :param dos_df: a dataframe of dos
        :param mode:
        :return: list: moment characteristics including filling, center, sigma, skewness, kurtosis
        '''

        energies = dos_df['energy'].values
        dos = dos_df['d'].values
        # smooth the noisy density state
        dos_rev = self.smooth(dos[:], 15)  # smoothing function
        # determine the index of first non-positive value
        Ind = np.argmax(energies > 0)
        # calculate the moment descriptors for the density of states
        if mode is None:
            filling = np.trapz(dos_rev[0:Ind], energies[0:Ind]) / np.trapz(dos_rev, energies)
            center = self.moment(energies[:], dos_rev[:], 1)
            sigma_c = np.sqrt(self.moment(energies[:] - center, dos_rev[:], 2))
            skewness = self.moment(energies[:] - center, dos_rev[:], 3) / sigma_c ** 3
            kurtosis = self.moment(energies[:] - center, dos_rev[:], 4) / sigma_c ** 4
            return [filling, center, sigma_c, skewness, kurtosis]
        if mode == 'center':
            return [self.moment(energies[:], dos_rev[:], 1)]
        


orbitals_all = ["s","py","pz","px","dxy","dyz","dz2","dxz","x2-y2", "f1","f2","f3","f4","f5","f6","f7"]
orbitals_s = ['s']
orbitals_p = ["py","pz","px"]
orbitals_d = ["dxy","dyz","dz2","dxz","x2-y2"]
orbitals_f = ["f1","f2","f3","f4","f5","f6","f7"]
orbitals = ["s","p","d"]

dos = Dos('vasprun.xml')
fermi = dos.get_system_mass()
print(dos.get_basis_vectors())
# print(fermi)
# a,b,c = dos.get_tdos(spin=2)

# df1 = pd.DataFrame(b)

# e = df1.loc[:,'energy']
# f = df1.loc[:,'total']
# import matplotlib.pyplot as plt
# plt.plot(e,f)
# plt.show()


# e.g. 1
'''
dos_types = {}
dos_types['tdos'] = dos.get_tdos()
for atom in dos.atomtype:
    dos_types[atom] = dos.get_pdos(atom=atom, orbitals=orbitals)[0]
dos.write_to_excel(filename='test6.xlsx', dos_types=dos_types)
'''
# e.g. 2
'''
dos_types = {'Mn':['d','dxz'],'O':['p']}
for atom, orbitals in dos_types.items():
    dos_types[atom] = dos.get_pdos(atom=atom, spin=2, orbitals=orbitals)[0]
dos.write_to_excel(filename='test3.xlsx', dos_types=dos_types)
'''
# e.g. 3 get d band center of Mn
'''
d_dos, d_up_dos, d_down_dos = dos.get_d_pdos(atom='Mn')
dos_types = {'d_dos':d_dos,'d_up':d_up_dos,'d_down':d_down_dos}
dos.write_to_excel(filename='d_band_test.xlsx',dos_types=dos_types)
d_up_center = dos.density_moments(dos_df=d_up_dos, mode='center')
d_down_center = dos.density_moments(dos_df=d_down_dos, mode='center')
'''

# def get_basis_vectors(file):
#     with open(os.path.join('.',
#                        file)) as f:
#         tree = ET.parse(f)
        
    
    
#         basis_set_path = './/structure[@name="finalpos"]/crystal/varray[@name="basis"]/v'
#         basis = tree.findall(basis_set_path)
#         print(basis)
#         basis_vectors = []
#         for b in basis:
        
#             data = b.text.strip()
            
#             basis_vector = [float(val) for val in data.split()]
#             basis_vectors.append(basis_vector)
        
#     return basis_vectors
# a = get_basis_vectors('vasprun.xml')
# print(a)
a = dos.get_cross_section_area()
print(a)
# print(a)
# print(c)