
"""
Created on Mon Jul 18 15:25:31 2022

@author: hxps
"""

from ase.db import connect
from ase.formula import Formula
slab = connect('RuNi-relaxed.db')
adsorbate = connect('top-adsorbates-relaxed.db')


    
def get_atom(name, idx):
    
    """
    Get Atom object according to id 
    
    """
    db = connect(name)
    atoms = db.get_atoms(id=idx)
    return atoms

def get_atomic_ratio(atoms,strip_elements,ratio_element):
    
    """
    Get Atomic ratio of atoms
    Usage: get_atomic_ratio(atoms,'H', 'Ni')
    Do not consider 'H', and return 'Ni' ratio in total atoms
    """
    formula = atoms.get_chemical_formula().strip(strip_elements)
    w = Formula(formula).count()
    new_dict = {k: v / total for total in (sum(w.values()),) for k, v in w.items()}
    
    return new_dict[ratio_element]

def get_adsorption_energy(adsorbates_db, slab_db, adsorbate_energy):
    
    for i, row in enumerate(adsorbates_db.select()):
        idx = row.label_id
    # delta_G = (row.energy-slab.get(id=idx).energy+3.34 + 
    #             row.helmholtz + row.zpe)
        delta_G = (row.energy-slab_db.get(id=idx).energy-adsorbate_energy)
        adsorbates_db.update(i+1, delta_G = delta_G)
        
   

# for i in range(115):
#     atoms = get_atom('top-adsorbates-relaxed.db', i+1)
#     ratio = get_atomic_ratio(atoms, 'H','Ni')
#     adsorbate.update(i+1, Ni_ratio = ratio)
    


# from ase.visualize import view
# atoms = get_atom('top-adsorbates-relaxed.db', 1)
# atoms.write('POSCAR')
# view(atoms)