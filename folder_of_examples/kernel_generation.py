import numpy as np
import ase 
import ase.io as aseio
import json
from ase.spacegroup import crystal
from ase.symbols import Symbols

import rascal
from rascal.representations import SphericalInvariants as SOAP
from rascal.neighbourlist.structure_manager import (
                mask_center_atoms_by_species)

import sys
from skcosmo.decomposition import PCovR
from sklearn.decomposition import PCA
from skcosmo.preprocessing import StandardFlexibleScaler
from zipfile import ZipFile
import os





#establish lists to hold the atoms objects and the energies
struct=[]
atom_energy=[]

#establish variables needed to convert from kJ/mol structure energies to a.u energy per atom
avagadro = 6.0221409E23
molecule_atoms=24
conversion = 6.242E24

########################################################################################################################################

#This section is used to convert our group's style res files into atoms objects, create the extended xyz file and get the atom energies

#########################################################################################################################################
with ZipFile('example_set.zip','r') as structures_zip:
    #get list of filenames
    structure_list = structures_zip.namelist()


    for structure in structure_list:
    #extract the spacegroup info from the filename
       spacegroup = structure.split('-')
       spacegrp = int(spacegroup[2])

    #unzip and open the file to be worked with
       structures_zip.extract(structure)
       file = open(structure, "r")

    #read the file and edit the lines to add zeros columns and remove the numbers from the element coloumn
       lines = file.readlines()
    #Read in total energy from top line, convert to meV per atom and add to list
       energy_line = lines[0].split()
       total_energy = energy_line[2]
       atoms_energy = (float(total_energy)/((avagadro)*(molecule_atoms)))*conversion
       atom_energy.append(atoms_energy)
    #identify where co-ordinate part begins-i.e which lines to edit
       search = 'SFAC'
       one_less = [lines.index(line) for line in lines if search in line]
       start = one_less[0] + 1
    #looping over coordinate section of lines to be edited
       for i in range(start,(len(lines))):
    
           col_list = (lines[i]).split()
        #take only the first character from column 1 - to get rid of the numbers
           col_list[0]= (col_list[0])[0:1]
        #add two columns of zeros
           col_list.append(0)
           col_list.append(0)
        #make a complete line again
           str_col_list = [str(element) for element in col_list]
           lines[i]= " ".join(str_col_list) + "\n"
       #write the new edited file content back into the original file
       file = open(structure, "w")
       file.writelines(lines)
       file.close()
       

       #Read the structure into ase to get atoms object
       atom_struc = aseio.read(structure)
       #add the spacegroup info and use crystal structure to apply bulk/crystal info to the releavant atoms object
       atom_struc = crystal(symbols=atom_struc,spacegroup=spacegrp,pbc=True)
       #write atoms object to the list of structures and to the main file
       struct.append(atom_struc)
       ase.io.write('test_set.xyz',images=atom_struc,append=True)
       #remove the extracted file to save space
       os.remove(structure)

structures_zip.close()

#Save the per atom energies to file
energy_array = np.array(atom_energy)
np.savetxt('test_energies.txt',energy_array)

####################################################################################################################################

#This section is used to get SOAP descriptors and calculate the kernel

####################################################################################################################################

#set up hyperparamaters
HYPERS = {'soap_type': 'PowerSpectrum','interaction_cutoff': 4.,'max_radial': 8,'max_angular': 6,'gaussian_sigma_constant': 0.3,
                             'gaussian_sigma_type': 'Constant',
                                 'cutoff_smooth_width': 0.5,
                                     'radial_basis': 'GTO',
                                         'inversion_symmetry': True,
                                             'normalize' : True
                                             }


#mask each atoms object in the list to know which atoms to do SOAP descriptors for
#wrap each object to ensure atoms are in unbit cell (this is for consistency with calculation of shaken kernels)

for s in struct:
    mask_center_atoms_by_species(s, ['C', 'H', 'O'])
    s.wrap(eps=1e-18)


#set up 'calculator' for the type of soap descriptors wanted
features = SOAP(**HYPERS)

#Transform the atoms objects to get the SOAP descriptors
struct = features.transform(struct)

#Set up 'calculator' for the kernel
kernel = rascal.models.Kernel(features, kernel_type='Full', target_type='Structure',zeta=1)


#Call the kernel to get the similarity matrix
results=kernel.__call__(struct)

#save unnormalised kernel as back-up
np.save('test_unnormalised_kernel.npy',results)
np.savetxt('test_unnormalised_kernel.ker',results)

#Normalise the kernel
normalised_results=results.copy()
for i in range(0,results.shape[0]):
    for j in range(0,results.shape[0]):
        normalised_results[i,j]=results[i,j]/((results[i,i]*results[j,j])**0.5)

#Save the kernel
np.save('test_normalised_kernel.npy',normalised_results)
np.savetxt('test_normalised_kernel.ker',normalised_results)
print(normalised_results)
