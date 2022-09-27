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
import os


########################################################################################################################################

#This section is used to read the xyz file into a list of ase atoms objects

#########################################################################################################################################
struct = aseio.read('example_set.xyz',index=':')

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
