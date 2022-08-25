from ase import io as aseio
from rascal.representations import SphericalInvariants as SOAP
from rascal.neighbourlist.structure_manager import (
                        mask_center_atoms_by_species)

import sys
from skcosmo.decomposition import PCovR
from sklearn.decomposition import PCA
from skcosmo.preprocessing import StandardFlexibleScaler
from ase.spacegroup import crystal
from ase.symbols import Symbols
import rascal
import numpy as np


#Read in necessaary structure sets
struct= aseio.read('example_set.xyz', index = ':')
struct_2 = aseio.read('example_outputs/shaketraj.xyz',index=':')

#fix cell wrapping issue
for shakenstruc in struct_2:
    shakenstruc.wrap(eps=1e-18)
for normstruc in struct:
    normstruc.wrap(eps=1e-18)


#mask all structure sets as necessary
for s in struct:
    mask_center_atoms_by_species(s, ['C', 'H', 'O'])
for s in struct_2:
    mask_center_atoms_by_species(s,['C','H','O'])

#set up SOAP calculator and kernel calculators
HYPERS = {'soap_type': 'PowerSpectrum','interaction_cutoff': 4.,'max_radial': 8,'max_angular': 6,'gaussian_sigma_constant': 0.3,
                                     'gaussian_sigma_type': 'Constant',
                                        'cutoff_smooth_width': 0.5,
                                        'radial_basis': 'GTO',
                                        'inversion_symmetry': True,
                                        'normalize' : True
 
                                        }


features = SOAP(**HYPERS)
kernel = rascal.models.Kernel(features, kernel_type='Full', target_type='Structure',zeta=1)


#Transform the structure sets
struct = features.transform(struct)
struct_2 = features.transform(struct_2)


#Generate unnormalised shaken kernel and save
unnormalised_shaken=kernel.__call__(struct_2,struct)

np.save('unnormalised_shaken_test.npy',unnormalised_shaken)
np.savetxt('unnormalised_shaken_test.ker',unnormalised_shaken)

#Regenerate unnormalised original kernel and save --y=You only need diagonal entries of this but for now I've re-generated all of them  
unnormalised_original= kernel.__call__(struct)

np.save('unnormalised_original_test.npy',unnormalised_original)
np.savetxt('unnormalised_original_test.ker',unnormalised_original)

#get rid of transformed original structures to save memory
struct=[]

#Generate shaken-shaken  kernel and save -- You only need diagonal entries of this but for now I've generated all
unnormalised_shaken_shaken = kernel.__call__(struct_2)
np.save('unnormalised_shaken_shaken_test.npy',unnormalised_shaken_shaken)
np.savetxt('unnormalised_shaken_shaken_test.ker',unnormalised_shaken_shaken)


#Get rid of transformed shaken structures to save memory
struct_2=[]


#Normalise the kernel, by dividing shaken by the necessary other kernel entries
normalised_shaken=unnormalised_shaken.copy()
for i in range(0,unnormalised_shaken.shape[0]):
    for j in range(0,unnormalised_shaken.shape[1]):
        normalised_shaken[i,j] = unnormalised_shaken[i,j]/((unnormalised_shaken_shaken[i,i]*unnormalised_original[j,j])**0.5)

#save the kernel
np.save('normalised_shaken_kernel_test.npy',normalised_shaken)
np.savetxt('normalised_shaken_kernel_test.ker',normalised_shaken)
print(normalised_shaken)
