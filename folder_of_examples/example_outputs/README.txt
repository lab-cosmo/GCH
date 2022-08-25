This folder contains example outputs from the script gch_init.py. It should contain:

1) 10 numbered 'shaketraj_i.xyz' files each containing the shaken versions of the corresponding reference structure as ase Atoms objects
2) 1 un-numbered 'shaketraj.xyz' file containing all the shaken versions of all reference structures as ase Atoms objects. This is the file normally used to generate the kernel between the shaken and original structure sets
3) A file 'input.json' containing the input parameters for gch construction and file paths for the input files. In an ordinary run, this will be read by gch_run.py in order to feed the files and paramaters needed into that second script.   In the case of running the examples, a new input.json file will be generated when running gch_init.py and it is that that will be read when testing gch_run.py 
4) A file 'refstruct.idx' containing the indices of the structures selected as reference structures 
5) A file 'labels.dat' containing the input paramaters in a seperate format to input.json
6) A file 'nrg-proj.npy' containing the input energies per atom and kPCA values for each structure. This data is used by gch_run.py to plot the GCH landscape and identify hull vertices. Again, as with 3)
   testing the examples will generate a new version of this file that will be used
