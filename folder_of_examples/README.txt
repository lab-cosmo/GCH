In this examples folder, you should find :

1) A set of example structures as a zip file containing the individual structure res files: 'example_set.zip'
2) An extended xyz file containing all the example structures: 'example_set.xyz'
3) A file containing an example NxN SOAP kernel matrix for the example set: 'example_kernel.ker'
4) A file containing an ordered list of all the enegies per atom for each structure in the example set (in meV): 'example_energies.txt'
5) A further folder containing example output files from the script 'gch_init.py' along with an explanation file: 'example_outputs'
6) A file containing an example ((Nshaken + 1)* Nref) x N SOAP kernel matrix comparing the example shaken structures to the original example set: 'example_shaken_kernel.ker'
7) A Python script 'kernel_generation.py' designed to show how the original kernel, structure sets, and energy files were generated. This is not a finalised generalisable script but is set up so if run, it will generate a set of test 
   files like those in this examples folder
8) A Python script 'shaken_kernel_generation.py' designed to show how the shaken kernel file was generated. Like 7, This is not a finalised generalisable script but is set up so if run, it will generate a test shaken kernel like that in    this example folder


USING THE EXAMPLE FILES
---------------------------------------------------------------------

Running gch_init.py : Files 2-4 are required when running gch_init.py

Running gch_run.py : File 6 is needed to run gch_run.py

Starting in a directory containing atleast the files 'gch_init.py', 'gch_run.py', 'gch_utils.py', and 'lib_gch.py' and using a Python environment set up as outlined on the GCH github:



Copy files 2-4 and file 6 into your starting directory and run the command: 

python3 gch_init.py example_kernel.ker -ixyz example_set.xyz -nrg example_energies.txt --nshake 10 --nref 10 -wdir [output folder name] 

 ([output folder name] is the name you wish to give to your output folder. This can be whatever you choose, but must not already be an existing sub-directory.)


When this has finished, go to the output folder to check that it is generating the expected files - it should have the appearance of 5) 

Return to your starting directory and run the command: python gch_run.py example_shaken_kernel.ker -wdir [output folder name]
 ([output folder name] must be the same as before)

When this has run, enter the output folder and read 'vlist.idx' this should give a list of the indices of structures selected as GCH vertices


NOTE
--------------------------------------------------------------------

Running gch_init.py using these examples will generate its own set of outputs. Due to the different random selection and rattling of reference structures,
these examples will differ from the example outputs contained here. As such, running the two scripts with the example files as outlined here is not a test of
a single run of the process from start to finish, but rather individual tests of the two scripts.

