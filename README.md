> [!CAUTION]
> This code is deprecated! If you want to construct GCH, please use `DirectionalConvexHull` in [skmatter](https://github.com/scikit-learn-contrib/scikit-matter)

------------------------------

The Generalised Convex Hull suite contains 4 scripts written in python, which depend on :

# DEPENDENCIES

- ASE (for the generation of the "rattled" structures
- scipy
- numpy 
- scikit-learn (only for a quick kernel PCA)
- multiprocessing (for Parallel gch : to be added to v2.0)
- random 
- argparse
- json
- shutil (can be removed, down the line)

# HOW TO USE

The current implementations expects the user to be able to produce autonomously similarity kernels
between structural datasets (A.xyz and B.xyz).
To create a generalised convex hull for a structure search set "Set.xyz", first feed the gch_init.py
an extended xyz formatted dataset library (check ASE.io documentation for correct formatting), a kernel
measure of pairwise similarity of the elements in Set.xyz (Set.ker). It is expected to be a square matrix
with N x N elements, with N being the number of structures at hand.
Finally, pass to gch_init.py is the list of energies/atom for each entry in the set, and the name of the directory
you will keep the (several) outputs in.
The remaining options determine the sensitivity of the construction and are hyperparameters explained in [ref].


Once you have generated the folder "wdir", check inside and find the shaketraj.xyz, you will need that to compute 
the last element needed to run the gch with : shaketraj.ker .
Use your kernel calculator to produce a kernel matrix of similarity between the nshaken*nref rattled structures and your
set structures, it will have (nshaken*nref) * N entries.

Pass only the wdir and the freshly calculated shaketraj.ker to the gch_run.py to start the gch sampling.
The outputs "vprobprune.dat" and "vlist.dat"  will be stored in wdir, and will contain respectively :

vproborune.dat : List of probabilities of each structures of being vertex at every iteration of the pruning (will contain N*n_pruning_steps elements)
vlist 	       : List of vertices found when all vertices have at least 50% chance of being a vertex 


# IN THE FOLDER

You will find, besides gch_init and gch_run , also lib_gch and gch_utils. They contain functions useful for the GCH framework constructions,
that you may want to explore in case you're feeling adventurous.

gch_utils : A set of simple, non gch related functions that we use throughout the code
	-FPS : finds the indices of the farthest point sampled structures in your set, given a kernel
	-skenter : uses scikit learn to center a kernel matrix 
	-kpca : builds a kernel principal component analyisis and extracts the kpca vectors
	-ookpca : embeds external points onto the reference kernel and gives their kpca vectors
	-extractsubm : extracts the subkernel relative to the plist indices


lib_gch : A set of interconnected functions that build the GCH framework	
	-get_refgch : build a single convex hull made on the cols columns of pfile, using 0 as energy. It returns the distance of each point from the hull and it's used as a preliminary screening step to remove excessively unstable compounds.
	-get_gch    : builds a quick convex hull and finds the vertex list and associates sigma_s.Used to perform the sampling during the pruning iterations.
	-eval_sampled_sigmaKPCA : uses the shaken kernel coordinates in the kPCA space to estimate how the cartesian uncertainty reflects on the kernel space. 
	-estimate_residual_sigmaE : it uses a ridge regression scheme to quantify the weight that every kpca component has on the total energy of the system


	-create_samples_sigmaKPCA : It creates a list of shaken structures from the nref selected during the initalisation procedure and saves them in shaketraj.xyz in the wdir folder.
	- initialize_random_sample_GCH : It selects randomly the nref points to build the shaketraj with and creates the folder and the input.json files needed to run the gch_run on.
	- initialize_fps_sample_GCH : Same as before, but finds the structures using farthest point sampling

	-sample_GCH : repeats Nsamples the GCH constructions to extract each points probability of becoming a vertex 
	-parallel_sample_GCH : sames as above, but parallelized over the threads. CURRENTLY NOT WORKING
	-prune_GCH : Continues calling Sample_GCH until minimum probability per vertex reaches the requires number. 
	-parallel_prune_GCH : Same as above, but parallelized over the threads . CURRENTLY NOT WORKING  
