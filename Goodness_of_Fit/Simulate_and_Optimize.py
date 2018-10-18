import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions_GOF

'''
Usage: python Simulate_and_Optimize.py

This is meant to be a general use script to run dadi to perform simulations and goodness
of fit tests for any model on an afs/jsfs with one to three populations. The user 
will have to edit information about their allele frequency spectrum and provide a 
custom model. The instructions are annotated below, with a #************** marking 
sections that will have to be edited.

This script must be in the same working directory as Optimize_Functions_GOF.py, which
contains all the functions necessary for generating simulations and optimizing the model.


General workflow:
 The user provides a model and the previously optimized parameters for their empirical 
 data. The model is fit using these parameters, and the resulting model SFS is used to
 generate a user-selected number of Poisson-sampled SFS (ie simulated data). If the SNPs
 are unlinked, these simulations represent parametric bootstraps. For each of 
 the simulated SFS, the optimization routine is performed and the best scoring replicate
 is saved. The important results for such replicates include the log-likelihood and 
 Pearson's chi-squared test statistic, which are used to generate a distribution of the
 simulated data values to compare the empirical values to. In addition, theta, the sum of 
 the SFS, and the optimized parameters are also saved (to calculate confidence intervals).
 
Outputs:
 For each simulation, there will be a log file (ex. Simulation_1.someModel.log.txt) 
 showing the optimization steps per replicate and a summary file 
 (ex. Simulation_1.someModel.optimized.txt) that has all the important information. 
 Here is an example of the output from a summary file, which will be in tab-delimited 
 format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params( )
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688

 The results are also summarized across all simulations, and the best-scoring replicate 
 for each is written to a tab-delimited file called Simulation_Results.txt. Here is an 
 example of the contents for this file:
 
 Simulation	Best_Replicate	log-likelihood	theta	sfs_sum	chi-squared	optimized_params
 1	Round_3_Replicate_2	-162.23	1205.59	1628.0	218.38	0.0333,0.0338,0.4525,0.9208
 2	Round_3_Replicate_2	-186.49	1079.25	1588.0	271.11	0.0431,0.0427,0.4401,1.3532
 3	Round_3_Replicate_3	-225.82	650.49	1565.0	560.88	0.0997,0.1129,0.1136,1.0299
 4	Round_3_Replicate_2	-361.47	1279.59	1567.0	1369.01	0.0101,0.0478,0.5306,4.452
 5	Round_3_Replicate_1	-221.84	73.7	1492.0	518.56	0.5149,0.8345,0.0214,23.2449


Notes/Caveats:
 The data are simulated using the fs.sample() method, which is equivalent to a 
 parametric boostrap only if SNPs are unlinked across loci. For ddRADseq data where a 
 single SNP is selected per locus, this is generally true, and this workflow is valid.

Citations:
 If you use these scripts for your work, please cite the following publications.
 
 For the general optimization routine:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 52455263.
    doi: 10.1111/mec.14266
    
 For the goodness of fit testing scripts:
 	Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., 
	Portik, D.M., Streicher, J.W., and Loader, S.P. Vanishing refuge? Testing the forest
	refuge hypothesis in coastal East Africa using genome-wide sequence data for seven
	amphibians. In Review, Molecular Ecology.
 	
-------------------------
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated September 2018
'''

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
snps = "/Users/dan/Dropbox/dadi_inputs/General_Script/dadi_2pops_North_South_snps.txt"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["North", "South"]

#**************
#projection sizes, in ALLELES not individuals
proj = [16,32]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'



#================================================================================
# Fit the empirical data based on prior optimization results, obtain model SFS
#================================================================================
''' 
 We will use a function from the Optimize_Functions_GOF.py script:
 	Optimize_Empirical(fs, pts, outfile, model_name, func, in_params, fs_folded)

Mandatory Arguments =
 	fs:  spectrum object name
 	pts: grid size for extrapolation, list of three values
 	outfile: prefix for output naming
 	model_name: a label to help name the output files; ex. "no_mig"
 	func: access the model function from within script
 	in_params: the previously optimized parameters to use
    fs_folded: A Boolean value indicating whether the empirical fs is folded (True) or not (False).
'''


#Let's start by defining our model
def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

#**************
#Make sure to define your extrapolation grid size.
pts = [50,60,70]

#**************
#Provide best optimized parameter set for empirical data.
#These will come from previous analyses you have already completed.
emp_params = [0.1487,0.1352,0.2477,0.1877]

#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True


#Fit the model using these parameters and return the folded or unfolded model SFS (scaled by theta).
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model, but
#everything else can stay as it is. See above for argument explanations.
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, pts, "Empirical", "sym_mig", sym_mig, emp_params, fs_folded=fs_folded)



#================================================================================
# Performing simulations and optimizing
#================================================================================
'''
 We will use a function from the Optimize_Functions_GOF.py script:
 	Perform_Sims(sim_number, model_fs, pts, model_name, func, rounds, param_number, 
 				reps=None, maxiters=None, folds=None)

Mandatory Arguments =
	sim_number: the number of simulations to perform
    model_fs:  the scaled model spectrum object name
    pts: grid size for extrapolation, list of three values
    model_name: a label to help label on the output files; ex. "no_mig"
    func: access the model function from within this script
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)

Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
'''

#**************
#Set the number of simulations to perform here. This should be ~100 or more.
sims = 100

#**************
#Enter the number of parameters found in the model to test.
p_num = 4

#**************
#Set the number of rounds here.
rounds = 3
#I strongly recommend defining the lists for optional arguments to control the settings 
#of the optimization routine for all the simulated data.
reps = [20,30,50]
maxiters = [5,10,20]
folds = [3,2,1]

#Execute the optimization routine for each of the simulated SFS.
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model, but
#everything else can stay as it is (as the actual values can be changed above).
Optimize_Functions_GOF.Perform_Sims(sims, scaled_fs, pts, "sym_mig", sym_mig, rounds, p_num, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds)

