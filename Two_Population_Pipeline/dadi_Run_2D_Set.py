import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions
import Models_2D

'''
Usage: python dadi_Run_2D_Set.py

This is a modified version of the 'dadi_Run_Optimizations.py' script in which
we run optimizations for 2D comparisons for a large set of models that have been
made available as part of published works. These models are stored in the
Models_2D.py script, and will be called directly here. The user can delete or
comment out models to analyze a subset of the models available. 

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  Models_2D.py script, which
has all the model definitions.


General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688


Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject. 

Citations:
 If you use these scripts or the main diversification models for your work, please
 cite the following publication:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266
 
 If you use the additional diversification models or the island models set please cite 
 the following publication:
    Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K.,
    Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. Sky, sea,
    and forest islands: diversification in the African leaf-folding frog
    Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). Early Access,
    Journal of Biogeography.
        
 If you are interesting in contributing your models to this workflow, please email me!

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
Updated May 2018
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
# Calling external 2D models from the Models_2D.py script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script for our optimization routines:
 
 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")
 
   Mandatory Arguments =
    fs:  spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to help label the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).

   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order

Below, I give all the necessary information to call each model available in the
Models_2D.py script. I have set the optimization routine to be the same for each
model using the optional lists below, which are included as optional arguments for
each model. This particular configuration will run 4 rounds as follows:
Round1 - 10 replicates, maxiter = 3, fold = 3
Round2 - 20 replicates, maxiter = 5, fold = 2
Round3 - 30 replicates, maxiter = 10, fold = 2
Round4 - 40 replicates, maxiter = 15, fold = 1

If this script was run as is, each model would be called and optimized sequentially;
this could take a very long time. For your actual analyses, I strongly recommend
creating multiple scripts with only a few models each and running them
independently. It is also not a good idea to mix models from the Diversification Set
and the Island Set, as each was meant to be mutually exclusive.

'''


#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
prefix = "_".join(pop_ids)

#**************
#make sure to define your extrapolation grid size (based on your projections)
pts = [50,60,70]

#**************
#Set the number of rounds here
rounds = 4

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]

#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True


'''
Diversification Model Set

This first set of models come from the following publication:

    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266

'''

# Split into two populations, no migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")


# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")


# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")


# Split with continuous symmetric migration, followed by isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")


# Split with continuous asymmetric migration, followed by isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_asym_mig", Models_2D.anc_asym_mig, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2")


# Split with no gene flow, followed by period of continuous symmetrical gene flow.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_sym_mig", Models_2D.sec_contact_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")


# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_asym_mig", Models_2D.sec_contact_asym_mig, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2")


# Split with no migration, then instantaneous size change with no migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig_size", Models_2D.no_mig_size, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, T1, T2")


# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig_size", Models_2D.sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")


# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig_size", Models_2D.asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")


# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig_size", Models_2D.anc_sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")


# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_asym_mig_size", Models_2D.anc_asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")


# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_sym_mig_size", Models_2D.sec_contact_sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")


# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_asym_mig_size", Models_2D.sec_contact_asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")

'''
The following 6 models were added to the Diversification Model Set by:

    Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K.,
    Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. Sky, sea,
    and forest islands: diversification in the African leaf-folding frog
    Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). Early Access,
    Journal of Biogeography.
'''
# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig_twoepoch", Models_2D.sym_mig_twoepoch, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m1, m2, T1, T2")


# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig_twoepoch", Models_2D.asym_mig_twoepoch, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12a, m21a, m12b, m21b, T1, T2")


# Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_sym_mig_three_epoch", Models_2D.sec_contact_sym_mig_three_epoch, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2, T3")


# Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_asym_mig_three_epoch", Models_2D.sec_contact_asym_mig_three_epoch, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2, T3")


# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_sym_mig_size_three_epoch", Models_2D.sec_contact_sym_mig_size_three_epoch, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3")


# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_asym_mig_size_three_epoch", Models_2D.sec_contact_asym_mig_size_three_epoch, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3")




'''
Island Model Set

This set of models come from the following publication:

    Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K.,
    Gvodk, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. Sky, sea,
    and forest islands: diversification in the African leaf-folding frog
    Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). Early Access,
    Journal of Biogeography.

For all the following models, pop2 is assumed to be the 'island' population,
and pop2=nuA(s), pop1=nuA(1-s), where NuA = ancestral population and 's' is
a fraction. Vicariance events involve no population size change, whereas
founder event models always enforce exponential growth in pop2
(the island population).

Note: For these models, we must include an upper bound optional argument in
order to constrain the parameter 's' to be no more than 0.5, and we must also
supply the starting parameters to keep this within these bounds. This is
similar for models with discrete admixture events, for which we need to
constrain 'f' to between 0-1.

In general, you would not want to run these with the above model set,
because they are constructed quite differently.

'''

# Island: Vicariance with no migration.
up = [20, 20, 20, 10, 0.5]
ps = [1,1,1,1,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "vic_no_mig", Models_2D.vic_no_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T, s", in_upper=up, in_params=ps)


# Island: Vicariance with with ancient continuous asymmetric migration.
up = [20, 20, 20, 10, 10, 10, 10, 0.5]
ps = [1,1,1,1,1,1,1,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "vic_anc_asym_mig", Models_2D.vic_anc_asym_mig, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, m12, m21, T1, T2, s", in_upper=up, in_params=ps)


# Island: Vicariance with no migration, secondary contact with continuous asymmetric migration
up = [20, 20, 20, 10, 10, 10, 10, 0.5]
ps = [1,1,1,1,1,1,1,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "vic_sec_contact_asym_mig", Models_2D.vic_sec_contact_asym_mig, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, m12, m21, T1, T2, s", in_upper=up, in_params=ps)


# Island: Founder event with no migration.
up = [20, 20, 20, 10, 0.5]
ps = [1,1,1,1,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "founder_nomig", Models_2D.founder_nomig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T, s", in_upper=up, in_params=ps)


# Island: Founder event with continuous symmetric migration.
up = [20, 20, 20, 20, 10, 0.5]
ps = [1,1,1,1,1,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "founder_sym", Models_2D.founder_sym, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, m, T, s", in_upper=up, in_params=ps)


# Island: Founder event with continuous asymmetric migration.
up = [20, 20, 20, 20, 20, 10, 0.5]
ps = [1,1,1,1,1,1,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "founder_asym", Models_2D.founder_asym, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, m12, m21, T, s", in_upper=up, in_params=ps)


# Island: Vicariance, early unidirectional discrete admixture event (before any drift).
up = [20, 20, 20, 10, 0.5, 0.99]
ps = [1,1,1,1,0.25,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "vic_no_mig_admix_early", Models_2D.vic_no_mig_admix_early, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T, s, f", in_upper=up, in_params=ps)


# Island: Vicariance, late unidirectional discrete admixture event (after any drift).
up = [20, 20, 20, 10, 0.5, 0.99]
ps = [1,1,1,1,0.25,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "vic_no_mig_admix_late", Models_2D.vic_no_mig_admix_late, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T, s, f", in_upper=up, in_params=ps)


# Island: Vicariance, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
up = [20, 20, 20, 10, 10, 0.5, 0.99]
ps = [1,1,1,1,1,0.25,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "vic_two_epoch_admix", Models_2D.vic_two_epoch_admix, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T1, T2, s, f", in_upper=up, in_params=ps)


# Founder event with no migration, early unidirectional discrete admixture event.
up = [20, 20, 10, 10, 0.5, 0.99]
ps = [1,1,1,1,0.25,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "founder_nomig_admix_early", Models_2D.founder_nomig_admix_early, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T, s, f", in_upper=up, in_params=ps)


# Founder event with no migration, late unidirectional discrete admixture event.
up = [20, 20, 10, 10, 0.5, 0.99]
ps = [1,1,1,1,0.25,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "founder_nomig_admix_late", Models_2D.founder_nomig_admix_late, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T, s, f", in_upper=up, in_params=ps)


# Island: Founder event, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
up = [20, 20, 10, 10, 10, 0.5, 0.99]
ps = [1,1,1,1,1,0.25,0.25]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "founder_nomig_admix_two_epoch", Models_2D.founder_nomig_admix_two_epoch, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuA, nu1, nu2, T1, T2, s, f", in_upper=up, in_params=ps)


