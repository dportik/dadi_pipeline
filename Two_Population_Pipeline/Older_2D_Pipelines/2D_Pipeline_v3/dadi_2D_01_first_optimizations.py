import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
usage: python dadi_2D_01_first_optimizations.py

Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 
This is where all the population models and functions are stored for this script. 

Script will perform optimizations from multiple starting points using a
3-fold perturbed set of random starting values for parameters. The output for
each model is a tab-delimited text file which can be opened and sorted to
find the best scoring replicate. The parameter values from this run should
be used as the starting values in the next script, "dadi_2D_02_second_optimizations.py".
The order of output parameters matches the input order, so they can be copied
and pasted into the next script within the appropriate list.

The outputs from here will be labeled according to model and a user
selected prefix name, and are written to the working directory:

"Round1_[prefix]_[model_name]_optimized.txt"

example output file from one model:

Model	param_set	Replicate	log-likelihood	theta	AIC	optimized_params
Divergence with no migration	parameter set = [nu1, nu2, T]	1	-443.41	317.15	892.82	0.6795	0.7315	0.141	
Divergence with no migration	parameter set = [nu1, nu2, T]	2	-179.8	141.7	365.6	4.3798	1.6121	0.7995	
Divergence with no migration	parameter set = [nu1, nu2, T]	3	-241.28	116.99	488.56	3.1119	17.7282	1.1496	
Divergence with no migration	parameter set = [nu1, nu2, T]	4	-1763.04	54.2	3532.08	9.0388	0.2684	5.3911	
Divergence with no migration	parameter set = [nu1, nu2, T]	5	-177.49	150.28	360.98	2.7819	2.6884	0.8431	

*Note the params are written to as many columns as there are values, without any column heading.


Requires user to edit sections of code marked with #**************

You'll absolutely need to provide the path to your SNPs input file
along with your specific projections and population labels. 


############################################
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-Matplotlib
-dadi
############################################

Dan Portik
daniel.portik@uta.edu -> danielportik@email.arizona.edu
October 2017
'''
#keep track of start time
t_begin = datetime.now()


#===========================================================================
#get snps file and convert into allele frequency spectrum object in dadi

#**************
snps1 = "/FULL PATH TO /dadi_2pops_North_South_snps.txt"

#Create python dictionary from snps file
dd1 = dadi.Misc.make_data_dict(snps1)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["North", "South"]
#projection sizes, in ALLELES not individuals
proj_1 = [16,32]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)

print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'


#======================================================================================
# Now prepare to run the optimization round 1 function, which is defined in the 
# script 'Optimize_Functions.py'.

# Function:
# Optimize_Round1(pts, fs, outfile, reps, maxiter, model_name)

# Argument definitions:
# pts:  grid choice (list of three numbers, ex. [20,30,40]
# fs:  spectrum object name
# outfile:  prefix for output naming -> "Round1_{0}_{1}_optimized.txt".format(outfile,model_name)
# reps:  integer to control number of replicates, ex. 10
# maxiter:  max number of iterations per optimization step (not intuitive! see dadi user group)
# model_name:  "no_divergence", "no_mig", "sym_mig", "asym_mig", "anc_sym_mig", "anc_asym_mig",
#        "sec_contact_sym_mig", "sec_contact_asym_mig", "no_mig_size", "sym_mig_size",
#        "asym_mig_size", "anc_sym_mig_size", "anc_asym_mig_size", "sec_contact_sym_mig_size",
#        "sec_contact_asym_mig_size", "sym_mig_twoepoch", "asym_mig_twoepoch", 
#		 "sec_contact_sym_mig_three_epoch", "sec_contact_asym_mig_three_epoch", 
#	     "sec_contact_sym_mig_size_three_epoch", "sec_contact_asym_mig_size_three_epoch", 
#		 "vic_no_mig", "vic_anc_asym_mig", "vic_sec_contact_asym_mig", "founder_nomig", 
#        "founder_sym", "founder_asym", "vic_no_mig_admix_early", "vic_no_mig_admix_late", 
#        "vic_two_epoch_admix", "founder_nomig_admix_early", "founder_nomig_admix_late",
#        "founder_nomig_admix_two_epoch"

#======================================================================================
#Input some of the basic reusable arguments here specific to your data set

#These are specific to your data set:
#**************
#grid choice
pts = [50,60,70]
#prefix for output file naming
outfile = "N_v_S"

#These can be left alone, unless you want more searches:
#spectrum object name (we defined this above)
fs = fs_1
#integer to control number of replicates per model
reps = int(50)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(20)

#======================================================================================
# Now call the function with the relevant arguments.

# There are many models to test here. A brief definition is given for each, but the actual
# models are defined in the Models_2D.py script. The first 15 were implemented in Portik 
# et al. 2016 (doi: 10.1111/mec.14266), the following 9 are newer for various projects.

# Here it is set up to call each model one by one sequentially, which could finish relatively quickly.
# If it takes too long, create multiple verisions of this script, block out some models (use hashes or delete),
# and execute one version for every core you have available. It will greatly speed up these steps,
# and sometimes if extrapolations fail the script will crash too and this could prevent it from
# happening too many times.


# Standard neutral model, populations never diverge
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_divergence")

# Split into two populations, no migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_mig")

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig")

# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig")

# Split with continuous symmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_sym_mig")

# Split with continuous asymmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_asym_mig")

# Split with no gene flow, followed by period of continuous symmetrical gene flow.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig")

# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig")

# Split with no migration, then instantaneous size change with no migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_mig_size")

# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig_size")

# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig_size")

# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.  
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_sym_mig_size")

# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_asym_mig_size")

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size")

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size")


# Newer Models.
# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig_twoepoch")

# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig_twoepoch")

# Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_three_epoch")

# Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_three_epoch")

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size_three_epoch")

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size_three_epoch")


###### 'Island' specific models.

# Island: Vicariance with no migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "vic_no_mig")

# Island: Vicariance with with ancient continuous asymmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "vic_anc_asym_mig")

# Island: Vicariance with no migration, secondary contact with continuous asymmetric migration
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "vic_sec_contact_asym_mig")

# Island: Founder event with no migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "founder_nomig")

# Island: Founder event with continuous symmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "founder_sym")

# Island: Founder event with continuous asymmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "founder_asym")

# Island: Vicariance, early unidirectional discrete admixture event (before any drift).
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "vic_no_mig_admix_early")

# Island: Vicariance, late unidirectional discrete admixture event (after any drift).
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "vic_no_mig_admix_late")

# Island: Vicariance, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "vic_two_epoch_admix")

# Founder event with no migration, early unidirectional discrete admixture event.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_early")

# Founder event with no migration, late unidirectional discrete admixture event.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_late")

# Island: Founder event, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_two_epoch")



#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================

