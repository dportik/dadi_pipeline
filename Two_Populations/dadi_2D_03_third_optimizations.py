import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
usage: python dadi_2D_03_third_optimizations.py

Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 
This is where all the population models and functions are stored for this script. 

Script will perform optimizations from multiple starting points using a
1-fold perturbed set of USER SELECTED starting values for parameters. The output for
each model is a tab-delimited text file which can be opened and sorted to
find the best scoring replicate based on ll or AIC.

The outputs from here will be labeled according to model and a user
selected prefix name, and are written to the working directory:

"Round3_[prefix]_[model_name]_optimized.txt"

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
along with your specific projections and population labels. Look through
the entire script for other sections requiring your attention.


############################################
Written for Python 2.7
Python dependencies required:
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

#print relevant info to screen
print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================
# Now prepare to run the optimization round 3 function, which is defined in the 
# script 'Optimize_Functions.py'.

# Function:
# Optimize_Round3(pts, fs, outfile, reps, maxiter, model_name, params)

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
# params:  list of best parameter values to perturb to start the optimizations from


#===========================================================================
# Here we need to enter the parameter values for the previous best scoring replicate, so we 
# can start with those values, perturb them, and run optimizations again. You should get these
# from the output of the previous script. The default order of the output file parameters 
# will be the same order they go in here, so they can just be copy pasted over.

#"no_divergence"
#leave blank, no parameters
no_divergence_params = []

#**************

#"no_mig"
#3 Values
no_mig_params = [0.1047,0.0995,0.1094]

#"sym_mig"
#4 Values
sym_mig_params = [0.1426,0.1329,0.2587,0.1807]

#"asym_mig"
#5 Values
asym_mig_params = [0.1309,0.1281,0.3641,0.2054,0.1666]

#"anc_sym_mig"
#5 Values
anc_sym_mig_params = [0.5781,0.4977,0.1879,1.9393,0.0248]

#"anc_asym_mig"
#6 Values
anc_asym_mig_params = [0.1804,0.1557,0.0879,0.999,0.1226,0.0839]

#"sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [0.5862,0.5004,0.1536,0.3196,1.8769]

#"sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.2655,0.2406,0.2704,0.2566,0.225,0.1914]

#"no_mig_size"
#6 Values
no_mig_size_params = [4.4381,0.3583,0.2522,0.246,0.0327,0.2786]

#"sym_mig_size"
#7 Values
sym_mig_size_params = [21.5615,0.5914,0.5348,0.4467,0.175,0.1188,1.6146]

#"asym_mig_size"
#8 Values
asym_mig_size_params = [3.6827,0.387,0.0896,0.1855,0.8176,0.2717,0.6448,0.0858]

#"anc_sym_mig_size"
#7 Values
anc_sym_mig_size_params = [1.5776,1.5117,0.2912,0.3076,0.1617,2.2767,0.1215]

#"anc_asym_mig_size"
#8 Values
anc_asym_mig_size_params = [0.8985,1.1107,0.2179,0.2933,0.5017,0.1264,2.6391,0.1036]

#"sec_contact_sym_mig_size"
#7 Values
sec_contact_sym_mig_size_params = [12.5722,1.712,0.1455,0.2442,0.1197,1.2068,0.156]

#"sec_contact_asym_mig_size"
#8 Values
sec_contact_asym_mig_size_params = [0.5375,0.339,0.1337,0.1501,0.0543,0.9015,0.4737,0.0389]

#"sym_mig_twoepoch" 
# 6 Values
sym_mig_twoepoch_params = [0.1526,0.1324,0.0651,0.2601,0.0603,0.1262]

#"asym_mig_twoepoch" 
# 8 Values
asym_mig_twoepoch_params = [0.0469,0.0795,1.5763,4.0853,3.6007,0.5837,0.3399,4.2534]

#"sec_contact_sym_mig_three_epoch" 
# 6 Values
sec_contact_sym_mig_three_epoch_params = [0.2553,0.2395,14.824,6.6129,0.5088,0.2637]

#"sec_contact_asym_mig_three_epoch" 
# 7 Values
sec_contact_asym_mig_three_epoch_params = [0.9464,0.5687,0.2956,2.1062,8.6712,0.5747,0.5933]

#"sec_contact_sym_mig_size_three_epoch" 
# 8 Values
sec_contact_sym_mig_size_three_epoch_params = [0.1305,5.8714,1.0551,1.0063,0.2353,7.094,3.4126,0.387]

#"sec_contact_asym_mig_size_three_epoch" 
# 9 Values
sec_contact_asym_mig_size_three_epoch_params = [0.0375,1.5796,0.5774,0.5646,0.7689,0.5669,3.0461,0.8893,0.3395]

#"vic_no_mig"
# 5 Values
vic_no_mig_params = []

#"vic_anc_asym_mig"
# 8 Values
vic_anc_asym_mig_params = []

#"vic_sec_contact_asym_mig"
# 8 Values
vic_sec_contact_asym_mig_params = []

#"founder_nomig"
# 5 Values
founder_nomig_params = []

#"founder_sym"
# 6 Values
founder_sym_params = []

#"founder_asym"
# 7 Values
founder_asym_params = []

#"vic_no_mig_admix_early"
# 6 Values
vic_no_mig_admix_early_params = []

#"vic_no_mig_admix_late"
# 6 Values
vic_no_mig_admix_late_params = []

#"vic_two_epoch_admix"
# 7 Values
vic_two_epoch_admix_params = []

#"founder_nomig_admix_early"
# 6 Values
founder_nomig_admix_early_params = []

#"founder_nomig_admix_late"
# 6 Values
founder_nomig_admix_late_params = []

#"founder_nomig_admix_two_epoch"
# 7 Values
founder_nomig_admix_two_epoch_params = []


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
reps = int(100)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(50)


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
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "no_divergence", no_divergence_params)

# Split into two populations, no migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)

# Split with continuous symmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_sym_mig", anc_sym_mig_params)

# Split with continuous asymmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_asym_mig", anc_asym_mig_params)

# Split with no gene flow, followed by period of continuous symmetrical gene flow.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig", sec_contact_sym_mig_params)

# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig", sec_contact_asym_mig_params)

# Split with no migration, then instantaneous size change with no migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "no_mig_size", no_mig_size_params)

# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sym_mig_size", sym_mig_size_params)

# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "asym_mig_size", asym_mig_size_params)

# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.  
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_sym_mig_size", anc_sym_mig_size_params)

# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_asym_mig_size", anc_asym_mig_size_params)

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size", sec_contact_sym_mig_size_params)

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size", sec_contact_asym_mig_size_params)


# Newer Models
# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sym_mig_twoepoch", sym_mig_twoepoch_params)

# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "asym_mig_twoepoch", asym_mig_twoepoch_params)

# Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_three_epoch", sec_contact_sym_mig_three_epoch_params)

# Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_three_epoch", sec_contact_asym_mig_three_epoch_params)

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size_three_epoch", sec_contact_sym_mig_size_three_epoch_params)

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size_three_epoch", sec_contact_asym_mig_size_three_epoch_params)


###### 'Island' specific models.

# Island: Vicariance with no migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "vic_no_mig", vic_no_mig_params)

# Island: Vicariance with with ancient continuous asymmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "vic_anc_asym_mig", vic_anc_asym_mig_params)

# Island: Vicariance with no migration, secondary contact with continuous asymmetric migration
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "vic_sec_contact_asym_mig", vic_sec_contact_asym_mig_params)

# Island: Founder event with no migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "founder_nomig", founder_nomig_params)

# Island: Founder event with continuous symmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "founder_sym", founder_sym_params)

# Island: Founder event with continuous asymmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "founder_asym", founder_asym_params)

# Island: Vicariance, early unidirectional discrete admixture event (before any drift).
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "vic_no_mig_admix_early", vic_no_mig_admix_early_params)

# Island: Vicariance, late unidirectional discrete admixture event (after any drift).
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "vic_no_mig_admix_late", vic_no_mig_admix_late_params)

# Island: Vicariance, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "vic_two_epoch_admix", vic_two_epoch_admix_params)

# Founder event with no migration, early unidirectional discrete admixture event.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_early", founder_nomig_admix_early_params)

# Founder event with no migration, late unidirectional discrete admixture event.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_late", founder_nomig_admix_late_params)

# Island: Founder event, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_two_epoch", founder_nomig_admix_two_epoch_params)

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


