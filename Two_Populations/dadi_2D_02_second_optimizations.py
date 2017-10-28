import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
usage: python dadi_2D_02_second_optimizations.py

Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 
This is where all the population models and functions are stored for this script. 

Script will perform optimizations from multiple starting points using a
2-fold perturbed set of USER SELECTED starting values for parameters. The output for
each model is a tab-delimited text file which can be opened and sorted to
find the best scoring replicate. The parameter values from this run can
be used as the starting values in the next script, "dadi_2D_03_third_optimizations.py"
if you feel that the likelihoods are still somewhat inconsistent and not finding
an optimal solution. That script will perturb input values of parameters only 1-fold,
so you should see more consistent results if your likelihood surface isn't flat.
The order of output parameters matches the input order, so they can be copied
and pasted into the next script within the appropriate list.

The outputs from here will be labeled according to model and a user
selected prefix name, and are written to the working directory:

"Round2_[prefix]_[model_name]_optimized.txt"

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

print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'


#======================================================================================
# Now prepare to run the optimization round 2 function, which is defined in the 
# script 'Optimize_Functions.py'.

# Function:
# Optimize_Round2(pts, fs, outfile, reps, maxiter, model_name, params)

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

#"no_mig"
#3 Values
no_mig_params = [0.0622,0.0716,0.0683]

#"sym_mig"
#4 Values
sym_mig_params = [0.3912,0.3509,0.1655,0.6334]

#"asym_mig"
#5 Values
asym_mig_params = [0.0881,0.0946,0.8871,0.438,0.1703]

#"anc_sym_mig"
#5 Values
anc_sym_mig_params = [0.7087,0.656,0.3311,5.3682,0.1292]

#"anc_asym_mig"
#6 Values
anc_asym_mig_params = [0.384,0.3265,0.2192,0.3929,0.2113,0.2249]

#"sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [0.8282,0.6621,0.1142,0.4408,2.7886]

#"sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.1871,0.263,1.1289,0.1161,0.4958,0.2885]

#"no_mig_size"
#6 Values
no_mig_size_params = [3.6872,0.9558,0.6311,0.6907,0.0935,0.8674]

#"sym_mig_size"
#7 Values
sym_mig_size_params = [11.0922,0.3099,0.6772,0.5948,0.1586,0.0934,2.8574]

#"asym_mig_size"
#8 Values
asym_mig_size_params = [1.7551,0.1893,0.2261,0.4507,0.3175,0.2954,0.2775,0.3068]

#"anc_sym_mig_size"
#7 Values
anc_sym_mig_size_params = [4.8808,3.2591,0.6068,0.4988,0.5057,5.9915,0.4379]

#"anc_asym_mig_size"
#8 Values
anc_asym_mig_size_params = [2.3122,0.4563,0.2102,0.4988,0.5503,0.1259,0.8031,0.1847]

#"sec_contact_sym_mig_size"
#7 Values
sec_contact_sym_mig_size_params = [6.1645,0.6576,0.2344,0.3341,0.2551,0.5079,0.2505]

#"sec_contact_asym_mig_size"
#8 Values
sec_contact_asym_mig_size_params = [0.6235,0.4559,0.4373,0.3229,0.1855,0.289,1.062,0.2822]

#"sym_mig_twoepoch" 
# 6 Values
sym_mig_twoepoch_params = [0.2425,0.2489,0.0964,0.3329,0.2021,0.3038]

#"asym_mig_twoepoch" 
# 8 Values
asym_mig_twoepoch_params = [0.1191,0.1429,0.4132,8.8809,2.3879,0.1926,0.2715,2.2659]

#"sec_contact_sym_mig_three_epoch" 
# 6 Values
sec_contact_sym_mig_three_epoch_params = [0.4486,0.35,5.8658,5.4086,0.3468,0.3818]

#"sec_contact_asym_mig_three_epoch" 
# 7 Values
sec_contact_asym_mig_three_epoch_params = [1.1834,1.1633,0.5618,0.5342,9.497,0.8519,1.021]

#"sec_contact_sym_mig_size_three_epoch" 
# 8 Values
sec_contact_sym_mig_size_three_epoch_params = [0.1463,3.2944,1.9881,2.253,0.1718,7.955,9.6366,0.7891]

#"sec_contact_asym_mig_size_three_epoch" 
# 9 Values
sec_contact_asym_mig_size_three_epoch_params = [0.1592,0.6018,1.5768,1.3381,0.4383,0.2937,2.3318,0.4379,0.5023]

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
reps = int(50)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(30)


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
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "no_divergence", no_divergence_params)

# Split into two populations, no migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)

# Split with continuous symmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_sym_mig", anc_sym_mig_params)

# Split with continuous asymmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_asym_mig", anc_asym_mig_params)

# Split with no gene flow, followed by period of continuous symmetrical gene flow.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig", sec_contact_sym_mig_params)

# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig", sec_contact_asym_mig_params)

# Split with no migration, then instantaneous size change with no migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "no_mig_size", no_mig_size_params)

# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sym_mig_size", sym_mig_size_params)

# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "asym_mig_size", asym_mig_size_params)

# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.  
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_sym_mig_size", anc_sym_mig_size_params)

# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_asym_mig_size", anc_asym_mig_size_params)

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size", sec_contact_sym_mig_size_params)

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size", sec_contact_asym_mig_size_params)


# Newer Models
# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sym_mig_twoepoch", sym_mig_twoepoch_params)

# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "asym_mig_twoepoch", asym_mig_twoepoch_params)

# Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_three_epoch", sec_contact_sym_mig_three_epoch_params)

# Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_three_epoch", sec_contact_asym_mig_three_epoch_params)

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size_three_epoch", sec_contact_sym_mig_size_three_epoch_params)

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size_three_epoch", sec_contact_asym_mig_size_three_epoch_params)


###### 'Island' specific models.

# Island: Vicariance with no migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "vic_no_mig", vic_no_mig_params)

# Island: Vicariance with with ancient continuous asymmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "vic_anc_asym_mig", vic_anc_asym_mig_params)

# Island: Vicariance with no migration, secondary contact with continuous asymmetric migration
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "vic_sec_contact_asym_mig", vic_sec_contact_asym_mig_params)

# Island: Founder event with no migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "founder_nomig", founder_nomig_params)

# Island: Founder event with continuous symmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "founder_sym", founder_sym_params)

# Island: Founder event with continuous asymmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "founder_asym", founder_asym_params)

# Island: Vicariance, early unidirectional discrete admixture event (before any drift).
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "vic_no_mig_admix_early", vic_no_mig_admix_early_params)

# Island: Vicariance, late unidirectional discrete admixture event (after any drift).
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "vic_no_mig_admix_late", vic_no_mig_admix_late_params)

# Island: Vicariance, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "vic_two_epoch_admix", vic_two_epoch_admix_params)

# Founder event with no migration, early unidirectional discrete admixture event.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_early", founder_nomig_admix_early_params)

# Founder event with no migration, late unidirectional discrete admixture event.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_late", founder_nomig_admix_late_params)

# Island: Founder event, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "founder_nomig_admix_two_epoch", founder_nomig_admix_two_epoch_params)


#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


