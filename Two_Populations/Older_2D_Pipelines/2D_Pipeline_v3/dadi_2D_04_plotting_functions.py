import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime
import matplotlib.pyplot as plt

'''
usage: python dadi_2D_04_plotting_functions.py

Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 
This is where all the population models and functions are stored for this script. 

Script will simulate the model using a fixed set of input parameters. That is,
no optimizations will occur here, and the assumption is you've searched
thoroughly for the parameter set with the highest likelihood value. The goal
of this script is to produce plots of the 2D-JSFS for the data and model,
plus the residuals. A pdf file will be written for each model to the
current working directory.

Requires user to edit sections of code marked with #**************

You'll absolutely need to provide the path to your SNPs input file
along with your specific projections and population labels. 

###########################################################
**Note, if you see this error when plotting:
"ValueError: Data has no positive values, and therefore can not be log-scaled."
You will need to change the vmin in the plotting routine from None to something 0<1.

I have changed vmin to 0.005-0.01 with good results. So, at the bottom of this script:

dadi.Plotting.plot_2d_comp_multinom(model, data, resid_range = 3)

becomes:

dadi.Plotting.plot_2d_comp_multinom(model, data, resid_range = 3, vmin=0.005)

This should eliminate that particular error.
###########################################################

###########################################################
*Note: This may be version specific to both dadi and numpy, but
I kept getting a plotting error and so changed some lines in the
dadi Plotting.py module starting with line 395:

    ax4 = pylab.subplot(2,2,4)
    flatresid = numpy.compress(numpy.logical_not(resid.mask.ravel()), 
                               resid.ravel())
    hister, bin_edges = numpy.histogram(flatresid, bins = 30, range=None)
    #ax4.hist(flatresid, bins=20, normed=True)
    ax4.bar(bin_edges[:-1], hister, width=0.2)
    ax4.set_title('residuals')
    ax4.set_yticks([])
    
This solved an issue with numpy incompatibility (instead of using the 
ax4.hist syntax). An updated numpy should however work fine as is.
###########################################################

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

#print relevant info to screen
print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'


#======================================================================================
# Now prepare to optimize the models, which is defined in the script 'Optimize_Functions.py'.

# Function usage:
# returned_model_name = Optimize_Single(pts, fs, model_name, params)

# Argument definitions:
# pts:  grid choice (list of three numbers, ex. [20,30,40]
# fs:  spectrum object name
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
# params:  list of best parameter values to optimize with

#===========================================================================
# Here we need to enter the parameter values for the previous best scoring replicate. Unlike
# before, these are the EXACT values that will be used to optimize. You should get these
# from the output of the previous script. The default order of the output file parameters 
# will be the same order they go in here, so they can just be copy pasted over.

#"no_divergence"
#leave blank, no parameters
no_divergence_params = []

#**************

#"no_mig"
#3 Values
no_mig_params = [0.1039,0.0992,0.109]

#"sym_mig"
#4 Values
sym_mig_params = [0.1502,0.1371,0.2473,0.188]

#"asym_mig"
#5 Values
asym_mig_params = [0.1493,0.1351,0.2418,0.2552,0.1864]

#"anc_sym_mig"
#5 Values
anc_sym_mig_params = [0.3381,0.2909,0.2607,0.6751,0.0001]

#"anc_asym_mig"
#6 Values
anc_asym_mig_params = [0.1035,0.0966,0.1321,0.8946,0.0309,0.0757]

#"sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [0.359,0.3356,0.1333,0.2051,0.3733]

#"sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.1713,0.1529,0.4419,0.6281,0.1888,0.0127]

#"no_mig_size"
#6 Values
no_mig_size_params = [3.9744,1.8151,0.0255,0.0333,0.8283,0.0201]

#"sym_mig_size"
#7 Values
sym_mig_size_params = [20.181,0.3609,0.2571,0.2502,0.1963,0.0118,0.3336]

#"asym_mig_size"
#8 Values
asym_mig_size_params = [2.1611,0.4093,0.0228,0.0856,0.9088,0.2608,0.8099,0.0212]

#"anc_sym_mig_size"
#7 Values
anc_sym_mig_size_params = [1.0925,1.029,0.0701,0.0817,0.11,2.661,0.0122]

#"anc_asym_mig_size"
#8 Values
anc_asym_mig_size_params = [1.7858,0.9298,0.0698,0.1742,0.1013,0.1285,3.9189,0.0213]

#"sec_contact_sym_mig_size"
#7 Values
sec_contact_sym_mig_size_params = [14.2991,1.281,0.0374,0.1584,0.9115,2.249,0.0334]

#"sec_contact_asym_mig_size"
#8 Values
sec_contact_asym_mig_size_params = [0.9971,0.4688,0.0202,0.0471,0.2187,1.4038,0.5891,0.0101]

#"sym_mig_twoepoch" 
# 6 Values
sym_mig_twoepoch_params = [0.1514,0.1393,0.066,0.269,0.1339,0.0509]

#"asym_mig_twoepoch" 
# 8 Values
asym_mig_twoepoch_params = [0.0879,0.1266,2.2299,9.0892,1.9294,0.4392,0.3187,2.7066]

#"sec_contact_sym_mig_three_epoch" 
# 6 Values
sec_contact_sym_mig_three_epoch_params = [0.2193,0.204,17.9698,7.5162,0.586,0.2243]

#"sec_contact_asym_mig_three_epoch" 
# 7 Values
sec_contact_asym_mig_three_epoch_params = [0.5657,0.3789,1.0008,6.5421,9.9127,0.4171,0.5922]

#"sec_contact_sym_mig_size_three_epoch" 
# 8 Values
sec_contact_sym_mig_size_three_epoch_params = [0.2275,11.3693,0.8588,0.8269,0.5425,9.089,1.9292,0.5098]

#"sec_contact_asym_mig_size_three_epoch" 
# 9 Values
sec_contact_asym_mig_size_three_epoch_params = [0.0322,1.06,0.4402,0.5182,1.0017,0.4149,5.7379,1.0344,0.2081]

#"founder_sym" 
# 6 Values
founder_sym_params = [1.5933,1.4055,0.1415,0.1483,0.4247,0.1949]

#"founder_asym" 
# 7 Values
founder_asym_params = [1.5758,2.6038,0.205,0.2628,0.1122,0.5803,0.2188]

#"founder_nomig" 
# 5 Values
founder_nomig_params = [1.9544,2.7918,0.1119,0.3278,0.151]


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

#===========================================================================
# Now call the function with the relevant arguments.

# There are several models to optimize.
# Each model is executed with one replicate using fixed parameter values.
# The simulated model is stored as an object to be called on in the plotting function.

# Standard neutral model, populations never diverge
no_divergence = Optimize_Functions.Optimize_Single(pts, fs, "no_divergence", no_divergence_params)

# Split into two populations, no migration.
no_mig = Optimize_Functions.Optimize_Single(pts, fs, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
sym_mig = Optimize_Functions.Optimize_Single(pts, fs, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
asym_mig = Optimize_Functions.Optimize_Single(pts, fs, "asym_mig", asym_mig_params)

# Split with continuous symmetric migration, followed by isolation.
anc_sym_mig = Optimize_Functions.Optimize_Single(pts, fs, "anc_sym_mig", anc_sym_mig_params)

# Split with continuous asymmetric migration, followed by isolation.
anc_asym_mig = Optimize_Functions.Optimize_Single(pts, fs, "anc_asym_mig", anc_asym_mig_params)

# Split with no gene flow, followed by period of continuous symmetrical gene flow.
sec_contact_sym_mig = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_sym_mig", sec_contact_sym_mig_params)

# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
sec_contact_asym_mig = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_asym_mig", sec_contact_asym_mig_params)

# Split with no migration, then instantaneous size change with no migration.
no_mig_size = Optimize_Functions.Optimize_Single(pts, fs, "no_mig_size", no_mig_size_params)

# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
sym_mig_size = Optimize_Functions.Optimize_Single(pts, fs, "sym_mig_size", sym_mig_size_params)

# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
asym_mig_size = Optimize_Functions.Optimize_Single(pts, fs, "asym_mig_size", asym_mig_size_params)

# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.  
anc_sym_mig_size = Optimize_Functions.Optimize_Single(pts, fs, "anc_sym_mig_size", anc_sym_mig_size_params)

# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
anc_asym_mig_size = Optimize_Functions.Optimize_Single(pts, fs, "anc_asym_mig_size", anc_asym_mig_size_params)

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
sec_contact_sym_mig_size = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_sym_mig_size", sec_contact_sym_mig_size_params)

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
sec_contact_asym_mig_size = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_asym_mig_size", sec_contact_asym_mig_size_params)


# Newer Models
# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
sym_mig_twoepoch = Optimize_Functions.Optimize_Single(pts, fs, "sym_mig_twoepoch", sym_mig_twoepoch_params)

# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
asym_mig_twoepoch = Optimize_Functions.Optimize_Single(pts, fs, "asym_mig_twoepoch", asym_mig_twoepoch_params)

# Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
sec_contact_sym_mig_three_epoch = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_sym_mig_three_epoch", sec_contact_sym_mig_three_epoch_params)

# Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
sec_contact_asym_mig_three_epoch = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_asym_mig_three_epoch", sec_contact_asym_mig_three_epoch_params)

# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
sec_contact_sym_mig_size_three_epoch = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_sym_mig_size_three_epoch", sec_contact_sym_mig_size_three_epoch_params)

# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.
sec_contact_asym_mig_size_three_epoch = Optimize_Functions.Optimize_Single(pts, fs, "sec_contact_asym_mig_size_three_epoch", sec_contact_asym_mig_size_three_epoch_params)



###### 'Island' specific models.

# Island: Vicariance with no migration.
vic_no_mig = Optimize_Functions.Optimize_Single(pts, fs, "vic_no_mig", vic_no_mig_params)

# Island: Vicariance with with ancient continuous asymmetric migration.
vic_anc_asym_mig = Optimize_Functions.Optimize_Single(pts, fs, "vic_anc_asym_mig", vic_anc_asym_mig_params)

# Island: Vicariance with no migration, secondary contact with continuous asymmetric migration
vic_sec_contact_asym_mig = Optimize_Functions.Optimize_Single(pts, fs, "vic_sec_contact_asym_mig", vic_sec_contact_asym_mig_params)

# Island: Founder event with no migration.
founder_nomig = Optimize_Functions.Optimize_Single(pts, fs, "founder_nomig", founder_nomig_params)

# Island: Founder event with continuous symmetric migration.
founder_sym = Optimize_Functions.Optimize_Single(pts, fs, "founder_sym", founder_sym_params)

# Island: Founder event with continuous asymmetric migration.
founder_asym = Optimize_Functions.Optimize_Single(pts, fs, "founder_asym", founder_asym_params)

# Island: Vicariance, early unidirectional discrete admixture event (before any drift).
vic_no_mig_admix_early = Optimize_Functions.Optimize_Single(pts, fs, "vic_no_mig_admix_early", vic_no_mig_admix_early_params)

# Island: Vicariance, late unidirectional discrete admixture event (after any drift).
vic_no_mig_admix_late = Optimize_Functions.Optimize_Single(pts, fs, "vic_no_mig_admix_late", vic_no_mig_admix_late_params)

# Island: Vicariance, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
vic_two_epoch_admix = Optimize_Functions.Optimize_Single(pts, fs, "vic_two_epoch_admix", vic_two_epoch_admix_params)

# Founder event with no migration, early unidirectional discrete admixture event.
founder_nomig_admix_early = Optimize_Functions.Optimize_Single(pts, fs, "founder_nomig_admix_early", founder_nomig_admix_early_params)

# Founder event with no migration, late unidirectional discrete admixture event.
founder_nomig_admix_late = Optimize_Functions.Optimize_Single(pts, fs, "founder_nomig_admix_late", founder_nomig_admix_late_params)

# Island: Founder event, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
founder_nomig_admix_two_epoch = Optimize_Functions.Optimize_Single(pts, fs, "founder_nomig_admix_two_epoch", founder_nomig_admix_two_epoch_params)


#======================================================================================
# Define an editable plotting function for data and model comparison

# See instructions at top for troubleshooting issues with plotting, here is where you
# would edit the arguments for dadi.Plotting.plot_2d_comp_multinom()

def plot(sim_model, data, outfile, model_name):
    print '{0}_{1}.pdf'.format(outfile,model_name), '\n'
    outname = '{0}_{1}.pdf'.format(outfile,model_name)
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(sim_model, data, resid_range = 3)
    fig.savefig(outname)
    
#======================================================================================
# Call the plotting function for all optimized models

# Function
# plot_all(sim_model, fs, outfile, model_name)

# Argument defintions:
# sim_model:  the simulated model returned from function Two_Pop_Plot
# fs:  spectrum object name (ie the "data")
# outfile:  same prefix for output naming of pdf files -> "{0}_{1}.pdf".format(outfile,model_name)
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

# A plot should pop up with the model for each, close the box and the next will load. They
# are all saved automatically to the working directory.

plot(no_divergence, fs, outfile, "no_divergence")
plot(no_mig, fs, outfile, "no_mig")
plot(sym_mig, fs, outfile, "sym_mig")
plot(asym_mig, fs, outfile, "asym_mig")
plot(anc_sym_mig, fs, outfile, "anc_sym_mig")
plot(anc_asym_mig, fs, outfile, "anc_asym_mig")
plot(sec_contact_sym_mig, fs, outfile, "sec_contact_sym_mig")
plot(sec_contact_asym_mig, fs, outfile, "sec_contact_asym_mig")
plot(no_mig_size, fs, outfile, "no_mig_size")
plot(sym_mig_size, fs, outfile, "sym_mig_size")
plot(asym_mig_size, fs, outfile, "asym_mig_size")
plot(anc_sym_mig_size, fs, outfile, "anc_sym_mig_size")
plot(anc_asym_mig_size, fs, outfile, "anc_asym_mig_size")
plot(sec_contact_sym_mig_size, fs, outfile, "sec_contact_sym_mig_size")
plot(sec_contact_asym_mig_size, fs, outfile, "sec_contact_asym_mig_size")
plot(sym_mig_twoepoch, fs, outfile, "sym_mig_twoepoch")
plot(asym_mig_twoepoch, fs, outfile, "asym_mig_twoepoch")
plot(sec_contact_sym_mig_three_epoch, fs, outfile, "sec_contact_sym_mig_three_epoch")
plot(sec_contact_asym_mig_three_epoch, fs, outfile, "sec_contact_asym_mig_three_epoch")
plot(sec_contact_sym_mig_size_three_epoch, fs, outfile, "sec_contact_sym_mig_size_three_epoch")
plot(sec_contact_asym_mig_size_three_epoch, fs, outfile, "sec_contact_asym_mig_size_three_epoch")
plot(vic_no_mig, fs, outfile, "vic_no_mig")
plot(vic_anc_asym_mig, fs, outfile, "vic_anc_asym_mig")
plot(vic_sec_contact_asym_mig, fs, outfile, "vic_sec_contact_asym_mig")
plot(founder_nomig, fs, outfile, "founder_nomig")
plot(founder_sym, fs, outfile, "founder_sym")
plot(founder_asym, fs, outfile, "founder_asym")
plot(vic_no_mig_admix_early, fs, outfile, "vic_no_mig_admix_early")
plot(vic_no_mig_admix_late, fs, outfile, "vic_no_mig_admix_late")
plot(vic_two_epoch_admix, fs, outfile, "vic_two_epoch_admix")
plot(founder_nomig_admix_early, fs, outfile, "founder_nomig_admix_early")
plot(founder_nomig_admix_late, fs, outfile, "founder_nomig_admix_late")
plot(founder_nomig_admix_two_epoch, fs, outfile, "founder_nomig_admix_two_epoch")

#======================================================================================
