import sys
import os
import numpy as np
import dadi
import pylab
import TwoD_models
import matplotlib.pyplot as plt
'''
usage: python dadi_2D_bestruns_plotting.py

Requires the TwoD_models.py script to be in same working directory. This is where
all the population model functions are stored.

Overview: input optimized parameters, return best model, plot, save.

Assumes you've optimized your parameters elsewhere. That is, no optimizations
are performed here and you need to use the exact values from you best scoring models.

Requires user to edit sections of code marked with #**************
You'll need to be familiar with dadi to do this.

###########################################################
*Note: I changed some lines in the dadi Plotting.py module
starting with line 395:

    ax4 = pylab.subplot(2,2,4)
    flatresid = numpy.compress(numpy.logical_not(resid.mask.ravel()), 
                               resid.ravel())
    hister, bin_edges = numpy.histogram(flatresid, bins = 30, range=None)
    #ax4.hist(flatresid, bins=20, normed=True)
    ax4.bar(bin_edges[:-1], hister, width=0.2)
    ax4.set_title('residuals')
    ax4.set_yticks([])
    
This solved an issue with numpy incompatibility (instead of using the 
ax4.hist syntax).
You'll need to find this in your python framework and edit the correct file
if you are having trouble with the standard code. Upgrading numpy may
resolve the issue.
###########################################################

**Note, if you see this error when plotting:
"ValueError: Data has no positive values, and therefore can not be log-scaled."
You will need to change the vmin in the plotting routine from None to something 0<1.
You'll need to track the installation path for dadi and change the code in the
Plotting.py script.

I used the following change with good results:
def plot_2d_comp_multinom(model, data, vmin=0.01, vmax=None,
                          resid_range=None, fig_num=None,
                          pop_ids=None, residual='Anscombe',
                          adjust=True):
                          
If you can't locate the Plotting.py script duplicate its contents in another
local script (put in same directory as this), call it Plotting2.py, and
edit that script with the above changes.

Then, at the top of this script also add 'import Plotting2' then at line 673
(in the plotter function) change: 
dadi.Plotting.plot_2d_comp_multinom(model, data, resid_range = 3)
to 
Plotting2.plot_2d_comp_multinom(model, data, resid_range = 3)
This provides a workaround so as to avoid editing the actual dadi scripts in your system.

############################################
Written for Python 2.7.3
External Dependencies:
-Numpy (Numerical Python)
-Scipy
-Matplotlib
-dadi
############################################

Dan Portik
daniel.portik@uta.edu
September 2016
'''

#===========================================================================
#get snps file 
#change to full path of your file
#**************
snps1 = "/Users/dan/Dropbox/dadi_inputs/Scotobleps_Final/2D/CVL_South/dadi_2pops_CVL_South_snps.txt"

#Create python dictionary from snps file
dd1 = dadi.Misc.make_data_dict(snps1)

#Convert this dictionary into folded AFS object

#**************
proj_1 = [40,18]
#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
#[polarized = False] creates folded spectrum object
fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=['CVL','South'], projections = proj_1, polarized = False)

print "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================
#create function to run models with user input parameters (presumably from best runs)
#return model for plotting the 2d spectrum in data vs. model comparisons

def Two_Pop_Models(pts, fs, params, modelchoice):
    if modelchoice == "snm":
        #####################################
        #Neutral model, no split
        #####################################
        print "---------------------------------------------------"
        print "Neutral Model, No Divergence",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_snm = TwoD_models.snm

        #create an extrapolating function 
        func_snm_exec = dadi.Numerics.make_extrap_log_func(func_snm)
        #run the model with the extrapolation call and all our parameters
        snm_model = func_snm_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_neutral = dadi.Inference.ll_multinom(snm_model, fs)
        ll_neutral = np.around(ll_neutral, 3)
        print "likelihood of Neutral Model, No Divergence model = ", ll_neutral

        #calculate theta
        theta_snm = dadi.Inference.optimal_sfs_scaling(snm_model, fs)
        theta_snm = np.around(theta_snm, 3)
        print "Theta = ", theta_snm

        #calculate AIC 
        #AIC = -2 ( ln ( likelihood )) + 2 K
        aic_neutral = ( -2*( float(ll_neutral)))
        print "AIC = ", aic_neutral, '\n', '\n'        
        print "---------------------------------------------------", '\n'
        return snm_model
        

    elif modelchoice == "bottlegrowth":
        #####################################
        #bottlegrowth model, no split
        #####################################
        print "---------------------------------------------------"
        print "Bottlegrowth Model, No Divergence",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_bottle = TwoD_models.bottlegrowth

        #create an extrapolating function 
        func_bottle_exec = dadi.Numerics.make_extrap_log_func(func_bottle)

        #simulate the model with the optimized parameters
        bottle_model = func_bottle_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_bottle = dadi.Inference.ll_multinom(bottle_model, fs)
        ll_bottle = np.around(ll_bottle, 3)
        print "likelihood of Neutral Model, No Divergence model = ", ll_bottle

        #calculate theta
        theta_bottle = dadi.Inference.optimal_sfs_scaling(bottle_model, fs)
        theta_bottle = np.around(theta_bottle, 3)
        print "Theta = ", theta_bottle

        #calculate AIC 
        #AIC = -2 ( ln ( likelihood )) + 2 K
        aic_bottle = ( -2*( float(ll_bottle))) + (2*3)
        print "AIC = ", aic_bottle, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return bottle_model

    elif modelchoice == "split_nomig":
        #####################################
        #Split with no migration
        #####################################
        print "---------------------------------------------------"
        print "Split with No Migration",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_nomig = TwoD_models.split_no_mig

        #create an extrapolating function 
        func_split_nomig_exec = dadi.Numerics.make_extrap_log_func(func_split_nomig)

        #simulate the model with the optimized parameters
        split_nomig_model = func_split_nomig_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_nomig = dadi.Inference.ll_multinom(split_nomig_model, fs)
        ll_split_nomig = np.around(ll_split_nomig, 3)
        print "likelihood of Split with No Migration model = ", ll_split_nomig

        #calculate theta
        theta_split_nomig = dadi.Inference.optimal_sfs_scaling(split_nomig_model, fs)
        theta_split_nomig = np.around(theta_split_nomig, 3)
        print "Theta = ", theta_split_nomig

        #calculate AIC 
        aic_split_nomig = ( -2*( float(ll_split_nomig))) + (2*3)
        print "AIC = ", aic_split_nomig, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_nomig_model
    
    elif modelchoice == "split_symmig":
        #####################################
        #Split with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Symmetric Migration",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_symmig = TwoD_models.split_mig

        #create an extrapolating function 
        func_split_symmig_exec = dadi.Numerics.make_extrap_log_func(func_split_symmig)

        #simulate the model with the optimized parameters
        split_symmig_model = func_split_symmig_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_symmig = dadi.Inference.ll_multinom(split_symmig_model, fs)
        ll_split_symmig = np.around(ll_split_symmig, 3)
        print "likelihood of Split with Symmetric Migration model = ", ll_split_symmig

        #calculate theta
        theta_split_symmig = dadi.Inference.optimal_sfs_scaling(split_symmig_model, fs)
        theta_split_symmig = np.around(theta_split_symmig, 3)
        print "Theta = ", theta_split_symmig

        #calculate AIC 
        aic_split_symmig = ( -2*( float(ll_split_symmig))) + (2*4)
        print "AIC = ", aic_split_symmig, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_symmig_model

    elif modelchoice == "split_asymmig":
        #####################################
        #Split with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Asymmetric Migration",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_asymmig = TwoD_models.split_asym_mig

        #create an extrapolating function 
        func_split_asymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_asymmig)

        #simulate the model with the optimized parameters
        split_asymmig_model = func_split_asymmig_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_asymmig = dadi.Inference.ll_multinom(split_asymmig_model, fs)
        ll_split_asymmig = np.around(ll_split_asymmig, 3)
        print "likelihood of Split with Asymmetric Migration model = ", ll_split_asymmig

        #calculate theta
        theta_split_asymmig = dadi.Inference.optimal_sfs_scaling(split_asymmig_model, fs)
        theta_split_asymmig = np.around(theta_split_asymmig, 3)
        print "Theta = ", theta_split_asymmig

        #calculate AIC 
        aic_split_asymmig = ( -2*( float(ll_split_asymmig))) + (2*5)
        print "AIC = ", aic_split_asymmig, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_asymmig_model

    elif modelchoice == "secasymmig":
        #####################################
        #Split with secondary contact, asymmetrical gene flow
        #####################################
        print "---------------------------------------------------"
        print "Split with Secondary Contact, Asymmetrical Gene Flow",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_secasymmig = TwoD_models.split_secondary_contact_asym

        #create an extrapolating function 
        func_split_secasymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_secasymmig)

        #simulate the model with the optimized parameters
        split_secasymmig_model = func_split_secasymmig_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secasymmig = dadi.Inference.ll_multinom(split_secasymmig_model, fs)
        ll_split_secasymmig = np.around(ll_split_secasymmig, 3)
        print "likelihood of Split with Secondary Contact and Asymmetric Migration model = ", ll_split_secasymmig

        #calculate theta
        theta_split_secasymmig = dadi.Inference.optimal_sfs_scaling(split_secasymmig_model, fs)
        theta_split_secasymmig = np.around(theta_split_secasymmig, 3)
        print "Theta = ", theta_split_secasymmig

        #calculate AIC 
        aic_split_secasymmig = ( -2*( float(ll_split_secasymmig))) + (2*6)
        print "AIC = ", aic_split_secasymmig, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_secasymmig_model

    elif modelchoice == "secsymmig":
        #####################################
        #Split with secondary contact, symmetrical gene flow
        #####################################
        print "---------------------------------------------------"
        print "Split with Secondary Contact, Symmetrical Gene Flow",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_secsymmig = TwoD_models.split_secondary_contact_sym

        #create an extrapolating function 
        func_split_secsymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_secsymmig)

        #simulate the model with the optimized parameters
        split_secsymmig_model = func_split_secsymmig_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secsymmig = dadi.Inference.ll_multinom(split_secsymmig_model, fs)
        ll_split_secsymmig = np.around(ll_split_secsymmig, 3)
        print "likelihood of Split with Secondary Contact and Symmetric Migration model = ", ll_split_secsymmig

        #calculate theta
        theta_split_secsymmig = dadi.Inference.optimal_sfs_scaling(split_secsymmig_model, fs)
        theta_split_secsymmig = np.around(theta_split_secsymmig, 3)
        print "Theta = ", theta_split_secsymmig

        #calculate AIC 
        aic_split_secsymmig = ( -2*( float(ll_split_secsymmig))) + (2*5)
        print "AIC = ", aic_split_secsymmig, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_secsymmig_model

    elif modelchoice == "ancasymmig":
        #####################################
        #Split with ancient Asymmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Asymmetrical Migration",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_asymancmig = TwoD_models.split_ancient_asymig

        #create an extrapolating function 
        func_split_asymancmig_exec = dadi.Numerics.make_extrap_log_func(func_split_asymancmig)

        #simulate the model with the optimized parameters
        split_asymancmig_model = func_split_asymancmig_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_asymancmig = dadi.Inference.ll_multinom(split_asymancmig_model, fs)
        ll_split_asymancmig = np.around(ll_split_asymancmig, 3)
        print "likelihood of Split with Ancient Symmetrical Migration model = ", ll_split_asymancmig

        #calculate theta
        theta_split_asymancmig = dadi.Inference.optimal_sfs_scaling(split_asymancmig_model, fs)
        theta_split_asymancmig = np.around(theta_split_asymancmig, 3)
        print "Theta = ", theta_split_asymancmig

        #calculate AIC 
        aic_split_asymancmig = ( -2*( float(ll_split_asymancmig))) + (2*6)
        print "AIC = ", aic_split_asymancmig, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_asymancmig_model

    elif modelchoice == "ancsymmig": 
        #####################################
        #Split with ancient symmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Symmetrical Migration",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_symancmig = TwoD_models.split_ancient_symig

        #create an extrapolating function 
        func_split_symancmig_exec = dadi.Numerics.make_extrap_log_func(func_split_symancmig)

        #simulate the model with the optimized parameters
        split_symancmig_model = func_split_symancmig_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_symancmig = dadi.Inference.ll_multinom(split_symancmig_model, fs)
        ll_split_symancmig = np.around(ll_split_symancmig, 3)
        print "likelihood of Split with Ancient Symmetrical Migration model = ", ll_split_symancmig

        #calculate theta
        theta_split_symancmig = dadi.Inference.optimal_sfs_scaling(split_symancmig_model, fs)
        theta_split_symancmig = np.around(theta_split_symancmig, 3)
        print "Theta = ", theta_split_symancmig

        #calculate AIC 
        aic_split_symancmig = ( -2*( float(ll_split_symancmig))) + (2*5)
        print "AIC = ", aic_split_symancmig, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_symancmig_model

    elif modelchoice == "split_no_mig_size":
        #####################################
        #Split with no migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Split with No Migration, size change",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_nomig_size = TwoD_models.split_no_mig_size

        #create an extrapolating function 
        func_split_nomig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_nomig_size)

        #simulate the model with the optimized parameters
        split_nomig_size_model = func_split_nomig_size_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_nomig_size = dadi.Inference.ll_multinom(split_nomig_size_model, fs)
        ll_split_nomig_size = np.around(ll_split_nomig_size, 3)
        print "likelihood of Split with No Migration Size Change model = ", ll_split_nomig_size

        #calculate theta
        theta_split_nomig_size = dadi.Inference.optimal_sfs_scaling(split_nomig_size_model, fs)
        theta_split_nomig_size = np.around(theta_split_nomig_size, 3)
        print "Theta = ", theta_split_nomig_size

        #calculate AIC 
        aic_split_nomig_size = ( -2*( float(ll_split_nomig_size))) + (2*6)
        print "AIC = ", aic_split_nomig_size, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_nomig_size_model

    elif modelchoice == "secasymmig_size":
        #####################################
        #Split with secondary contact, asymmetrical gene flow size change
        #####################################
        print "---------------------------------------------------"
        print "Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_secasymmig_size = TwoD_models.split_secondary_contact_asym_size

        #create an extrapolating function 
        func_split_secasymmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_secasymmig_size)

        #simulate the model with the optimized parameters
        split_secasymmig_size_model = func_split_secasymmig_size_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secasymmig_size = dadi.Inference.ll_multinom(split_secasymmig_size_model, fs)
        ll_split_secasymmig_size = np.around(ll_split_secasymmig_size, 3)
        print "likelihood of Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow = ", ll_split_secasymmig_size

        #calculate theta
        theta_split_secasymmig_size = dadi.Inference.optimal_sfs_scaling(split_secasymmig_size_model, fs)
        theta_split_secasymmig_size = np.around(theta_split_secasymmig_size, 3)
        print "Theta = ", theta_split_secasymmig_size

        #calculate AIC 
        aic_split_secasymmig_size = ( -2*( float(ll_split_secasymmig_size))) + (2*8)
        print "AIC = ", aic_split_secasymmig_size, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_secasymmig_size_model

    elif modelchoice == "secsymmig_size":
        #####################################
        #Split with secondary contact, symmetrical gene flow size change
        #####################################
        print "---------------------------------------------------"
        print "Split with No Mig, Size Change, Secondary Contact, Symmetrical Gene Flow",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_secsymmig_size = TwoD_models.split_secondary_contact_sym_size

        #create an extrapolating function 
        func_split_secsymmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_secsymmig_size)

        #simulate the model with the optimized parameters
        split_secsymmig_size_model = func_split_secsymmig_size_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secsymmig_size = dadi.Inference.ll_multinom(split_secsymmig_size_model, fs)
        ll_split_secsymmig_size = np.around(ll_split_secsymmig_size, 3)
        print "likelihood of Split with No Mig, Size Change, Secondary Contact, Symmetrical Gene Flow = ", ll_split_secsymmig_size

        #calculate theta
        theta_split_secsymmig_size = dadi.Inference.optimal_sfs_scaling(split_secsymmig_size_model, fs)
        theta_split_secsymmig_size = np.around(theta_split_secsymmig_size, 3)
        print "Theta = ", theta_split_secsymmig_size

        #calculate AIC 
        aic_split_secsymmig_size = ( -2*( float(ll_split_secsymmig_size))) + (2*7)
        print "AIC = ", aic_split_secsymmig_size, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_secsymmig_size_model

    
    elif modelchoice == "ancasymmig_size":
        #####################################
        #Split with ancient Asymmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Asymmetrical Migration, Size Change",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_asymancmig_size = TwoD_models.split_ancient_asymmig_size

        #create an extrapolating function 
        func_split_asymancmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_asymancmig_size)

        #simulate the model with the optimized parameters
        split_asymancmig_size_model = func_split_asymancmig_size_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_asymancmig_size = dadi.Inference.ll_multinom(split_asymancmig_size_model, fs)
        ll_split_asymancmig_size = np.around(ll_split_asymancmig_size, 3)
        print "likelihood of Split with Ancient Asymmetrical Migration, Size Change = ", ll_split_asymancmig_size

        #calculate theta
        theta_split_asymancmig_size = dadi.Inference.optimal_sfs_scaling(split_asymancmig_size_model, fs)
        theta_split_asymancmig_size = np.around(theta_split_asymancmig_size, 3)
        print "Theta = ", theta_split_asymancmig_size

        #calculate AIC 
        aic_split_asymancmig_size = ( -2*( float(ll_split_asymancmig_size))) + (2*8)
        print "AIC = ", aic_split_asymancmig_size, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_asymancmig_size_model

            
    elif modelchoice == "ancsymmig_size":
        #####################################
        #Split with ancient symmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Symmetrical Migration, Size Change",'\n','\n'
        print "parameters = {}".format(params)

        #first call a predefined model
        func_split_symancmig_size = TwoD_models.split_ancient_symmig_size

        #create an extrapolating function 
        func_split_symancmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_symancmig_size)

        #simulate the model with the optimized parameters
        split_symancmig_size_model = func_split_symancmig_size_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_symancmig_size = dadi.Inference.ll_multinom(split_symancmig_size_model, fs)
        ll_split_symancmig_size = np.around(ll_split_symancmig_size, 3)
        print "likelihood of Split with Ancient Symmetrical Migration, Size Change model = ", ll_split_symancmig_size

        #calculate theta
        theta_split_symancmig_size = dadi.Inference.optimal_sfs_scaling(split_symancmig_size_model, fs)
        theta_split_symancmig_size = np.around(theta_split_symancmig_size, 3)
        print "Theta = ", theta_split_symancmig_size

        #calculate AIC 
        aic_split_symancmig_size = ( -2*( float(ll_split_symancmig_size))) + (2*7)
        print "AIC = ", aic_split_symancmig_size, '\n', '\n'
        print "---------------------------------------------------", '\n'
        return split_symancmig_size_model
        

#======================================================================================
# Finally, execute model with appropriate arguments
# Two_Pop_Models(pts, fs, params, modelchoice):

# pts = grid choice (list of three numbers, ex. [20,30,40]
# fs = spectrum object name
# params = list of best parameters to start optimizations from
# modelchoice = snm, bottlegrowth, split_nomig, split_symmig, split_asymmig, secasymmig, secsymmig, ancasymmig,
# 				ancsymmig, split_no_mig_size, secasymmig_size, secsymmig_size, ancasymmig_size, ancsymmig_size

#**************
#Input some of the basic arguments here
pts = [50,60,70]
fs = fs_1

print "Creating best models:"

#Please note that NO OPTIMIZATIONS are happening here, so input your final
#parameter estimates from the best scoring rep for each model listed

#**************
#neutral, no divergence: no parameters for this model
params1 = []
best_model1 = Two_Pop_Models(pts, fs, params1, "snm")

#**************
#bottlegrowth: params list must contain 3 values
params2 = [0.0102,0.7498,0.1359]
best_model2 = Two_Pop_Models(pts, fs, params2, "bottlegrowth")

#**************
#split_nomig: params list must contain 3 values
params3 = [3.9198,1.6146,1.7241]
best_model3 = Two_Pop_Models(pts, fs, params3, "split_nomig")

#**************
#split_symmig: params list must contain 4 values
params4 = [8.4227,3.6922,5.6419,0.0115]
best_model4 = Two_Pop_Models(pts, fs, params4, "split_symmig")

#**************
#split_asymmig: params list must contain 5 values
params5 = [5.2636,2.188,3.1079,0.0163,0.0213]
best_model5 = Two_Pop_Models(pts, fs, params5, "split_asymmig")

#**************
#secasymmig: params list must contain 6 values
params6 = [3.8791,1.428,0.8541,1.3946,1.9231,0.0058]
best_model6 = Two_Pop_Models(pts, fs, params6, "secasymmig")

#**************
#secsymmig: params list must contain 5 values
params7 = [5.746,2.5323,0.0352,3.1766,0.3148]
best_model7 = Two_Pop_Models(pts, fs, params7, "secsymmig")

#**************
#ancasymmig: params list must contain 6 values
params8 = [3.8398,1.5847,0.0102,6.4057,0.0307,1.6629]
best_model8 = Two_Pop_Models(pts, fs, params8, "ancasymmig")

#**************
#ancsymmig: params list must contain 5 values
params9 = [8.959,3.9256,0.0113,6.0356,0.0206]
best_model9 = Two_Pop_Models(pts, fs, params9, "ancsymmig")

#**************
#split_no_mig_size: params list must contain 6 values
params10 = [0.159,0.1847,3.1508,1.3121,0.1176,0.5567]
best_model10 = Two_Pop_Models(pts, fs, params10, "split_no_mig_size")

#**************
#secasymmig_size: params list must contain 8 values
params11 = [0.896,0.9287,5.769,1.9553,1.4981,0.7833,0.0203,0.08]
best_model11 = Two_Pop_Models(pts, fs, params11, "secasymmig_size")

#**************
#secsymmig_size: params list must contain 7 values
params12 = [0.2588,0.2473,3.1207,1.2905,0.3018,0.5559,0.0491]
best_model12 = Two_Pop_Models(pts, fs, params12, "secsymmig_size")

#**************
#ancasymmig_size: params list must contain 8 values
params13 = [0.5885,1.0205,5.1647,1.9756,1.4529,0.6219,0.3315,0.2186]
best_model13 = Two_Pop_Models(pts, fs, params13, "ancasymmig_size")

#**************
#ancsymmig_size: params list must contain 7 values
params14 = [2.7982,2.9218,18.7601,7.7038,7.5219,2.213,0.0621]
best_model14 = Two_Pop_Models(pts, fs, params14, "ancsymmig_size")


#======================================================================================
#generic plotting function for data and model comparison

def plotter(model, data, outprefix, model_name):
    print '{0}_{1}.pdf'.format(outprefix,model_name), '\n'
    outname = '{0}_{1}.pdf'.format(outprefix,model_name)
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, data, resid_range = 3)
    fig.savefig(outname)

#======================================================================================
#Execute the plotting function for each of the model sets
#There are four arguments: (model, data, outprefix, model_name)
#1- model object returned from Two_Pop_Models function
#2- frequency spectrum object
#3- prefix for output file name = ('prefix_modelsuffix.pdf')
#4- model suffix for output file name = ('prefix_modelsuffix.pdf')

#**************
#basic arguments
fs = fs_1
outprefix = "Northern_Southern"
#**************
plotter(best_model1, fs, outprefix, "snm")
plotter(best_model2, fs, outprefix, "bottlegrowth")
plotter(best_model3, fs, outprefix, "split_nomig")
plotter(best_model4, fs, outprefix, "split_symmig")
plotter(best_model5, fs, outprefix, "split_asymmig")
plotter(best_model6, fs, outprefix, "secasymmig")
plotter(best_model7, fs, outprefix, "secsymmig")
plotter(best_model8, fs, outprefix, "ancasymmig")
plotter(best_model9, fs, outprefix, "ancsymmig")
plotter(best_model10, fs, outprefix, "split_no_mig_size")
plotter(best_model11, fs, outprefix, "secasymmig_size")
plotter(best_model12, fs, outprefix, "secsymmig_size")
plotter(best_model13, fs, outprefix, "ancasymmig_size")
plotter(best_model14, fs, outprefix, "ancsymmig_size")
#======================================================================================
