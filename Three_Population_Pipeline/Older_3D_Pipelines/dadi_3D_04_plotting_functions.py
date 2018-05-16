import sys
import os
import numpy as np
import dadi
from datetime import datetime
import Models_3D
import pylab
import matplotlib.pyplot as plt

'''
usage: python dadi_3D_04_plotting_functions.py

Requires the Models_3D.py script to be in same working directory.
This is where all the population model functions are stored. 
This is written for the models specifically found in that script. 

Script will simulate the model using a fixed set of input parameters. That is,
no optimizations will occur here, and the assumption is you've searched
thoroughly for the parameter set with the highest likelihood value. The goal
of this script is to produce plots of the 3D-JSFS for the data and model,
plus the residuals. This will take the form of three two-population comparisons (2D-JSFS plots)
per plot. A pdf file will be written for each model to the current working directory.

Requires user to edit sections of code marked with #**************

You'll absolutely need to provide the path to your SNPs input file
along with your specific projections and population labels. 

###########################################################
**Note, if you see this error when plotting:
"ValueError: Data has no positive values, and therefore can not be log-scaled."
You will need to change the vmin in the plotting routine from None to something 0<1.

I have changed vmin to 0.005-0.01 with good results. So, at the bottom of this script:

dadi.Plotting.plot_3d_comp_multinom(model, data, resid_range = 3)

becomes:

dadi.Plotting.plot_3d_comp_multinom(model, data, resid_range = 3, vmin=0.005)

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
daniel.portik@uta.edu
April 2017
'''
t_begin = datetime.now()

#===========================================================================
#get snps file 
#**************
snps1 = "/FULL PATH TO/dadi_3pops_CVLS_CVLN_Cross_snps.txt"

#Create python dictionary from snps file
dd1 = dadi.Misc.make_data_dict(snps1)

#**************
#projection sizes, in ALLELES not individuals
proj_1 = [14,30,18]
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=['CVLS','CVLN','Cross']

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)

print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'


#======================================================================================
#create function to run models with fixed user input parameters (presumably from best runs)

def Three_Pop_Plot(pts, fs, model_name, params):
    
    ######################################################################################################################
    #Basic models
    ######################################################################################################################
    
    if model_name == "split_nomig":
        
        #####################################
        #Split with no migration
        #####################################
        print "---------------------------------------------------"
        print "Split with No Migration",'\n'
        
        #first call a predefined model
        model_call = Models_3D.split_nomig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params
        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        
        return sim_model


    elif model_name == "split_symmig_all":
        
        #####################################
        #Split with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Symmetric Migration",'\n'
        
        #first call a predefined model
        model_call = Models_3D.split_symmig_all

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*10)
        print "AIC = ", aic, '\n', '\n'

        return sim_model


    elif model_name == "split_symmig_adjacent":
        
        #####################################
        #Split with adjacent symmetric migration 
        #####################################
        print "---------------------------------------------------"
        print "Split with Adjacent Symmetric Migration",'\n'
        
        #first call a predefined model
        model_call = Models_3D.split_symmig_adjacent

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*9)
        print "AIC = ", aic, '\n', '\n'
        
        return sim_model

    
    ######################################################################################################################
    #Refugial models
    ######################################################################################################################

    elif model_name == "refugia_1":
        
        #####################################
        #Refugia 1 with secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Refugia 1 with secondary contact",'\n'

        #first call a predefined model
        model_call = Models_3D.refugia_1

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*9)
        print "AIC = ", aic, '\n', '\n'
            
        return sim_model

    elif model_name == "refugia_2":

        #####################################
        #Refugia 2 with secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Refugia 2 with secondary contact",'\n'
        
        #first call a predefined model
        model_call = Models_3D.refugia_2

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*8)
        print "AIC = ", aic, '\n', '\n'
            
        return sim_model

    elif model_name == "refugia_3":

        #####################################
        #Refugia with secondary contact, shortest isolation
        #####################################
        print "---------------------------------------------------"
        print "Refugia 3 with secondary contact",'\n'
        
        #first call a predefined model
        model_call = Models_3D.refugia_3

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*10)
        print "AIC = ", aic, '\n', '\n'
            
        return sim_model

    
    ######################################################################################################################
    #Ancient migration models
    ######################################################################################################################

    elif model_name == "ancmig_3":

        #####################################
        #Ancient migration 3, long isolation
        #####################################
        print "---------------------------------------------------"
        print "Ancient migration 3",'\n'

        #first call a predefined model
        model_call = Models_3D.ancmig_3

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*8)
        print "AIC = ", aic, '\n', '\n'

        return sim_model

    elif model_name == "ancmig_2":

        #####################################
        #Ancient migration 2, shorter isolation
        #####################################
        print "---------------------------------------------------"
        print "Ancient migration 2",'\n'

        #first call a predefined model
        model_call = Models_3D.ancmig_2

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*7)
        print "AIC = ", aic, '\n', '\n'
            
        return sim_model

        
    elif model_name == "ancmig_1":

        #####################################
        #Ancient migration 1, shortest isolation
        #####################################
        print "---------------------------------------------------"
        print "Ancient migration 1",'\n'
        
        #first call a predefined model
        model_call = Models_3D.ancmig_1

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*10)
        print "AIC = ", aic, '\n', '\n'

        return sim_model


        
#======================================================================================
# Finally, execute model with appropriate arguments, including name for returned model
# returned_model = Three_Pop_Plot(pts, fs, model_name, params):

# returned_model:  name the simulated model returned by function (to store for plotting later)

# pts:  grid choice (list of three numbers, ex. [20,30,40]

# fs:  spectrum object name

# model_name:  "split_nomig", "split_symmig_all", "split_symmig_adjacent", "refugia_1",
#        "refugia_2", "refugia_3", "ancmig_3", "ancmig_2", "ancmig_1", 

# params:  list of best parameters to start optimizations from
#		ex. some_params = [4.787,0.465,7.071,1.879,0.181,0.635]


#===========================================================================
# enter best param values for each model here, presumably you will get these
# from the outputs of the previous script

#************** "split_nomig"
# 6 Values
split_nomig_params = [2.1127,0.9016,7.4879,0.8309,0.0737,0.5939]

#************** "split_symmig_all"
# 10 Values
split_symmig_all_params = [0.6047,1.8755,2.7982,0.1294,1.1115,0.4329,1.9536,0.5993,1.5957,0.3824]

#************** "split_symmig_adjacent"
# 9 Values
split_symmig_adjacent_params = [2.911,6.4277,4.4721,1.6473,3.1932,0.2043,0.2706,6.2816,1.7949]

#************** "refugia_1"
# 9 Values
refugia_1_params = [1.8977,1.1783,3.6668,1.2239,0.4916,0.1053,0.1639,0.8435,0.1292]

#************** "refugia_2"
# 8 Values
refugia_2_params = [1.5786,15.1636,1.9566,0.3874,0.4503,0.467,0.2157,0.7116]

#************** "refugia_3"
# 10 Values
refugia_3_params = [6.3651,9.2116,16.4184,0.8546,0.9772,0.1257,0.1512,0.7008,2.4985,1.5766]

#************** "ancmig_3"
# 8 Values
ancmig_3_params = [5.0892,11.3673,4.7975,1.3108,0.3493,3.0766,0.2625,0.2897]

#************** "ancmig_2"
# 7 Values
ancmig_2_params = [0.6107,22.3261,5.7073,0.1822,6.9999,0.6566,0.1537]

#************** "ancmig_1"
# 10 Values
ancmig_1_params = [2.193,0.323,5.482,3.195,0.170,1.745,0.134,2.742,5.314,0.184]



#**************
#Input some of the basic reusable arguments here
pts = [50,60,70]
fs = fs_1
outfile = "CVLS_CVLN_Cross"


#===========================================================================
# returned_model = Three_Pop_Plot(pts, fs, model_name, params):

# Each model is executed with one replicate using fixed parameter values.
# The simulated model is stored as an object to be called on in the
# subsequent plotting function.

#**************
#Models from the Models_3D.py script
split_nomig = Three_Pop_Plot(pts, fs, "split_nomig", split_nomig_params)
split_symmig_all = Three_Pop_Plot(pts, fs, "split_symmig_all", split_symmig_all_params)
split_symmig_adjacent = Three_Pop_Plot(pts, fs, "split_symmig_adjacent", split_symmig_adjacent_params)
refugia_1 = Three_Pop_Plot(pts, fs, "refugia_1", refugia_1_params)
refugia_2 = Three_Pop_Plot(pts, fs, "refugia_2", refugia_2_params)
refugia_3 = Three_Pop_Plot(pts, fs, "refugia_3", refugia_3_params)
ancmig_3 = Three_Pop_Plot(pts, fs, "ancmig_3", ancmig_3_params)
ancmig_2 = Three_Pop_Plot(pts, fs, "ancmig_2", ancmig_2_params)
ancmig_1 = Three_Pop_Plot(pts, fs, "ancmig_1", ancmig_1_params)

#======================================================================================
#write a plotting function for data and model comparison

def plot_all(sim_model, data, outfile, model_name):
    print '{0}_{1}.pdf'.format(outfile,model_name), '\n'
    outname = '{0}_{1}.pdf'.format(outfile,model_name)
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_3d_comp_multinom(sim_model, data, resid_range = 3)
    fig.savefig(outname)
    
#======================================================================================
# plot_all(sim_model, fs, outfile, model_name)

# sim_model:  the simulated model returned from function Three_Pop_Plot

# fs:  spectrum object name (ie the "data")

# outfile:  same prefix for output naming of pdf files -> "{0}_{1}.pdf".format(outfile,model_name)

# model_name:  "split_nomig", "split_symmig_all", "split_symmig_adjacent", "refugia_1",
#        "refugia_2", "refugia_3", "ancmig_3", "ancmig_2", "ancmig_1", "refugia_1_size",
#        "refugia_2_size", "refugia_3_size", "island_exchange", "island_exchange_size","refugia2_island",
#        "refugia3_island", "island_exchange_size", "refugia2_island_size", "refugia3_island_size",
#		 "island_exchange_variation_size"

# call plot function for each model

plot_all(split_nomig, fs, outfile, "split_nomig")
plot_all(split_symmig_all, fs, outfile, "split_symmig_all")
plot_all(split_symmig_adjacent, fs, outfile, "split_symmig_adjacent")
plot_all(refugia_1, fs, outfile, "refugia_1")
plot_all(refugia_2, fs, outfile, "refugia_2")
plot_all(refugia_3, fs, outfile, "refugia_3")
plot_all(ancmig_3, fs, outfile, "ancmig_3")
plot_all(ancmig_2, fs, outfile, "ancmig_2")
plot_all(ancmig_1, fs, outfile, "ancmig_1")


#===========================================================================
#clock it!

t_finish = datetime.now()
elapsed = t_finish - t_begin

print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'

#===========================================================================
