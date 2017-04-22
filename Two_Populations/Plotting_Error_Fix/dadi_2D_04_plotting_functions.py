import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import matplotlib.pyplot as plt

'''
usage: python dadi_2D_04_plotting_functions.py

Requires the Models_2D.py script to be in same working directory. This is where
all the population model functions are stored for this script. 

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
daniel.portik@uta.edu
April 2017
'''

#===========================================================================
#get snps file 

#**************
snps1 = "/FULL PATH TO /dadi_2pops_Cameroon_South_snps.txt"

#Create python dictionary from snps file
dd1 = dadi.Misc.make_data_dict(snps1)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["Cameroon", "South"]
#projection sizes, in ALLELES not individuals
proj_1 = [16,28]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)

print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================
#create function to run models with fixed user input parameters (presumably from best runs)

def Two_Pop_Plot(pts, fs, model_name, params):

    if model_name == "no_divergence":
        #####################################
        #No Divergence
        #####################################
        print "---------------------------------------------------"
        print "No Divergence",'\n'

        #first call a predefined model
        model_call = Models_2D.no_divergence
        params = []
        print "parameter set = [none]"
        
        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

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
        aic = ( -2*( float(ll))) + (2*3)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "no_mig":
        #####################################
        # Divergence with no migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration",'\n'

        #first call a predefined model
        model_call = Models_2D.no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "base parameters = ", params

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
        aic = ( -2*( float(ll))) + (2*3)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "sym_mig":
        #####################################
        #Divergence with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration",'\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "base parameters = ", params
        
        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*4)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "asym_mig":
        #####################################
        #Divergence with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration",'\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*5)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "anc_sym_mig":
        #####################################
        #Divergence with ancient symmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration",'\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*5)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "anc_asym_mig":
        #####################################
        #Divergence with ancient asymmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration",'\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "sec_contact_sym_mig":
        #####################################
        #Divergence and symmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact",'\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*5)
        print "AIC = ", aic, '\n', '\n'
        return sim_model
        
        
    elif model_name == "sec_contact_asym_mig":
        #####################################
        #Divergence and asymmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact",'\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params
        
        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


    ##########################################################################
    #Similar models but with size change allowed
    ##########################################################################
    
    elif model_name == "no_mig_size":
        #####################################
        #Divergence with no migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration, size change",'\n'

        #first call a predefined model
        model_call = Models_2D.no_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model
        
        
    elif model_name == "sym_mig_size":
        #####################################
        #Divergence with symmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, size change",'\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*7)
        print "AIC = ", aic, '\n', '\n'
        return sim_model
        

    elif model_name == "asym_mig_size":
        #####################################
        #Divergence with asymmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, size change",'\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*8)
        print "AIC = ", aic, '\n', '\n'
        return sim_model
 
 
    elif model_name == "anc_sym_mig_size":
        #####################################
        #Divergence with ancient symmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration, size change",'\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*7)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "anc_asym_mig_size":
        #####################################
        #Divergence with ancient asymmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration, size change",'\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*8)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


    elif model_name == "sec_contact_sym_mig_size":
        #####################################
        #Divergence and symmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change",'\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*7)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

           
    elif model_name == "sec_contact_asym_mig_size":
        #####################################
        #Divergence and asymmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change",'\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        print "base parameters = ", params

        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*8)
        print "AIC = ", aic, '\n', '\n'
        return sim_model
        

#======================================================================================
# Finally, execute model with appropriate arguments, including name for returned model
# returned_model = Two_Pop_Plot(pts, fs, model_name, params):

# returned_model:  name the simulated model returned by function (to store for plotting later)

# pts:  grid choice (list of three numbers, ex. [20,30,40]

# fs:  spectrum object name

# model_name:  "no_divergence", "no_mig", "sym_mig", "asym_mig", "anc_sym_mig", "anc_asym_mig",
#        "sec_contact_sym_mig", "sec_contact_asym_mig", "no_mig_size", "sym_mig_size",
#        "asym_mig_size", "anc_sym_mig_size", "anc_asym_mig_size", "sec_contact_sym_mig_size",
#        "sec_contact_asym_mig_size"

# params:  list of best parameters to start optimizations from

#===========================================================================
# enter best param values for each model here, presumably you will get these
# from the outputs of the previous script

#************** "no_divergence"
#leave blank, no parameters
no_divergence_params = []

#************** "no_mig"
#3 Values
no_mig_params = [1.6827,3.8731,2.2595]

#************** "sym_mig"
#4 Values
sym_mig_params = [2.6261,5.932,0.0107,4.5966]

#************** "asym_mig"
#5 Values
asym_mig_params = [2.432,5.4933,0.0121,0.0109,4.1582]

#************** "anc_sym_mig"
#5 Values
anc_sym_mig_params = [1.7253,3.8058,0.0251,2.4635,0.2359]

#************** "anc_asym_mig"
#6 Values
anc_asym_mig_params = [1.6635,4.3311,0.0557,0.6685,0.0815,2.3774]

#************** "sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [2.7735,6.5051,0.0218,4.3453,0.4201]

#************** "sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [1.4153,3.1098,0.2271,0.0837,1.9919,0.0014]

#************** "no_mig_size"
#6 Values
no_mig_size_params = [0.2166,0.1304,1.8324,3.7265,0.179,0.4558]

#************** "sym_mig_size"
#7 Values
sym_mig_size_params = [0.3127,0.2534,1.456,3.8703,0.0243,0.4483,0.3563]

#************** "asym_mig_size"
#8 Values
asym_mig_size_params = [0.3101,0.3323,1.2428,2.0772,0.0985,0.0285,0.9026,0.4849]

#************** "anc_sym_mig_size"
#7 Values
anc_sym_mig_size_params = [0.2692,0.1682,2.0483,3.8678,0.6095,0.5148,0.6448]

#************** "anc_asym_mig_size"
#8 Values
anc_asym_mig_size_params = [0.4814,0.5822,3.9615,8.851,0.4163,0.1136,9.2139,1.1795]

#************** "sec_contact_sym_mig_size"
#7 Values
sec_contact_sym_mig_size_params = [0.3245,0.4643,2.3541,5.5444,0.0199,1.0492,1.0659]

#************** "sec_contact_asym_mig_size"
#8 Values
sec_contact_asym_mig_size_params = [1.1401,0.1186,0.2519,2.4771,0.8829,0.0836,0.4183,0.26]


#===========================================================================
#**************
#Input some of the basic reusable arguments here
pts = [50,60,70]
fs = fs_1
outfile = "N_v_S"

#===========================================================================
# returned_model = Two_Pop_Plot(pts, fs, model_name, params):

# Each model is executed with one replicate using fixed parameter values.
# The simulated model is stored as an object to be called on in the
# subsequent plotting function.

no_divergence = Two_Pop_Plot(pts, fs, "no_divergence", no_divergence_params)
no_mig = Two_Pop_Plot(pts, fs, "no_mig", no_mig_params)
sym_mig = Two_Pop_Plot(pts, fs, "sym_mig", sym_mig_params)
asym_mig = Two_Pop_Plot(pts, fs, "asym_mig", asym_mig_params)
anc_sym_mig = Two_Pop_Plot(pts, fs, "anc_sym_mig", anc_sym_mig_params)
anc_asym_mig = Two_Pop_Plot(pts, fs, "anc_asym_mig", anc_asym_mig_params)
sec_contact_sym_mig = Two_Pop_Plot(pts, fs, "sec_contact_sym_mig", sec_contact_sym_mig_params)
sec_contact_asym_mig = Two_Pop_Plot(pts, fs, "sec_contact_asym_mig", sec_contact_asym_mig_params)
no_mig_size = Two_Pop_Plot(pts, fs, "no_mig_size", no_mig_size_params)
sym_mig_size = Two_Pop_Plot(pts, fs, "sym_mig_size", sym_mig_size_params)
asym_mig_size = Two_Pop_Plot(pts, fs, "asym_mig_size", asym_mig_size_params)
anc_sym_mig_size = Two_Pop_Plot(pts, fs, "anc_sym_mig_size", anc_sym_mig_size_params)
anc_asym_mig_size = Two_Pop_Plot(pts, fs, "anc_asym_mig_size", anc_asym_mig_size_params)
sec_contact_sym_mig_size = Two_Pop_Plot(pts, fs, "sec_contact_sym_mig_size", sec_contact_sym_mig_size_params)
sec_contact_asym_mig_size = Two_Pop_Plot(pts, fs, "sec_contact_asym_mig_size", sec_contact_asym_mig_size_params)

#======================================================================================
#write a plotting function for data and model comparison

def plot_all(sim_model, data, outfile, model_name):
    print '{0}_{1}.pdf'.format(outfile,model_name), '\n'
    outname = '{0}_{1}.pdf'.format(outfile,model_name)
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(sim_model, data, resid_range = 3, vmin= 0.005)
    fig.savefig(outname)
    
#======================================================================================
# plot_all(sim_model, fs, outfile, model_name)

# sim_model:  the simulated model returned from function Two_Pop_Plot

# fs:  spectrum object name (ie the "data")

# outfile:  same prefix for output naming of pdf files -> "{0}_{1}.pdf".format(outfile,model_name)

# model_name:  "no_divergence", "no_mig", "sym_mig", "asym_mig", "anc_sym_mig", "anc_asym_mig",
#        "sec_contact_sym_mig", "sec_contact_asym_mig", "no_mig_size", "sym_mig_size",
#        "asym_mig_size", "anc_sym_mig_size", "anc_asym_mig_size", "sec_contact_sym_mig_size",
#        "sec_contact_asym_mig_size"

# call plot function for each model 

plot_all(no_divergence, fs, outfile, "no_divergence")
plot_all(no_mig, fs, outfile, "no_mig")
plot_all(sym_mig, fs, outfile, "sym_mig")
plot_all(asym_mig, fs, outfile, "asym_mig")
plot_all(anc_sym_mig, fs, outfile, "anc_sym_mig")
plot_all(anc_asym_mig, fs, outfile, "anc_asym_mig")
plot_all(sec_contact_sym_mig, fs, outfile, "sec_contact_sym_mig")
plot_all(sec_contact_asym_mig, fs, outfile, "sec_contact_asym_mig")
plot_all(no_mig_size, fs, outfile, "no_mig_size")
plot_all(sym_mig_size, fs, outfile, "sym_mig_size")
plot_all(asym_mig_size, fs, outfile, "asym_mig_size")
plot_all(anc_sym_mig_size, fs, outfile, "anc_sym_mig_size")
plot_all(anc_asym_mig_size, fs, outfile, "anc_asym_mig_size")
plot_all(sec_contact_sym_mig_size, fs, outfile, "sec_contact_sym_mig_size")
plot_all(sec_contact_asym_mig_size, fs, outfile, "sec_contact_asym_mig_size")

#======================================================================================
