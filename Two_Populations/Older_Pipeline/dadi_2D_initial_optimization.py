import sys
import os
import numpy as np
import dadi
import matplotlib
import TwoD_models
'''
usage: python dadi_2D_initial_optimization.py

Requires the TwoD_models.py script to be in same working directory. This is where
all the population model functions are stored. This is written for the
models specifically found in that script. 

Script will perform initial optimizations from multiple starting points using a
perturbed set of basic starting values for parameters. The output is written as 
"{outfile}_INITIAL_model_results.txt", and the best scoring replicate for 
each model can be inspected. The parameter values of these models can be used
to run the dadi_2D_final_optimization.py script, which requires input values
for the parameters of each model to start the optimizations.

Requires user to edit sections of code marked with #**************


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
#you should test a few projections and see what gives the best results
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
#create function to take in grid size and freq spectrum for a population, then run all 
#models 'x' times from optimized parameters from heavily perturbed basic starting values
#write output of parameter optimizations, theta, likelihood, and AIC for each replicate for 
#each model to file

def Two_Pop_Models(pts, fs, outfile, x):
    print '\n',"============================================================================"
    print "Beginning analysis of {}".format(outfile)
    print "============================================================================"

    #create output file
    outname = outfile+"INITIAL_model_results.txt"
    fh_out = open(outname, 'a')
    fh_out.close()
    
    #variable to control number of loops per model (1 to x)
    x = int(x) + int(1)
    #set variable to control maxiter argument
    y = int(100)
    
    
    fh_out = open(outname, 'a')
    #####################################
    #Neutral model, no split
    #####################################
    print "---------------------------------------------------"
    print "Neutral Model, No Divergence",'\n','\n'
    fh_out.write("Neutral model, no split:"+'\n')

    #first call a predefined model
    func_snm = TwoD_models.snm

    #create empty parameter list for neutral model (there are no params here)
    no_params = []
    print "parameter set = [none]"
    fh_out.write("parameter set = [none]"+'\n')

    #create an extrapolating function 
    func_snm_exec = dadi.Numerics.make_extrap_log_func(func_snm)
    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "neutral no split replicate {}:".format(i)
        
        #run the model with the extrapolation call and all our parameters
        snm_model = func_snm_exec(no_params, fs.sample_sizes, pts)

        #calculate likelihood
        ll_neutral = dadi.Inference.ll_multinom(snm_model, fs)
        ll_neutral = np.around(ll_neutral, 3)
        print "likelihood of Neutral Model, No Divergence model = ", ll_neutral
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_neutral)+'\n')
        
        #calculate theta
        theta_snm = dadi.Inference.optimal_sfs_scaling(snm_model, fs)
        theta_snm = np.around(theta_snm, 3)
        print "Theta = ", theta_snm
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_snm)+'\n')

        #calculate AIC 
        #AIC = -2 ( ln ( likelihood )) + 2 K
        aic_neutral = ( -2*( float(ll_neutral)))
        print "AIC = ", aic_neutral, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_neutral)+'\n'+'\n')
        
    print "---------------------------------------------------", '\n'
    fh_out.close()

    
    fh_out = open(outname, 'a') 
    #####################################
    #bottlegrowth model, no split
    #####################################
    print "---------------------------------------------------"
    print "Bottlegrowth Model, No Divergence",'\n','\n'
    fh_out.write("Bottlegrowth Model, No Divergence:"+'\n')

    #first call a predefined model
    func_bottle = TwoD_models.bottlegrowth

    #create an extrapolating function 
    func_bottle_exec = dadi.Numerics.make_extrap_log_func(func_bottle)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0]
    upper_bound = [30, 30, 10]
    print "parameter set = [nu1, nu2, T]"
    fh_out.write("parameter set = [nu1, nu2, T]"+'\n')
    
    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "bottlegrowth no split replicate {}:".format(i)
        
        #perturb initial guesses
        param_bottle = dadi.Misc.perturb_params([1,1,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_bottle_opt = dadi.Inference.optimize_log_fmin(param_bottle, fs, func_bottle_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_bottle_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_bottle_opt)+'\n')


        #simulate the model with the optimized parameters
        bottle_model = func_bottle_exec(param_bottle_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_bottle = dadi.Inference.ll_multinom(bottle_model, fs)
        ll_bottle = np.around(ll_bottle, 3)
        print "likelihood of Bottlegrowth Model, No Divergence model = ", ll_bottle
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_bottle)+'\n')
        
        #calculate theta
        theta_bottle = dadi.Inference.optimal_sfs_scaling(bottle_model, fs)
        theta_bottle = np.around(theta_bottle, 3)
        print "Theta = ", theta_bottle
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_bottle)+'\n')

        #calculate AIC 
        #AIC = -2 ( ln ( likelihood )) + 2 K
        aic_bottle = ( -2*( float(ll_bottle))) + (2*3)
        print "AIC = ", aic_bottle, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_bottle)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()
    
    
    fh_out = open(outname, 'a') 
    #####################################
    #Split with no migration
    #####################################
    print "---------------------------------------------------"
    print "Split with No Migration",'\n','\n'
    fh_out.write("Split with No Migration"+'\n')

    #first call a predefined model
    func_split_nomig = TwoD_models.split_no_mig

    #create an extrapolating function 
    func_split_nomig_exec = dadi.Numerics.make_extrap_log_func(func_split_nomig)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0]
    upper_bound = [30, 30, 10]
    print "parameter set = [nu1, nu2, T]"
    fh_out.write("parameter set = [nu1, nu2, T]"+'\n')
    
    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_nomig replicate {}:".format(i)
        
        #perturb initial guesses
        param_split_nomig = dadi.Misc.perturb_params([1,1,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_nomig_opt = dadi.Inference.optimize_log_fmin(param_split_nomig, fs, func_split_nomig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_nomig_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_nomig_opt)+'\n')

        #simulate the model with the optimized parameters
        split_nomig_model = func_split_nomig_exec(param_split_nomig_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_nomig = dadi.Inference.ll_multinom(split_nomig_model, fs)
        ll_split_nomig = np.around(ll_split_nomig, 3)
        print "likelihood of Split with No Migration model = ", ll_split_nomig
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_nomig)+'\n')

        #calculate theta
        theta_split_nomig = dadi.Inference.optimal_sfs_scaling(split_nomig_model, fs)
        theta_split_nomig = np.around(theta_split_nomig, 3)
        print "Theta = ", theta_split_nomig
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_nomig)+'\n')

        #calculate AIC 
        aic_split_nomig = ( -2*( float(ll_split_nomig))) + (2*3)
        print "AIC = ", aic_split_nomig, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_nomig)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()

    
    fh_out = open(outname, 'a') 
    #####################################
    #Split with symmetric migration
    #####################################
    print "---------------------------------------------------"
    print "Split with Symmetric Migration",'\n','\n'
    fh_out.write("Split with Symmetric Migration:"+'\n')

    #first call a predefined model
    func_split_symmig = TwoD_models.split_mig

    #create an extrapolating function 
    func_split_symmig_exec = dadi.Numerics.make_extrap_log_func(func_split_symmig)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0, 0.01]
    upper_bound = [30, 30, 10, 20]
    print "parameter set = [nu1, nu2, T, m]"
    fh_out.write("parameter set = [nu1, nu2, T, m]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_symmig replicate {}:".format(i)
        
        #perturb initial guesses
        param_split_symmig = dadi.Misc.perturb_params([1,1,0.8,1], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_symmig_opt = dadi.Inference.optimize_log_fmin(param_split_symmig, fs, func_split_symmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_symmig_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_symmig_opt)+'\n')

        #simulate the model with the optimized parameters
        split_symmig_model = func_split_symmig_exec(param_split_symmig_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_symmig = dadi.Inference.ll_multinom(split_symmig_model, fs)
        ll_split_symmig = np.around(ll_split_symmig, 3)
        print "likelihood of Split with Symmetric Migration model = ", ll_split_symmig
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_symmig)+'\n')

        #calculate theta
        theta_split_symmig = dadi.Inference.optimal_sfs_scaling(split_symmig_model, fs)
        theta_split_symmig = np.around(theta_split_symmig, 3)
        print "Theta = ", theta_split_symmig
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_symmig)+'\n')

        #calculate AIC 
        aic_split_symmig = ( -2*( float(ll_split_symmig))) + (2*4)
        print "AIC = ", aic_split_symmig, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_symmig)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()

    
    fh_out = open(outname, 'a') 
    #####################################
    #Split with asymmetric migration
    #####################################
    print "---------------------------------------------------"
    print "Split with Asymmetric Migration",'\n','\n'
    fh_out.write("Split with Asymmetric Migration:"+'\n')

    #first call a predefined model
    func_split_asymmig = TwoD_models.split_asym_mig

    #create an extrapolating function 
    func_split_asymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_asymmig)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0, 0.01, 0.01]
    upper_bound = [30, 30, 10, 20, 20]
    print "parameter set = [nu1, nu2, T, m12, m21]"
    fh_out.write("parameter set = [nu1, nu2, T, m12, m21]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_asymmig replicate {}:".format(i)
        
        #perturb initial guesses
        param_split_asymmig = dadi.Misc.perturb_params([1,1,1,0.8,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_asymmig_opt = dadi.Inference.optimize_log_fmin(param_split_asymmig, fs, func_split_asymmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_asymmig_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_asymmig_opt)+'\n')

        #simulate the model with the optimized parameters
        split_asymmig_model = func_split_asymmig_exec(param_split_asymmig_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_asymmig = dadi.Inference.ll_multinom(split_asymmig_model, fs)
        ll_split_asymmig = np.around(ll_split_asymmig, 3)
        print "likelihood of Split with Asymmetric Migration model = ", ll_split_asymmig
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_asymmig)+'\n')

        #calculate theta
        theta_split_asymmig = dadi.Inference.optimal_sfs_scaling(split_asymmig_model, fs)
        theta_split_asymmig = np.around(theta_split_asymmig, 3)
        print "Theta = ", theta_split_asymmig
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_asymmig)+'\n')

        #calculate AIC 
        aic_split_asymmig = ( -2*( float(ll_split_asymmig))) + (2*5)
        print "AIC = ", aic_split_asymmig, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_asymmig)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()

    
    fh_out = open(outname, 'a') 
    #####################################
    #Split with secondary contact, asymmetrical gene flow
    #####################################
    print "---------------------------------------------------"
    print "Split with Secondary Contact, Asymmetrical Gene Flow",'\n','\n'
    fh_out.write("Split with Secondary Contact, Asymmetrical Gene Flow:"+'\n')
    
    #first call a predefined model
    func_split_secasymmig = TwoD_models.split_secondary_contact_asym

    #create an extrapolating function 
    func_split_secasymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_secasymmig)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
    upper_bound = [30, 30, 20, 20, 10, 10]
    print "parameter set = [nu1, nu2, m12, m21, T, Tsc]"
    fh_out.write("parameter set = [nu1, nu2, m12, m21, T, Tsc]"+'\n')
    
    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_secasymmig replicate {}:".format(i)
        #perturb initial guesses
        param_split_secasymmig = dadi.Misc.perturb_params([1,1,0.8,0.8,1,1], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_secasymmig_opt = dadi.Inference.optimize_log_fmin(param_split_secasymmig, fs, func_split_secasymmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_secasymmig_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_secasymmig_opt)+'\n')

        #simulate the model with the optimized parameters
        split_secasymmig_model = func_split_secasymmig_exec(param_split_secasymmig_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secasymmig = dadi.Inference.ll_multinom(split_secasymmig_model, fs)
        ll_split_secasymmig = np.around(ll_split_secasymmig, 3)
        print "likelihood of Split with Secondary Contact and Asymmetric Migration model = ", ll_split_secasymmig
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_secasymmig)+'\n')

        #calculate theta
        theta_split_secasymmig = dadi.Inference.optimal_sfs_scaling(split_secasymmig_model, fs)
        theta_split_secasymmig = np.around(theta_split_secasymmig, 3)
        print "Theta = ", theta_split_secasymmig
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_secasymmig)+'\n')

        #calculate AIC 
        aic_split_secasymmig = ( -2*( float(ll_split_secasymmig))) + (2*6)
        print "AIC = ", aic_split_secasymmig, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_secasymmig)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()


    fh_out = open(outname, 'a') 
    #####################################
    #Split with secondary contact, symmetrical gene flow
    #####################################
    print "---------------------------------------------------"
    print "Split with Secondary Contact, Symmetrical Gene Flow",'\n','\n'
    fh_out.write("Split with Secondary Contact, Symmetrical Gene Flow:"+'\n')
    
    #first call a predefined model
    func_split_secsymmig = TwoD_models.split_secondary_contact_sym

    #create an extrapolating function 
    func_split_secsymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_secsymmig)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0, 0]
    upper_bound = [30, 30, 20, 10, 10]
    print "parameter set = [nu1, nu2, m, T, Tsc]"
    fh_out.write("parameter set = [nu1, nu2, m1, T, Tsc]"+'\n')
    
    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_secsymmig replicate {}:".format(i)
        #perturb initial guesses
        param_split_secsymmig = dadi.Misc.perturb_params([1,1,0.8,1,1], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_secsymmig_opt = dadi.Inference.optimize_log_fmin(param_split_secsymmig, fs, func_split_secsymmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_secsymmig_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_secsymmig_opt)+'\n')

        #simulate the model with the optimized parameters
        split_secsymmig_model = func_split_secsymmig_exec(param_split_secsymmig_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secsymmig = dadi.Inference.ll_multinom(split_secsymmig_model, fs)
        ll_split_secsymmig = np.around(ll_split_secsymmig, 3)
        print "likelihood of Split with Secondary Contact and Symmetric Migration model = ", ll_split_secsymmig
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_secsymmig)+'\n')

        #calculate theta
        theta_split_secsymmig = dadi.Inference.optimal_sfs_scaling(split_secsymmig_model, fs)
        theta_split_secsymmig = np.around(theta_split_secsymmig, 3)
        print "Theta = ", theta_split_secsymmig
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_secsymmig)+'\n')

        #calculate AIC 
        aic_split_secsymmig = ( -2*( float(ll_split_secsymmig))) + (2*5)
        print "AIC = ", aic_split_secsymmig, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_secsymmig)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()

    
    fh_out = open(outname, 'a') 
    #####################################
    #Split with ancient Asymmetrical migration
    #####################################
    print "---------------------------------------------------"
    print "Split with Ancient Asymmetrical Migration",'\n','\n'
    fh_out.write("Split with Ancient Asymmetrical Migration:"+'\n')

    #first call a predefined model
    func_split_asymancmig = TwoD_models.split_ancient_asymig

    #create an extrapolating function 
    func_split_asymancmig_exec = dadi.Numerics.make_extrap_log_func(func_split_asymancmig)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
    upper_bound = [30, 30, 20, 20, 10, 10]
    print "parameter set = [nu1, nu2, m12, m21, T, Tam]"
    fh_out.write("parameter set = [nu1, nu2, m12, m21, T, Tam]"+'\n')
    
    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_asymancmig replicate {}:".format(i)
        
        #perturb initial guesses
        param_split_asymancmig = dadi.Misc.perturb_params([1,1,0.8,0.8,1,1], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_asymancmig_opt = dadi.Inference.optimize_log_fmin(param_split_asymancmig, fs, func_split_asymancmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_asymancmig_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_asymancmig_opt)+'\n')

        #simulate the model with the optimized parameters
        split_asymancmig_model = func_split_asymancmig_exec(param_split_asymancmig_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_asymancmig = dadi.Inference.ll_multinom(split_asymancmig_model, fs)
        ll_split_asymancmig = np.around(ll_split_asymancmig, 3)
        print "likelihood of Split with Ancient Symmetrical Migration model = ", ll_split_asymancmig
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_asymancmig)+'\n')

        #calculate theta
        theta_split_asymancmig = dadi.Inference.optimal_sfs_scaling(split_asymancmig_model, fs)
        theta_split_asymancmig = np.around(theta_split_asymancmig, 3)
        print "Theta = ", theta_split_asymancmig
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_asymancmig)+'\n')

        #calculate AIC 
        aic_split_asymancmig = ( -2*( float(ll_split_asymancmig))) + (2*6)
        print "AIC = ", aic_split_asymancmig, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_asymancmig)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()    


    fh_out = open(outname, 'a')
    #####################################
    #Split with ancient symmetrical migration
    #####################################
    print "---------------------------------------------------"
    print "Split with Ancient Symmetrical Migration",'\n','\n'
    fh_out.write("Split with Ancient Symmetrical Migration:"+'\n')

    #first call a predefined model
    func_split_symancmig = TwoD_models.split_ancient_symig

    #create an extrapolating function 
    func_split_symancmig_exec = dadi.Numerics.make_extrap_log_func(func_split_symancmig)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0, 0]
    upper_bound = [30, 30, 20, 10, 10]
    print "parameter set = [nu1, nu2, m, T, Tam]"
    fh_out.write("parameter set = [nu1, nu2, m, T, Tam]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_symancmig replicate {}:".format(i)

        #perturb initial guesses
        param_split_symancmig = dadi.Misc.perturb_params([1,1,0.8,1,1], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_symancmig_opt = dadi.Inference.optimize_log_fmin(param_split_symancmig, fs, func_split_symancmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_symancmig_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_symancmig_opt)+'\n')

        #simulate the model with the optimized parameters
        split_symancmig_model = func_split_symancmig_exec(param_split_symancmig_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_symancmig = dadi.Inference.ll_multinom(split_symancmig_model, fs)
        ll_split_symancmig = np.around(ll_split_symancmig, 3)
        print "likelihood of Split with Ancient Symmetrical Migration model = ", ll_split_symancmig
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_symancmig)+'\n')

        #calculate theta
        theta_split_symancmig = dadi.Inference.optimal_sfs_scaling(split_symancmig_model, fs)
        theta_split_symancmig = np.around(theta_split_symancmig, 3)
        print "Theta = ", theta_split_symancmig
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_symancmig)+'\n')

        #calculate AIC 
        aic_split_symancmig = ( -2*( float(ll_split_symancmig))) + (2*5)
        print "AIC = ", aic_split_symancmig, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_symancmig)+'\n'+'\n')
    print "---------------------------------------------------", '\n'

    fh_out.close()    

    ########################################################################################################
    ########################################################################################################

    fh_out = open(outname, 'a') 
    #####################################
    #Split with no migration, size change
    #####################################
    print "---------------------------------------------------"
    print "Split with No Migration, size change",'\n','\n'
    fh_out.write("Split with No Migration, size change:"+'\n')

    #first call a predefined model
    func_split_nomig_size = TwoD_models.split_no_mig_size

    #create an extrapolating function 
    func_split_nomig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_nomig_size)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
    upper_bound = [30, 30, 30, 30, 10, 10]
    print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2]"
    fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_nomig_size replicate {}:".format(i)

        #perturb initial guesses
        param_split_nomig_size = dadi.Misc.perturb_params([1,1,1,1,0.8,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_nomig_size_opt = dadi.Inference.optimize_log_fmin(param_split_nomig_size, fs, func_split_nomig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_nomig_size_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_nomig_size_opt)+'\n')

        #simulate the model with the optimized parameters
        split_nomig_size_model = func_split_nomig_size_exec(param_split_nomig_size_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_nomig_size = dadi.Inference.ll_multinom(split_nomig_size_model, fs)
        ll_split_nomig_size = np.around(ll_split_nomig_size, 3)
        print "likelihood of Split with No Migration Size Change model = ", ll_split_nomig_size
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_nomig_size)+'\n')

        #calculate theta
        theta_split_nomig_size = dadi.Inference.optimal_sfs_scaling(split_nomig_size_model, fs)
        theta_split_nomig_size = np.around(theta_split_nomig_size, 3)
        print "Theta = ", theta_split_nomig_size
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_nomig_size)+'\n')

        #calculate AIC 
        aic_split_nomig_size = ( -2*( float(ll_split_nomig_size))) + (2*6)
        print "AIC = ", aic_split_nomig_size, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_nomig_size)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()


    fh_out = open(outname, 'a') 
    #####################################
    #Split with secondary contact, asymmetrical gene flow size change
    #####################################
    print "---------------------------------------------------"
    print "Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow",'\n','\n'
    fh_out.write("Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow:"+'\n')

    #first call a predefined model
    func_split_secasymmig_size = TwoD_models.split_secondary_contact_asym_size

    #create an extrapolating function 
    func_split_secasymmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_secasymmig_size)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01, 0.01]
    upper_bound = [30, 30, 30, 30, 10, 10, 20, 20]
    print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"
    fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_secasymmig_size replicate {}:".format(i)
        #perturb initial guesses
        param_split_secasymmig_size = dadi.Misc.perturb_params([1,1,1,1,1,1,0.8,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_secasymmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_secasymmig_size, fs, func_split_secasymmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_secasymmig_size_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_secasymmig_size_opt)+'\n')

        #simulate the model with the optimized parameters
        split_secasymmig_size_model = func_split_secasymmig_size_exec(param_split_secasymmig_size_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secasymmig_size = dadi.Inference.ll_multinom(split_secasymmig_size_model, fs)
        ll_split_secasymmig_size = np.around(ll_split_secasymmig_size, 3)
        print "likelihood of Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow = ", ll_split_secasymmig_size
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_secasymmig_size)+'\n')

        #calculate theta
        theta_split_secasymmig_size = dadi.Inference.optimal_sfs_scaling(split_secasymmig_size_model, fs)
        theta_split_secasymmig_size = np.around(theta_split_secasymmig_size, 3)
        print "Theta = ", theta_split_secasymmig_size
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_secasymmig_size)+'\n')

        #calculate AIC 
        aic_split_secasymmig_size = ( -2*( float(ll_split_secasymmig_size))) + (2*8)
        print "AIC = ", aic_split_secasymmig_size, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_secasymmig_size)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()


    fh_out = open(outname, 'a') 
    #####################################
    #Split with secondary contact, symmetrical gene flow, size change
    #####################################
    print "---------------------------------------------------"
    print "Split with No Mig, Size Change, Secondary Contact, Symmetrical Gene Flow",'\n','\n'
    fh_out.write("Split with No Mig, Size Change, Secondary Contact, Symmetrical Gene Flow:"+'\n')

    #first call a predefined model
    func_split_secsymmig_size = TwoD_models.split_secondary_contact_sym_size

    #create an extrapolating function 
    func_split_secsymmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_secsymmig_size)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01]
    upper_bound = [30, 30, 30, 30, 10, 10, 20]
    print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"
    fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_secsymmig_size replicate {}:".format(i)
        #perturb initial guesses
        param_split_secsymmig_size = dadi.Misc.perturb_params([1,1,1,1,1,1,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_secsymmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_secsymmig_size, fs, func_split_secsymmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_secsymmig_size_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_secsymmig_size_opt)+'\n')

        #simulate the model with the optimized parameters
        split_secsymmig_size_model = func_split_secsymmig_size_exec(param_split_secsymmig_size_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_secsymmig_size = dadi.Inference.ll_multinom(split_secsymmig_size_model, fs)
        ll_split_secsymmig_size = np.around(ll_split_secsymmig_size, 3)
        print "likelihood of Split with No Mig, Size Change, Secondary Contact, Symmetric Migration model = ", ll_split_secsymmig_size
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_secsymmig_size)+'\n')

        #calculate theta
        theta_split_secsymmig_size = dadi.Inference.optimal_sfs_scaling(split_secsymmig_size_model, fs)
        theta_split_secsymmig_size = np.around(theta_split_secsymmig_size, 3)
        print "Theta = ", theta_split_secsymmig_size
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_secsymmig_size)+'\n')

        #calculate AIC 
        aic_split_secsymmig_size = ( -2*( float(ll_split_secsymmig_size))) + (2*7)
        print "AIC = ", aic_split_secsymmig_size, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_secsymmig_size)+'\n'+'\n')
    print "---------------------------------------------------", '\n'
    fh_out.close()


    fh_out = open(outname, 'a') 
    #####################################
    #Split with ancient Asymmetrical migration, size change
    #####################################
    print "---------------------------------------------------"
    print "Split with Ancient Asymmetrical Migration, Size Change",'\n','\n'
    fh_out.write("Split with Ancient Asymmetrical Migration, Size Change:"+'\n')

    #first call a predefined model
    func_split_asymancmig_size = TwoD_models.split_ancient_asymmig_size

    #create an extrapolating function 
    func_split_asymancmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_asymancmig_size)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01, 0.01]
    upper_bound = [30, 30, 30, 30, 10, 10, 20, 20]
    print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"
    fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_asymancmig_size replicate {}:".format(i)

        #perturb initial guesses
        param_split_asymancmig_size = dadi.Misc.perturb_params([1,1,1,1,1,1,0.8,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_asymancmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_asymancmig_size, fs, func_split_asymancmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_asymancmig_size_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_asymancmig_size_opt)+'\n')

        #simulate the model with the optimized parameters
        split_asymancmig_size_model = func_split_asymancmig_size_exec(param_split_asymancmig_size_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_asymancmig_size = dadi.Inference.ll_multinom(split_asymancmig_size_model, fs)
        ll_split_asymancmig_size = np.around(ll_split_asymancmig_size, 3)
        print "likelihood of Split with Ancient Asymmetrical Migration, Size Change = ", ll_split_asymancmig_size
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_asymancmig_size)+'\n')

        #calculate theta
        theta_split_asymancmig_size = dadi.Inference.optimal_sfs_scaling(split_asymancmig_size_model, fs)
        theta_split_asymancmig_size = np.around(theta_split_asymancmig_size, 3)
        print "Theta = ", theta_split_asymancmig_size
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_asymancmig_size)+'\n')

        #calculate AIC 
        aic_split_asymancmig_size = ( -2*( float(ll_split_asymancmig_size))) + (2*8)
        print "AIC = ", aic_split_asymancmig_size, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_asymancmig_size)+'\n'+'\n')
    print "---------------------------------------------------", '\n'

    fh_out.close()    

    fh_out = open(outname, 'a')
    #####################################
    #Split with ancient symmetrical migration, size change
    #####################################
    print "---------------------------------------------------"
    print "Split with Ancient Symmetrical Migration, Size Change",'\n','\n'
    fh_out.write("Split with Ancient Symmetrical Migration, Size Change:"+'\n')

    #first call a predefined model
    func_split_symancmig_size = TwoD_models.split_ancient_symmig_size

    #create an extrapolating function 
    func_split_symancmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_symancmig_size)

    #create parameter list for optimization, set bounds for search
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01]
    upper_bound = [30, 30, 30, 30, 10, 10, 20]
    print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"
    fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"+'\n')

    for i in range(1,x):
        fh_out.write("\tReplicate {}: ".format(i)+'\n')
        print '\n', "split_symancmig_size replicate {}:".format(i)

        #perturb initial guesses
        param_split_symancmig_size = dadi.Misc.perturb_params([1,1,1,1,1,1,0.8], fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

        #run optimization from different starting points, don't print to screen
        param_split_symancmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_symancmig_size, fs, func_split_symancmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
        print '\n',"optimized parameters = ", param_split_symancmig_size_opt
        fh_out.write("\t\tparameters = \t\t\t\t{}".format(param_split_symancmig_size_opt)+'\n')

        #simulate the model with the optimized parameters
        split_symancmig_size_model = func_split_symancmig_size_exec(param_split_symancmig_size_opt, fs.sample_sizes, pts)

        #calculate likelihood
        ll_split_symancmig_size = dadi.Inference.ll_multinom(split_symancmig_size_model, fs)
        ll_split_symancmig_size = np.around(ll_split_symancmig_size, 3)
        print "likelihood of Split with Ancient Symmetrical Migration, Size Change model = ", ll_split_symancmig_size
        fh_out.write("\t\tlog-likelihood = \t\t\t{}".format(ll_split_symancmig_size)+'\n')

        #calculate theta
        theta_split_symancmig_size = dadi.Inference.optimal_sfs_scaling(split_symancmig_size_model, fs)
        theta_split_symancmig_size = np.around(theta_split_symancmig_size, 3)
        print "Theta = ", theta_split_symancmig_size
        fh_out.write("\t\tTheta = \t\t\t\t\t{}".format(theta_split_symancmig_size)+'\n')

        #calculate AIC 
        aic_split_symancmig_size = ( -2*( float(ll_split_symancmig_size))) + (2*7)
        print "AIC = ", aic_split_symancmig_size, '\n', '\n'
        fh_out.write("\t\tAIC = \t\t\t\t\t\t{}".format(aic_split_symancmig_size)+'\n'+'\n')
    print "---------------------------------------------------", '\n'

    fh_out.close()    


#======================================================================================
#Finally, execute function with appropriate arguments
#Two_Pop_Models takes four arguments: (pts, fs, outfile, reps)

# pts = grid choice (list of three numbers, ex. [20,30,40]
# fs = spectrum object name
# outfile = prefix for output naming "{outfile}_INITIAL_model_results.txt"
# reps = integer to control number of replicates, ex. 20

# Ex. Two_Pop_Models([40,50,60], fs_1, "West_East", 10)

#**************
Two_Pop_Models([50,60,70], fs_1, "CVL_South", 20)

#======================================================================================
