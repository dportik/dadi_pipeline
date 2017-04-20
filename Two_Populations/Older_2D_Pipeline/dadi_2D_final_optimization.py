import sys
import os
import numpy as np
import dadi
import pylab
import TwoD_models
import matplotlib.pyplot as plt
'''
usage: python dadi_2D_final_optimization.py

Requires the TwoD_models.py script to be in same working directory. This is where
all the population model functions are stored. This is written for the
models specifically found in that script.

Script will perform final optimizations from multiple starting points using the
perturbed best input parameters selected by the user. The output is written as 
"{outfile}_OPTIMIZED_model_results.txt". You'll want to match projections and grid
sizes as used in the initial optimization script.


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
#create function to run models with user input parameters (presumably from best runs)
#run optimization 'x' times on perturbed best starting params

def Two_Pop_Models(pts, fs, params, modelchoice, outfile, x):
    #output file name
    outname = outfile+"_OPTIMIZED_model_results.txt"

    x = int(x) + int(1)
    #set variable to control maxiter argument
    y = int(60)
        
    if modelchoice == "snm":
        fh_out = open(outname, 'a')
        #####################################
        #Neutral model, no split
        #####################################
        print "---------------------------------------------------"
        print "Neutral Model, No Divergence",'\n','\n'

        #first call a predefined model
        func_snm = TwoD_models.snm

        #create empty parameter list for neutral model (there are no params here)
        no_params = []
        print "parameter set = [none]"
        
        #create an extrapolating function 
        func_snm_exec = dadi.Numerics.make_extrap_log_func(func_snm)
        
        for i in range(1,x):
            fh_out.write("Neutral model, no split"+'\t')
            fh_out.write("parameter set = [none]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "neutral no split replicate {}:".format(i)

            #run the model with the extrapolation call and all our parameters
            snm_model = func_snm_exec(no_params, fs.sample_sizes, pts)

            #calculate likelihood
            ll_neutral = dadi.Inference.ll_multinom(snm_model, fs)
            ll_neutral = np.around(ll_neutral, 2)
            print "likelihood of Neutral Model, No Divergence model = ", ll_neutral
            fh_out.write("{}\t".format(ll_neutral))

            #calculate theta
            theta_snm = dadi.Inference.optimal_sfs_scaling(snm_model, fs)
            theta_snm = np.around(theta_snm, 2)
            print "Theta = ", theta_snm
            fh_out.write("{}\t".format(theta_snm))

            #calculate AIC 
            #AIC = -2 ( ln ( likelihood )) + 2 K
            aic_neutral = ( -2*( float(ll_neutral)))
            print "AIC = ", aic_neutral, '\n', '\n'
            fh_out.write("{}\t".format(aic_neutral)+'\n')
        print "---------------------------------------------------", '\n'
        fh_out.close()


    elif modelchoice == "bottlegrowth":
        fh_out = open(outname, 'a')
        #####################################
        #bottlegrowth model, no split
        #####################################
        print "---------------------------------------------------"
        print "Bottlegrowth Model, No Divergence",'\n','\n'

        #first call a predefined model
        func_bottle = TwoD_models.bottlegrowth

        #create an extrapolating function 
        func_bottle_exec = dadi.Numerics.make_extrap_log_func(func_bottle)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0]
        upper_bound = [30, 30, 10]
        print "parameter set = [nu1, nu2, T]"

        for i in range(1,x):
            fh_out.write("Bottlegrowth Model, No Divergence"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "bottlegrowth no split replicate {}:".format(i)
            
            #perturb best input params
            param_bottle = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_bottle_opt = dadi.Inference.optimize_log_fmin(param_bottle, fs, func_bottle_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_bottle_opt

            #simulate the model with the optimized parameters
            bottle_model = func_bottle_exec(param_bottle_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_bottle = dadi.Inference.ll_multinom(bottle_model, fs)
            ll_bottle = np.around(ll_bottle, 2)
            print "likelihood of Bottlegrowth Model, No Divergence model = ", ll_bottle
            fh_out.write("{}\t".format(ll_bottle))

            #calculate theta
            theta_bottle = dadi.Inference.optimal_sfs_scaling(bottle_model, fs)
            theta_bottle = np.around(theta_bottle, 2)
            print "Theta = ", theta_bottle
            fh_out.write("{}\t".format(theta_bottle))

            #calculate AIC 
            #AIC = -2 ( ln ( likelihood )) + 2 K
            aic_bottle = ( -2*( float(ll_bottle))) + (2*3)
            print "AIC = ", aic_bottle, '\n', '\n'
            fh_out.write("{}\t".format(aic_bottle))
            
            for p in param_bottle_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')


        print "---------------------------------------------------", '\n'
        fh_out.close()
          
    elif modelchoice == "split_nomig":
        fh_out = open(outname, 'a')
        #####################################
        #Split with no migration
        #####################################
        print "---------------------------------------------------"
        print "Split with No Migration",'\n','\n'

        #first call a predefined model
        func_split_nomig = TwoD_models.split_no_mig

        #create an extrapolating function 
        func_split_nomig_exec = dadi.Numerics.make_extrap_log_func(func_split_nomig)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0]
        upper_bound = [30, 30, 10]
        print "parameter set = [nu1, nu2, T]"

        for i in range(1,x):
            fh_out.write("Split with No Migration:"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_nomig replicate {}:".format(i)

            #perturb input params
            param_split_nomig = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_nomig_opt = dadi.Inference.optimize_log_fmin(param_split_nomig, fs, func_split_nomig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_nomig_opt

            #simulate the model with the optimized parameters
            split_nomig_model = func_split_nomig_exec(param_split_nomig_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_nomig = dadi.Inference.ll_multinom(split_nomig_model, fs)
            ll_split_nomig = np.around(ll_split_nomig, 2)
            print "likelihood of Split with No Migration model = ", ll_split_nomig
            fh_out.write("{}\t".format(ll_split_nomig))

            #calculate theta
            theta_split_nomig = dadi.Inference.optimal_sfs_scaling(split_nomig_model, fs)
            theta_split_nomig = np.around(theta_split_nomig, 2)
            print "Theta = ", theta_split_nomig
            fh_out.write("{}\t".format(theta_split_nomig))

            #calculate AIC 
            aic_split_nomig = ( -2*( float(ll_split_nomig))) + (2*3)
            print "AIC = ", aic_split_nomig, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_nomig))
            
            for p in param_split_nomig_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()
    
    elif modelchoice == "split_symmig":
        fh_out = open(outname, 'a')
        #####################################
        #Split with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Symmetric Migration",'\n','\n'

        #first call a predefined model
        func_split_symmig = TwoD_models.split_mig

        #create an extrapolating function 
        func_split_symmig_exec = dadi.Numerics.make_extrap_log_func(func_split_symmig)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0, 0.01]
        upper_bound = [30, 30, 10, 20]
        print "parameter set = [nu1, nu2, T, m]"

        for i in range(1,x):
            fh_out.write("Split with Symmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T, m]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_symmig replicate {}:".format(i)

            ##perturb best input params
            param_split_symmig = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_symmig_opt = dadi.Inference.optimize_log_fmin(param_split_symmig, fs, func_split_symmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_symmig_opt

            #simulate the model with the optimized parameters
            split_symmig_model = func_split_symmig_exec(param_split_symmig_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_symmig = dadi.Inference.ll_multinom(split_symmig_model, fs)
            ll_split_symmig = np.around(ll_split_symmig, 2)
            print "likelihood of Split with Symmetric Migration model = ", ll_split_symmig
            fh_out.write("{}\t".format(ll_split_symmig))

            #calculate theta
            theta_split_symmig = dadi.Inference.optimal_sfs_scaling(split_symmig_model, fs)
            theta_split_symmig = np.around(theta_split_symmig, 2)
            print "Theta = ", theta_split_symmig
            fh_out.write("{}\t".format(theta_split_symmig))

            #calculate AIC 
            aic_split_symmig = ( -2*( float(ll_split_symmig))) + (2*4)
            print "AIC = ", aic_split_symmig, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_symmig))
            
            for p in param_split_symmig_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()

    elif modelchoice == "split_asymmig":
        fh_out = open(outname, 'a')
        #####################################
        #Split with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Asymmetric Migration",'\n','\n'

        #first call a predefined model
        func_split_asymmig = TwoD_models.split_asym_mig

        #create an extrapolating function 
        func_split_asymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_asymmig)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0, 0.01, 0.01]
        upper_bound = [30, 30, 10, 20, 20]
        print "parameter set = [nu1, nu2, T, m12, m21]"

        for i in range(1,x):
            fh_out.write("Split with Asymmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T, m12, m21]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_asymmig replicate {}:".format(i)

            #perturb best input params
            param_split_asymmig = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            param_split_asymmig_opt = dadi.Inference.optimize_log_fmin(param_split_asymmig, fs, func_split_asymmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_asymmig_opt

            #simulate the model with the optimized parameters
            split_asymmig_model = func_split_asymmig_exec(param_split_asymmig_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_asymmig = dadi.Inference.ll_multinom(split_asymmig_model, fs)
            ll_split_asymmig = np.around(ll_split_asymmig, 2)
            print "likelihood of Split with Asymmetric Migration model = ", ll_split_asymmig
            fh_out.write("{}\t".format(ll_split_asymmig))

            #calculate theta
            theta_split_asymmig = dadi.Inference.optimal_sfs_scaling(split_asymmig_model, fs)
            theta_split_asymmig = np.around(theta_split_asymmig, 2)
            print "Theta = ", theta_split_asymmig
            fh_out.write("{}\t".format(theta_split_asymmig))

            #calculate AIC 
            aic_split_asymmig = ( -2*( float(ll_split_asymmig))) + (2*5)
            print "AIC = ", aic_split_asymmig, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_asymmig))
            
            for p in param_split_asymmig_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()

    elif modelchoice == "secasymmig":
        fh_out = open(outname, 'a')
        #####################################
        #Split with secondary contact, asymmetrical gene flow
        #####################################
        print "---------------------------------------------------"
        print "Split with Secondary Contact, Asymmetrical Gene Flow",'\n','\n'

        #first call a predefined model
        func_split_secasymmig = TwoD_models.split_secondary_contact_asym

        #create an extrapolating function 
        func_split_secasymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_secasymmig)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T, Tsc]"
        
        for i in range(1,x):
            fh_out.write("Split with Secondary Contact, Asymmetrical Gene Flow"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T, Tsc]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_secasymmig replicate {}:".format(i)

            #perturb best input params
            param_split_secasymmig = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_secasymmig_opt = dadi.Inference.optimize_log_fmin(param_split_secasymmig, fs, func_split_secasymmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_secasymmig_opt

            #simulate the model with the optimized parameters
            split_secasymmig_model = func_split_secasymmig_exec(param_split_secasymmig_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_secasymmig = dadi.Inference.ll_multinom(split_secasymmig_model, fs)
            ll_split_secasymmig = np.around(ll_split_secasymmig, 2)
            print "likelihood of Split with Secondary Contact and Asymmetric Migration model = ", ll_split_secasymmig
            fh_out.write("{}\t".format(ll_split_secasymmig))

            #calculate theta
            theta_split_secasymmig = dadi.Inference.optimal_sfs_scaling(split_secasymmig_model, fs)
            theta_split_secasymmig = np.around(theta_split_secasymmig, 2)
            print "Theta = ", theta_split_secasymmig
            fh_out.write("{}\t".format(theta_split_secasymmig))

            #calculate AIC 
            aic_split_secasymmig = ( -2*( float(ll_split_secasymmig))) + (2*6)
            print "AIC = ", aic_split_secasymmig, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_secasymmig))
            for p in param_split_secasymmig_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()

    elif modelchoice == "secsymmig":
        fh_out = open(outname, 'a')
        #####################################
        #Split with secondary contact, symmetrical gene flow
        #####################################
        print "---------------------------------------------------"
        print "Split with Secondary Contact, Symmetrical Gene Flow",'\n','\n'

        #first call a predefined model
        func_split_secsymmig = TwoD_models.split_secondary_contact_sym

        #create an extrapolating function 
        func_split_secsymmig_exec = dadi.Numerics.make_extrap_log_func(func_split_secsymmig)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T, Tsc]"

        for i in range(1,x):
            fh_out.write("Split with Secondary Contact, Symmetrical Gene Flow"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m1, T, Tsc]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_secsymmig replicate {}:".format(i)

            #perturb best input params
            param_split_secsymmig = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_secsymmig_opt = dadi.Inference.optimize_log_fmin(param_split_secsymmig, fs, func_split_secsymmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_secsymmig_opt

            #simulate the model with the optimized parameters
            split_secsymmig_model = func_split_secsymmig_exec(param_split_secsymmig_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_secsymmig = dadi.Inference.ll_multinom(split_secsymmig_model, fs)
            ll_split_secsymmig = np.around(ll_split_secsymmig, 2)
            print "likelihood of Split with Secondary Contact and Symmetric Migration model = ", ll_split_secsymmig
            fh_out.write("{}\t".format(ll_split_secsymmig))

            #calculate theta
            theta_split_secsymmig = dadi.Inference.optimal_sfs_scaling(split_secsymmig_model, fs)
            theta_split_secsymmig = np.around(theta_split_secsymmig, 2)
            print "Theta = ", theta_split_secsymmig
            fh_out.write("{}\t".format(theta_split_secsymmig))

            #calculate AIC 
            aic_split_secsymmig = ( -2*( float(ll_split_secsymmig))) + (2*5)
            print "AIC = ", aic_split_secsymmig, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_secsymmig))
            for p in param_split_secsymmig_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()

    elif modelchoice == "ancasymmig":
        fh_out = open(outname, 'a')
        #####################################
        #Split with ancient Asymmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Asymmetrical Migration",'\n','\n'

        #first call a predefined model
        func_split_asymancmig = TwoD_models.split_ancient_asymig

        #create an extrapolating function 
        func_split_asymancmig_exec = dadi.Numerics.make_extrap_log_func(func_split_asymancmig)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T, Tam]"

        for i in range(1,x):
            fh_out.write("Split with Ancient Asymmetrical Migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T, Tam]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_asymancmig replicate {}:".format(i)

            #perturb best input params
            param_split_asymancmig = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_asymancmig_opt = dadi.Inference.optimize_log_fmin(param_split_asymancmig, fs, func_split_asymancmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_asymancmig_opt

            #simulate the model with the optimized parameters
            split_asymancmig_model = func_split_asymancmig_exec(param_split_asymancmig_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_asymancmig = dadi.Inference.ll_multinom(split_asymancmig_model, fs)
            ll_split_asymancmig = np.around(ll_split_asymancmig, 2)
            print "likelihood of Split with Ancient Symmetrical Migration model = ", ll_split_asymancmig
            fh_out.write("{}\t".format(ll_split_asymancmig))

            #calculate theta
            theta_split_asymancmig = dadi.Inference.optimal_sfs_scaling(split_asymancmig_model, fs)
            theta_split_asymancmig = np.around(theta_split_asymancmig, 2)
            print "Theta = ", theta_split_asymancmig
            fh_out.write("{}\t".format(theta_split_asymancmig))

            #calculate AIC 
            aic_split_asymancmig = ( -2*( float(ll_split_asymancmig))) + (2*6)
            print "AIC = ", aic_split_asymancmig, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_asymancmig))
            
            for p in param_split_asymancmig_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()    

    elif modelchoice == "ancsymmig":
        fh_out = open(outname, 'a') 
        #####################################
        #Split with ancient symmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Symmetrical Migration",'\n','\n'

        #first call a predefined model
        func_split_symancmig = TwoD_models.split_ancient_symig

        #create an extrapolating function 
        func_split_symancmig_exec = dadi.Numerics.make_extrap_log_func(func_split_symancmig)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T, Tam]"

        for i in range(1,x):
            fh_out.write("Split with Ancient Symmetrical Migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T, Tam]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_symancmig replicate {}:".format(i)

            #perturb best input params
            param_split_symancmig = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            param_split_symancmig_opt = dadi.Inference.optimize_log_fmin(param_split_symancmig, fs, func_split_symancmig_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_symancmig_opt

            #simulate the model with the optimized parameters
            split_symancmig_model = func_split_symancmig_exec(param_split_symancmig_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_symancmig = dadi.Inference.ll_multinom(split_symancmig_model, fs)
            ll_split_symancmig = np.around(ll_split_symancmig, 2)
            print "likelihood of Split with Ancient Symmetrical Migration model = ", ll_split_symancmig
            fh_out.write("{}\t".format(ll_split_symancmig))

            #calculate theta
            theta_split_symancmig = dadi.Inference.optimal_sfs_scaling(split_symancmig_model, fs)
            theta_split_symancmig = np.around(theta_split_symancmig, 2)
            print "Theta = ", theta_split_symancmig
            fh_out.write("{}\t".format(theta_split_symancmig))

            #calculate AIC 
            aic_split_symancmig = ( -2*( float(ll_split_symancmig))) + (2*5)
            print "AIC = ", aic_split_symancmig, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_symancmig))
            
            for p in param_split_symancmig_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            
        print "---------------------------------------------------", '\n'
        fh_out.close()    

    elif modelchoice == "split_no_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Split with no migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Split with No Migration, size change",'\n','\n'

        #first call a predefined model
        func_split_nomig_size = TwoD_models.split_no_mig_size

        #create an extrapolating function 
        func_split_nomig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_nomig_size)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 10, 10]
        print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2]"
 
        for i in range(1,x):
            fh_out.write("Split with No Migration, size change"+'\t')
            fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_nomig_size replicate {}:".format(i)

            #perturb best input params
            param_split_nomig_size = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            param_split_nomig_size_opt = dadi.Inference.optimize_log_fmin(param_split_nomig_size, fs, func_split_nomig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_nomig_size_opt

            #simulate the model with the optimized parameters
            split_nomig_size_model = func_split_nomig_size_exec(param_split_nomig_size_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_nomig_size = dadi.Inference.ll_multinom(split_nomig_size_model, fs)
            ll_split_nomig_size = np.around(ll_split_nomig_size, 2)
            print "likelihood of Split with No Migration Size Change model = ", ll_split_nomig_size
            fh_out.write("{}\t".format(ll_split_nomig_size))

            #calculate theta
            theta_split_nomig_size = dadi.Inference.optimal_sfs_scaling(split_nomig_size_model, fs)
            theta_split_nomig_size = np.around(theta_split_nomig_size, 2)
            print "Theta = ", theta_split_nomig_size
            fh_out.write("{}\t".format(theta_split_nomig_size))

            #calculate AIC 
            aic_split_nomig_size = ( -2*( float(ll_split_nomig_size))) + (2*6)
            print "AIC = ", aic_split_nomig_size, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_nomig_size))
            
            for p in param_split_nomig_size_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()

    elif modelchoice == "secasymmig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Split with secondary contact, asymmetrical gene flow size change
        #####################################
        print "---------------------------------------------------"
        print "Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow",'\n','\n'

        #first call a predefined model
        func_split_secasymmig_size = TwoD_models.split_secondary_contact_asym_size

        #create an extrapolating function 
        func_split_secasymmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_secasymmig_size)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01, 0.01]
        upper_bound = [30, 30, 30, 30, 10, 10, 20, 20]
        print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"
    
        for i in range(1,x):
            fh_out.write("Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow"+'\t')
            fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_secasymmig_size replicate {}:".format(i)

            #perturb best input params
            param_split_secasymmig_size = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_secasymmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_secasymmig_size, fs, func_split_secasymmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_secasymmig_size_opt

            #simulate the model with the optimized parameters
            split_secasymmig_size_model = func_split_secasymmig_size_exec(param_split_secasymmig_size_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_secasymmig_size = dadi.Inference.ll_multinom(split_secasymmig_size_model, fs)
            ll_split_secasymmig_size = np.around(ll_split_secasymmig_size, 2)
            print "likelihood of Split with No Mig, Size Change, Secondary Contact, Asymmetrical Gene Flow = ", ll_split_secasymmig_size
            fh_out.write("{}\t".format(ll_split_secasymmig_size))

            #calculate theta
            theta_split_secasymmig_size = dadi.Inference.optimal_sfs_scaling(split_secasymmig_size_model, fs)
            theta_split_secasymmig_size = np.around(theta_split_secasymmig_size, 2)
            print "Theta = ", theta_split_secasymmig_size
            fh_out.write("{}\t".format(theta_split_secasymmig_size))

            #calculate AIC 
            aic_split_secasymmig_size = ( -2*( float(ll_split_secasymmig_size))) + (2*8)
            print "AIC = ", aic_split_secasymmig_size, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_secasymmig_size))
            
            for p in param_split_secasymmig_size_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()

    elif modelchoice == "secsymmig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Split with secondary contact, symmetrical gene flow size change
        #####################################
        print "---------------------------------------------------"
        print "Split with No Mig, Size Change, Secondary Contact, Symmetrical Gene Flow",'\n','\n'

        #first call a predefined model
        func_split_secsymmig_size = TwoD_models.split_secondary_contact_sym_size

        #create an extrapolating function 
        func_split_secsymmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_secsymmig_size)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01]
        upper_bound = [30, 30, 30, 30, 10, 10, 20]
        print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"

        for i in range(1,x):
            fh_out.write("Split with No Mig, Size Change, Secondary Contact, Symmetrical Gene Flow"+'\t')
            fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_secsymmig_size replicate {}:".format(i)

            #perturb best input params
            param_split_secsymmig_size = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_secsymmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_secsymmig_size, fs, func_split_secsymmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_secsymmig_size_opt

            #simulate the model with the optimized parameters
            split_secsymmig_size_model = func_split_secsymmig_size_exec(param_split_secsymmig_size_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_secsymmig_size = dadi.Inference.ll_multinom(split_secsymmig_size_model, fs)
            ll_split_secsymmig_size = np.around(ll_split_secsymmig_size, 2)
            print "likelihood of Split with No Mig, Size Change, Secondary Contact, Symmetric Migration model = ", ll_split_secsymmig_size
            fh_out.write("{}\t".format(ll_split_secsymmig_size))

            #calculate theta
            theta_split_secsymmig_size = dadi.Inference.optimal_sfs_scaling(split_secsymmig_size_model, fs)
            theta_split_secsymmig_size = np.around(theta_split_secsymmig_size, 2)
            print "Theta = ", theta_split_secsymmig_size
            fh_out.write("{}\t".format(theta_split_secsymmig_size))

            #calculate AIC 
            aic_split_secsymmig_size = ( -2*( float(ll_split_secsymmig_size))) + (2*7)
            print "AIC = ", aic_split_secsymmig_size, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_secsymmig_size))
            
            for p in param_split_secsymmig_size_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()
    
    elif modelchoice == "ancasymmig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Split with ancient Asymmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Asymmetrical Migration, Size Change",'\n','\n'

        #first call a predefined model
        func_split_asymancmig_size = TwoD_models.split_ancient_asymmig_size

        #create an extrapolating function 
        func_split_asymancmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_asymancmig_size)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01, 0.01]
        upper_bound = [30, 30, 30, 30, 10, 10, 20, 20]
        print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"

        for i in range(1,x):
            fh_out.write("Split with Ancient Asymmetrical Migration, Size Change"+'\t')
            fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m12, m21]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_asymancmig_size replicate {}:".format(i)

            #perturb best input params
            param_split_asymancmig_size = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_asymancmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_asymancmig_size, fs, func_split_asymancmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_asymancmig_size_opt

            #simulate the model with the optimized parameters
            split_asymancmig_size_model = func_split_asymancmig_size_exec(param_split_asymancmig_size_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_asymancmig_size = dadi.Inference.ll_multinom(split_asymancmig_size_model, fs)
            ll_split_asymancmig_size = np.around(ll_split_asymancmig_size, 2)
            print "likelihood of Split with Ancient Asymmetrical Migration, Size Change = ", ll_split_asymancmig_size
            fh_out.write("{}\t".format(ll_split_asymancmig_size))

            #calculate theta
            theta_split_asymancmig_size = dadi.Inference.optimal_sfs_scaling(split_asymancmig_size_model, fs)
            theta_split_asymancmig_size = np.around(theta_split_asymancmig_size, 2)
            print "Theta = ", theta_split_asymancmig_size
            fh_out.write("{}\t".format(theta_split_asymancmig_size))

            #calculate AIC 
            aic_split_asymancmig_size = ( -2*( float(ll_split_asymancmig_size))) + (2*8)
            print "AIC = ", aic_split_asymancmig_size, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_asymancmig_size))
            
            for p in param_split_asymancmig_size_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()    

            
    elif modelchoice == "ancsymmig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Split with ancient symmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Split with Ancient Symmetrical Migration, Size Change",'\n','\n'

        #first call a predefined model
        func_split_symancmig_size = TwoD_models.split_ancient_symmig_size

        #create an extrapolating function 
        func_split_symancmig_size_exec = dadi.Numerics.make_extrap_log_func(func_split_symancmig_size)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0.01]
        upper_bound = [30, 30, 30, 30, 10, 10, 20]
        print "parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"

        for i in range(1,x):
            fh_out.write("Split with Ancient Symmetrical Migration, Size Change"+'\t')
            fh_out.write("parameter set = [nu1b, nu2b, nu1r, nu2r, T1, T2, m]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "split_symancmig_size replicate {}:".format(i)

            #perturb best input params
            param_split_symancmig_size = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            param_split_symancmig_size_opt = dadi.Inference.optimize_log_fmin(param_split_symancmig_size, fs, func_split_symancmig_size_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"starting parameters = ", params
            print "optimized parameters = ", param_split_symancmig_size_opt

            #simulate the model with the optimized parameters
            split_symancmig_size_model = func_split_symancmig_size_exec(param_split_symancmig_size_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll_split_symancmig_size = dadi.Inference.ll_multinom(split_symancmig_size_model, fs)
            ll_split_symancmig_size = np.around(ll_split_symancmig_size, 2)
            print "likelihood of Split with Ancient Symmetrical Migration, Size Change model = ", ll_split_symancmig_size
            fh_out.write("{}\t".format(ll_split_symancmig_size))

            #calculate theta
            theta_split_symancmig_size = dadi.Inference.optimal_sfs_scaling(split_symancmig_size_model, fs)
            theta_split_symancmig_size = np.around(theta_split_symancmig_size, 2)
            print "Theta = ", theta_split_symancmig_size
            fh_out.write("{}\t".format(theta_split_symancmig_size))

            #calculate AIC 
            aic_split_symancmig_size = ( -2*( float(ll_split_symancmig_size))) + (2*7)
            print "AIC = ", aic_split_symancmig_size, '\n', '\n'
            fh_out.write("{}\t".format(aic_split_symancmig_size))
            
            for p in param_split_symancmig_size_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')

        print "---------------------------------------------------", '\n'
        fh_out.close()    
        

#======================================================================================
# Finally, execute model with appropriate arguments
# Two_Pop_Models(pts, fs, params, modelchoice, outfile, reps):

# pts = grid choice (list of three numbers, ex. [20,30,40]
# fs = spectrum object name
# params = list of best parameters to start optimizations from
# modelchoice = snm, bottlegrowth, split_nomig, split_symmig, split_asymmig, secasymmig, secsymmig, ancasymmig,
# 				ancsymmig, split_no_mig_size, secasymmig_size, secsymmig_size, ancasymmig_size, ancsymmig_size
# outfile = prefix for output naming
# reps = integer to control number of replicates, ex. 10


#**************
#Input some of the basic arguments here
pts = [50,60,70]
outfile = "CVL_South"
reps = 20

outname = outfile+"_OPTIMIZED_model_results.txt"
fh_out = open(outname, 'a')
fh_out.write("Model"+'\t'+"param_set"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"AIC"+'\t'+"optimized_params"+'\n')
fh_out.close()


#If you used the initial optimization script, the parameters are in the same order
#so you can copy and paste results from "{outfile}_INITIAL_model_results.txt"
#to here, making sure to put in list format

#Below, you really just need to change the params values for each model

#**************
#neutral, no divergence: no parameters for this model
params1 = []
best_model1 = Two_Pop_Models(pts, fs_1, params1, "snm", outfile, reps)

#**************
#bottlegrowth: params list must contain 3 values
params2 = [0.01000154,0.73999852,0.13444088]
best_model2 = Two_Pop_Models(pts, fs_1, params2, "bottlegrowth", outfile, reps)

#**************
#split_nomig: params list must contain 3 values
params3 = [3.92796509,1.61440453,1.72422307]
best_model3 = Two_Pop_Models(pts, fs_1, params3, "split_nomig", outfile, reps)

#**************
#split_symmig: params list must contain 4 values
params4 = [8.30892475, 3.63535352, 5.53927482, 0.01171039]
best_model4 = Two_Pop_Models(pts, fs_1, params4, "split_symmig", outfile, reps)

#**************
#split_asymmig: params list must contain 5 values
params5 = [3.27477801,1.36630643,2.08793036,0.03903948,0.05170418]
best_model5 = Two_Pop_Models(pts, fs_1, params5, "split_asymmig", outfile, reps)

#**************
#secasymmig: params list must contain 6 values
params6 = [3.22736357,1.25843554,0.23032225,1.63451371,1.48353851,0.00717123]
best_model6 = Two_Pop_Models(pts, fs_1, params6, "secasymmig", outfile, reps)

#**************
#secsymmig: params list must contain 5 values
params7 = [6.62221727,2.90739022,0.0299099,3.7105653,0.36040096]
best_model7 = Two_Pop_Models(pts, fs_1, params7, "secsymmig", outfile, reps)

#**************
#ancasymmig: params list must contain 6 values
params8 = [3.95371971,1.6589575,0.03451684,10.70408597,0.10203438,1.70989959]
best_model8 = Two_Pop_Models(pts, fs_1, params8, "ancasymmig", outfile, reps)

#**************
#ancsymmig: params list must contain 5 values
params9 = [8.99944794,4.05098536,0.01008406,6.28042677,0.02785271]
best_model9 = Two_Pop_Models(pts, fs_1, params9, "ancsymmig", outfile, reps)

#**************
#split_no_mig_size: params list must contain 6 values
params10 = [0.46105388,0.47056775,3.78052092,1.7902496,0.3765935,0.44682741]
best_model10 = Two_Pop_Models(pts, fs_1, params10, "split_no_mig_size", outfile, reps)

#**************
#secasymmig_size: params list must contain 8 values
params11 = [2.26398181,1.08482008,8.8516358,2.32630116,1.42718774,0.08638007,0.07176869,0.24119363]
best_model11 = Two_Pop_Models(pts, fs_1, params11, "secasymmig_size", outfile, reps)

#**************
#secsymmig_size: params list must contain 7 values
params12 = [0.39795623,1.13919786,4.40443606,1.41800793,0.59967262,0.71371194, 0.02354581]
best_model12 = Two_Pop_Models(pts, fs_1, params12, "secsymmig_size", outfile, reps)

#**************
#ancasymmig_size: params list must contain 8 values
params13 = [0.44345171,1.13116023,7.26993936,2.42418427,3.16343886,1.23189833,0.37458392,0.58187062]
best_model13 = Two_Pop_Models(pts, fs_1, params13, "ancasymmig_size", outfile, reps)

#**************
#ancsymmig_size: params list must contain 7 values
params14 = [0.94242724,2.58511376,10.87045377,3.86404326,5.05774359,1.99828185,0.26165181]
best_model14 = Two_Pop_Models(pts, fs_1, params14, "ancsymmig_size", outfile, reps)
