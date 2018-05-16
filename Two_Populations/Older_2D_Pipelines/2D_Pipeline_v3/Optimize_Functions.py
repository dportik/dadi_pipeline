import dadi
import numpy as np
import pylab
import Models_2D
import matplotlib.pyplot as plt

#======================================================================================
#======================================================================================
#======================================================================================
#======================================================================================
#Optimizations Round 1 function (run models with basic starting params)

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
#        "sec_contact_asym_mig_size"

def Optimize_Round1(pts, fs, outfile, reps, y, model_name):
    print '\n',"============================================================================"
    print "Beginning analysis of {}".format(model_name)
    print "============================================================================"

    #create output file
    outname = "Round1_{0}_{1}_optimized.txt".format(outfile,model_name)
    fh_out = open(outname, 'a')
    fh_out.write("Model"+'\t'+"param_set"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"AIC"+'\t'+"optimized_params"+'\n')
    fh_out.close()
    
    #variable to control number of loops per model (1 to x)
    x = int(reps) + int(1)
        
    if model_name == "no_divergence":
        #####################################
        #No Divergence
        #####################################
        print "---------------------------------------------------"
        print "No Divergence",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_divergence
        params = []
        print "parameter set = [none]"
        
        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("No divergence model"+'\t')
            fh_out.write("parameter set = [none]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #simulate the model with the optimized parameters
            sim_model = func_exec(params, fs.sample_sizes, pts)
            
            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))
            
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'

        
    elif model_name == "no_mig":
        fh_out = open(outname, 'a')
        #####################################
        # Divergence with no migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0]
        upper_bound = [30, 30, 10]
        params = [1,1,1]
        print "parameter set = [nu1, nu2, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 10]
        params = [1,1,1,1]
        print "parameter set = [nu1, nu2, m, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*4)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 20, 10]
        params = [1,1,1,1,1]
        print "parameter set = [nu1, nu2, m12, m21, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_sym_mig":
        fh_out = open(outname, 'a') 
        #####################################
        #Divergence with ancient symmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        params = [1,1,1,1,1]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            
            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        params = [1,1,1,1,1,1]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sec_contact_sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        params = [1,1,1,1,1]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sec_contact_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        params = [1,1,1,1,1,1]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    ##########################################################################
    #Similar models but with size change allowed
    ##########################################################################
    
    elif model_name == "no_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with no migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 10, 10]
        params = [1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        params = [1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        

    elif model_name == "asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        params = [1,1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
 
 
    elif model_name == "anc_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient symmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        params = [1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        params = [1,1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    elif model_name == "sec_contact_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        params = [1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

           
    elif model_name == "sec_contact_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        params = [1,1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
    
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    ########################################################################################################################################
    #Newer Models added here below
    ########################################################################################################################################

    elif model_name == "sym_mig_twoepoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_twoepoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        params = [1,1,1,1,1,1]
        print "parameter set = [nu1, nu2, m1, m2, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration, two epochs"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m1, m2, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "asym_mig_twoepoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_twoepoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 20, 20, 10, 10]
        params = [1,1,1,1,1,1,1,1]
        print "parameter set = [nu1, nu2, m12a, m21a, m12b, m21b, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration, two epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12a, m21a, m12b, m21b, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    #####################################################################
        
    elif model_name == "sec_contact_sym_mig_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 20, 10, 10, 10]
        params = [1,1,1,1,1,1]
        print "parameter set = [nu1, nu2, m, T1, T2, T3]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, three epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sec_contact_asym_mig_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10, 10]
        params = [1,1,1,1,1,1,1]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2, T3]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, three epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    elif model_name == "sec_contact_sym_mig_size_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10, 10]
        params = [1,1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, size change, three epoch"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

           
    elif model_name == "sec_contact_asym_mig_size_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10, 10]
        params = [1,1,1,1,1,1,1,1,1]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3]"
    
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, size change, three epoch"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*9)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    #####################################################################
    #ISLAND MODELS
    #####################################################################

    elif model_name == "vic_no_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and no migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99]
        params = [1,1,1,1,0.25]
        print "parameter set = [nuA, nu1, nu2, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and no migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
    
    elif model_name == "vic_anc_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and ancient asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and ancient asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 10, 10, 0.99]
        params = [1,1,1,1,1,1,1,0.25]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and ancient asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "vic_sec_contact_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and secondary contact with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and secondary contact with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 10, 10, 0.99]
        params = [1,1,1,1,1,1,1,0.25]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and secondary contact with asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_sym":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_sym

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 20, 10, 0.99]
        params = [1,1,1,1,1,0.25]
        print "parameter set = [nuA, nu1, nu2, m, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and symmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    elif model_name == "founder_asym":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_asym

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [30, 30, 30, 20, 20, 10, 0.99]
        params = [1,1,1,1,1,1,0.25]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and no migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [30, 30, 30, 10, 0.99]
        params = [1,1,1,1,0.25]
        print "parameter set = [nuA, nu1, nu2, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and no migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    ###########################################################################
    #ISLAND Discrete admixture models
    ###########################################################################

    
    elif model_name == "vic_no_mig_admix_early":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance with early discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with early discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_early

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99, 0.99]
        params = [1,1,1,1,0.25,0.25]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with early discrete admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
    
    elif model_name == "vic_no_mig_admix_late":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with late discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with late discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_late

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99, 0.99]
        params = [1,1,1,1,0.25,0.25]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with late discrete admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "vic_two_epoch_admix":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance with discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_two_epoch_admix

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 0.99, 0.99]
        params = [1,1,1,1,1,0.25,0.25]
        print "parameter set = [nuA, nu1, nu2, T1, T2, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with discrete admixture, two epoch"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T1, T2, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    elif model_name == "founder_nomig_admix_early":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete early admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete early admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_early

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 0.99, 0.99]
        params = [1,1,1,1,0.25,0.25]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete early admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig_admix_late":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete late admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete early admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_late

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 0.99, 0.99]
        params = [1,1,1,1,0.25,0.25]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete late admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig_admix_two_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_two_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 10, 0.99, 0.99]
        params = [1,1,1,1,1,0.25,0.25]
        print "parameter set = [nuA, nu1, nu2, T1, T2, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete admixture, two epoch"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T1, T2, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'



                
#======================================================================================
#======================================================================================
#======================================================================================
#======================================================================================
#Optimizations Round 2 function (run models with user input parameters, presumably from best runs)

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
#        "sec_contact_asym_mig_size"
# params:  list of best parameter values to perturb to start the optimizations from

def Optimize_Round2(pts, fs, outfile, reps, y, model_name, params):
    print '\n',"============================================================================"
    print "Beginning analysis of {}".format(model_name)
    print "============================================================================"

    #create output file
    outname = "Round2_{0}_{1}_optimized.txt".format(outfile,model_name)
    fh_out = open(outname, 'a')
    fh_out.write("Model"+'\t'+"param_set"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"AIC"+'\t'+"optimized_params"+'\n')
    fh_out.close()
    
    #variable to control number of loops per model (1 to x)
    x = int(reps) + int(1)
        
    if model_name == "no_divergence":
        #####################################
        #No Divergence
        #####################################
        print "---------------------------------------------------"
        print "No Divergence",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_divergence
        params = []
        print "parameter set = [none]"
        
        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("No divergence model"+'\t')
            fh_out.write("parameter set = [none]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #simulate the model with the optimized parameters
            sim_model = func_exec(params, fs.sample_sizes, pts)
            
            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))
            
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'

        
    elif model_name == "no_mig":
        fh_out = open(outname, 'a')
        #####################################
        # Divergence with no migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0]
        upper_bound = [30, 30, 10]
        print "parameter set = [nu1, nu2, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 10]
        print "parameter set = [nu1, nu2, m, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*4)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 20, 10]
        print "parameter set = [nu1, nu2, m12, m21, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_sym_mig":
        fh_out = open(outname, 'a') 
        #####################################
        #Divergence with ancient symmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sec_contact_sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sec_contact_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    ##########################################################################
    #Similar models but with size change allowed
    ##########################################################################
    
    elif model_name == "no_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with no migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        

    elif model_name == "asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
 
 
    elif model_name == "anc_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient symmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    elif model_name == "sec_contact_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

           
    elif model_name == "sec_contact_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
    
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    ########################################################################################################################################
    #Newer Models added here below
    ########################################################################################################################################

    elif model_name == "sym_mig_twoepoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_twoepoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m1, m2, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration, two epochs"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m1, m2, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "asym_mig_twoepoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_twoepoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12a, m21a, m12b, m21b, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration, two epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12a, m21a, m12b, m21b, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    #####################################################################
        
    elif model_name == "sec_contact_sym_mig_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 20, 10, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2, T3]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, three epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sec_contact_asym_mig_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2, T3]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, three epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    elif model_name == "sec_contact_sym_mig_size_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, size change, three epoch"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

           
    elif model_name == "sec_contact_asym_mig_size_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3]"
    
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, size change, three epoch"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*9)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    #####################################################################
    #ISLAND MODELS
    #####################################################################

    elif model_name == "vic_no_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and no migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and no migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
    
    elif model_name == "vic_anc_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and ancient asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and ancient asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 10, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and ancient asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "vic_sec_contact_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and secondary contact with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and secondary contact with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 10, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and secondary contact with asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_sym":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_sym

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 20, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and symmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    elif model_name == "founder_asym":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_asym

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [30, 30, 30, 20, 20, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and no migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [30, 30, 30, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and no migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    ###########################################################################
    #ISLAND Discrete admixture models
    ###########################################################################

    
    elif model_name == "vic_no_mig_admix_early":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance with early discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with early discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_early

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with early discrete admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
    
    elif model_name == "vic_no_mig_admix_late":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with late discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with late discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_late

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with late discrete admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "vic_two_epoch_admix":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance with discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_two_epoch_admix

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T1, T2, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with discrete admixture, two epoch"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T1, T2, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    elif model_name == "founder_nomig_admix_early":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete early admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete early admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_early

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete early admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig_admix_late":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete late admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete late admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_late

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete late admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig_admix_two_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_two_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T1, T2, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete admixture, two epoch"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T1, T2, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


        
#======================================================================================
#======================================================================================
#======================================================================================
#======================================================================================
#======================================================================================
#Optimizations Round 3 function (run models with user input parameters, presumably from best runs)

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
#        "sec_contact_asym_mig_size"
# params:  list of best parameter values to perturb to start the optimizations from

def Optimize_Round3(pts, fs, outfile, reps, y, model_name, params):
    print '\n',"============================================================================"
    print "Beginning analysis of {}".format(model_name)
    print "============================================================================"

    #create output file
    outname = "Round3_{0}_{1}_optimized.txt".format(outfile,model_name)
    fh_out = open(outname, 'a')
    fh_out.write("Model"+'\t'+"param_set"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"AIC"+'\t'+"optimized_params"+'\n')
    fh_out.close()
    
    #variable to control number of loops per model (1 to x)
    x = int(reps) + int(1)
        
    if model_name == "no_divergence":
        #####################################
        #No Divergence
        #####################################
        print "---------------------------------------------------"
        print "No Divergence",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_divergence
        params = []
        print "parameter set = [none]"
        
        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("No divergence model"+'\t')
            fh_out.write("parameter set = [none]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #simulate the model with the optimized parameters
            sim_model = func_exec(params, fs.sample_sizes, pts)
            
            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))
            
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'

        
    elif model_name == "no_mig":
        fh_out = open(outname, 'a')
        #####################################
        # Divergence with no migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0]
        upper_bound = [30, 30, 10]
        print "parameter set = [nu1, nu2, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 10]
        print "parameter set = [nu1, nu2, m, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*4)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 20, 10]
        print "parameter set = [nu1, nu2, m12, m21, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_sym_mig":
        fh_out = open(outname, 'a') 
        #####################################
        #Divergence with ancient symmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sec_contact_sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sec_contact_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    ##########################################################################
    #Similar models but with size change allowed
    ##########################################################################
    
    elif model_name == "no_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with no migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        

    elif model_name == "asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
 
 
    elif model_name == "anc_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient symmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    elif model_name == "sec_contact_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

           
    elif model_name == "sec_contact_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
    
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    ########################################################################################################################################
    #Newer Models added here below
    ########################################################################################################################################

    elif model_name == "sym_mig_twoepoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_twoepoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m1, m2, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration, two epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m1, m2, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "asym_mig_twoepoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_twoepoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12a, m21a, m12b, m21b, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration, two epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12a, m21a, m12b, m21b, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    #####################################################################
        
    elif model_name == "sec_contact_sym_mig_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 20, 10, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2, T3]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, three epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sec_contact_asym_mig_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2, T3]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, three epoch"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    elif model_name == "sec_contact_sym_mig_size_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, size change, three epoch"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

           
    elif model_name == "sec_contact_asym_mig_size_three_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size_three_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3]"
    
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, size change, three epoch"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*9)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    #####################################################################
    #ISLAND MODELS
    #####################################################################

    elif model_name == "vic_no_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and no migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and no migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
    
    elif model_name == "vic_anc_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and ancient asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and ancient asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 10, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and ancient asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "vic_sec_contact_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance and secondary contact with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and secondary contact with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 10, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance and secondary contact with asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T1, T2, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_sym":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_sym

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [20, 20, 20, 20, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and symmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    elif model_name == "founder_asym":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_asym

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [30, 30, 30, 20, 20, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, m12, m21, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and asymmetric migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, m12, m21, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and no migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001]
        upper_bound = [30, 30, 30, 10, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and no migration"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    ###########################################################################
    #ISLAND Discrete admixture models
    ###########################################################################

    
    elif model_name == "vic_no_mig_admix_early":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance with early discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with early discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_early

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with early discrete admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
    
    elif model_name == "vic_no_mig_admix_late":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with late discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with late discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_late

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with late discrete admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "vic_two_epoch_admix":
        fh_out = open(outname, 'a')
        #####################################
        #Vicariance with discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_two_epoch_admix

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 20, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T1, T2, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Vicariance with discrete admixture, two epoch"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T1, T2, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    elif model_name == "founder_nomig_admix_early":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete early admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete early admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_early

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete early admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig_admix_late":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete late admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete late admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_late

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete late admixture"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

    elif model_name == "founder_nomig_admix_two_epoch":
        fh_out = open(outname, 'a')
        #####################################
        #Founder event and discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_two_epoch

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        upper_bound = [20, 20, 10, 10, 10, 0.99, 0.99]
        print "parameter set = [nuA, nu1, nu2, T1, T2, s, f]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Founder event and discrete admixture, two epoch"+'\t')
            fh_out.write("parameter set = [nuA, nu1, nu2, T1, T2, s, f]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'



        
#======================================================================================
#======================================================================================
#======================================================================================
#======================================================================================        
#Run models for plotting function (only runs once per model with exact input params)

# Function usage:
# returned_model_name = Optimize_Single(pts, fs, model_name, params)

# Argument definitions:
# pts:  grid choice (list of three numbers, ex. [20,30,40]
# fs:  spectrum object name
# model_name:  "no_divergence", "no_mig", "sym_mig", "asym_mig", "anc_sym_mig", "anc_asym_mig",
#        "sec_contact_sym_mig", "sec_contact_asym_mig", "no_mig_size", "sym_mig_size",
#        "asym_mig_size", "anc_sym_mig_size", "anc_asym_mig_size", "sec_contact_sym_mig_size",
#        "sec_contact_asym_mig_size"
# params:  list of best parameter values to optimize with

def Optimize_Single(pts, fs, model_name, params):

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
    
    elif model_name == "sym_mig_twoepoch":
        #####################################
        #Divergence with symmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_twoepoch

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
        
    elif model_name == "asym_mig_twoepoch":
        #####################################
        #Divergence with asymmetric migration, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_twoepoch

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
        
    #####################################################################
        
    elif model_name == "sec_contact_sym_mig_three_epoch":
        #####################################
        #Divergence and symmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_three_epoch

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
        
        
    elif model_name == "sec_contact_asym_mig_three_epoch":
        #####################################
        #Divergence and asymmetrical secondary contact, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_three_epoch

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
        
    elif model_name == "sec_contact_sym_mig_size_three_epoch":
        #####################################
        #Divergence and symmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size_three_epoch

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

           
    elif model_name == "sec_contact_asym_mig_size_three_epoch":
        #####################################
        #Divergence and asymmetrical secondary contact, size change, three epoch
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change, three epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size_three_epoch

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
        aic = ( -2*( float(ll))) + (2*9)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

    #####################################################################
    #ISLAND MODELS
    #####################################################################

    elif model_name == "vic_no_mig":
        #####################################
        #Vicariance and no migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig

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
        aic = ( -2*( float(ll))) + (2*5)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

    
    elif model_name == "vic_anc_asym_mig":
        #####################################
        #Vicariance and ancient asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and ancient asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_anc_asym_mig

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
        aic = ( -2*( float(ll))) + (2*8)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


    elif model_name == "vic_sec_contact_asym_mig":
        #####################################
        #Vicariance and secondary contact with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Vicariance and secondary contact with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_sec_contact_asym_mig

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
        aic = ( -2*( float(ll))) + (2*8)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


    elif model_name == "founder_sym":
        #####################################
        #Founder event and symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_sym

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
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

        
    elif model_name == "founder_asym":
        #####################################
        #Founder event and asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_asym

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
        aic = ( -2*( float(ll))) + (2*7)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


    elif model_name == "founder_nomig":
        #####################################
        #Founder event and no migration
        #####################################
        print "---------------------------------------------------"
        print "Founder event and no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig

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
        aic = ( -2*( float(ll))) + (2*5)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


        
    ###########################################################################
    #ISLAND Discrete admixture models
    ###########################################################################

    
    elif model_name == "vic_no_mig_admix_early":
        #####################################
        #Vicariance with early discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with early discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_early

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
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

    
    elif model_name == "vic_no_mig_admix_late":
        #####################################
        #Divergence with late discrete admixture
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with late discrete admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_no_mig_admix_late

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
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


    elif model_name == "vic_two_epoch_admix":
        #####################################
        #Vicariance with discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Vicariance with discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.vic_two_epoch_admix

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
        aic = ( -2*( float(ll))) + (2*7)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


    elif model_name == "founder_nomig_admix_early":
        #####################################
        #Founder event and discrete early admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete early admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_early

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
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model
    
    elif model_name == "founder_nomig_admix_late":
        #####################################
        #Founder event and discrete late admixture
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete late admixture",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_late

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
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        return sim_model

    elif model_name == "founder_nomig_admix_two_epoch":
        #####################################
        #Founder event and discrete admixture, two epoch
        #####################################
        print "---------------------------------------------------"
        print "Founder event and discrete admixture, two epoch",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.founder_nomig_admix_two_epoch

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
        aic = ( -2*( float(ll))) + (2*7)
        print "AIC = ", aic, '\n', '\n'
        return sim_model


#======================================================================================
#======================================================================================
#======================================================================================
#======================================================================================
