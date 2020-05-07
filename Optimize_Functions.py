'''
-------------------------
Written for Python 2.7 and 3.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated September 2019
'''
import sys
import os
import numpy
import dadi
from datetime import datetime

def parse_params(param_number, in_params=None, in_upper=None, in_lower=None):
    """    
    Function to correctly deal with parameters and bounds, and if none were provided, 
    to generate them automatically.
    
    Arguments
    param_number: number of parameters in the model selected (can count in params line for the model)
    in_params: a list of parameter values 
    in_upper: a list of upper bound values
    in_lower: a list of lower bound values
    """
    param_number = int(param_number)
    
    #param set
    if in_params is None:
        params = [1] * param_number
    elif len(in_params) != param_number:
        raise ValueError("Set of input parameters does not contain the correct number of values: {}".format(param_number))
    else:
        params = in_params
        
    #upper bound    
    if in_upper is None:
        upper_bound = [30] * param_number
    elif len(in_upper) != param_number:
        raise ValueError("Upper bound set for parameters does not contain the correct number of values: {}".format(param_number))
    else:
        upper_bound = in_upper
        
    #lower bounds
    if in_lower is None:
        lower_bound = [0.01] * param_number
    elif len(in_lower) != param_number:
        raise ValueError("Lower bound set for parameters does not contain the correct number of values: {}".format(param_number))
    else:
        lower_bound = in_lower
        
    return params, upper_bound, lower_bound

def parse_opt_settings(rounds, reps=None, maxiters=None, folds=None):
    """    
    Function to correctly deal with replicate numbers, maxiter and fold args.
    
    Arguments
    rounds: number of optimization rounds to perform
    reps: a list of integers controlling the number of replicates in each of three optimization rounds
    maxiters: a list of integers controlling the maxiter argument in each of three optimization rounds
    folds: a list of integers controlling the fold argument when perturbing input parameter values
    """    
    rounds = int(rounds)
    
    #rep set
    #create scheme where final replicates will be 20, and all previous 10
    if reps is None:
        if rounds >= 2:
            reps_list = [10] * (rounds-1)
            reps_list.insert(len(reps_list),20)
        else:
            reps_list = [10] * rounds
    elif len(reps) != rounds:
        raise ValueError("List length of replicate values does match the number of rounds: {}".format(rounds))
    else:
        reps_list = reps
        
    #maxiters   
    if maxiters is None:
        maxiters_list = [5] * rounds
    elif len(maxiters) != rounds:
        raise ValueError("List length of maxiter values does match the number of rounds: {}".format(rounds))
    else:
        maxiters_list = maxiters
        
    #folds
    #create scheme so if rounds is greater than three, will always end with two fold and then one fold
    if folds is None:
        if rounds >= 3:
            folds_list = [3] * (rounds-2)
            folds_list.insert(len(folds_list),2)
            folds_list.insert(len(folds_list),1)
        elif rounds == 2:
            folds_list = [2] * (rounds-1)
            folds_list.insert(len(folds_list),1)
        else:
            folds_list = [2] * rounds
    elif len(folds) != rounds:
        raise ValueError("List length of fold values does match the number of rounds: {}".format(rounds))
    else:
        folds_list = folds
        
    return reps_list, maxiters_list, folds_list

def collect_results(fs, sim_model, params_opt, roundrep, fs_folded):
    """    
    Gather up a bunch of results, return a list with following elements: 
    [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values] 
    
    Arguments
    fs: spectrum object name
    sim_model: model fit with optimized parameters
    params_opt: list of the optimized parameters
    fs_folded: a Boolean (True, False) for whether empirical spectrum is folded or not
    """    

    #calculate theta
    theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
    theta = numpy.around(theta, 2)
    print("\t\t\tTheta = {:,}".format(theta))
    
    #calculate likelihood
    ll = dadi.Inference.ll_multinom(sim_model, fs)
    ll = numpy.around(ll, 2)
    print("\t\t\tLikelihood = {:,}".format(ll))

    #calculate AIC 
    aic = ( -2*( float(ll))) + (2*len(params_opt))
    print("\t\t\tAIC = {:,}".format(aic))

    #get Chi^2
    scaled_sim_model = sim_model*theta
    if fs_folded is True:
        #calculate Chi^2 statistic for folded
        folded_sim_model = scaled_sim_model.fold()
        chi2 = numpy.sum((folded_sim_model - fs)**2/folded_sim_model)
        chi2 = numpy.around(chi2, 2)
    elif fs_folded is False:
        #calculate Chi^2 statistic for unfolded
        chi2 = numpy.sum((scaled_sim_model - fs)**2/scaled_sim_model)
        chi2 = numpy.around(chi2, 2)
    print("\t\t\tChi-Squared = {:,}".format(chi2))        

    #store key results in temporary sublist, append to larger results list
    temp_results = [roundrep, ll, aic, chi2, theta, params_opt]

    return temp_results

def write_log(outfile, model_name, rep_results, roundrep):
    """    
    Reproduce replicate log to bigger log file, because constantly re-written.
    
    Arguments
    outfile: prefix for output naming
    model_name: a label to slap on the output files; ex. "no_mig"
    rep_results: the list returned by collect_results function: 
                 [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values]
    roundrep: name of replicate (ex, "Round_1_Replicate_10")
    """    
    fh_log = open("{0}.{1}.log.txt".format(outfile, model_name), 'a')
    fh_log.write("\n{}\n".format(roundrep))
    templogname = "{}.log.txt".format(model_name)
    try:
        fh_templog = open(templogname, 'r')
        for line in fh_templog:
            fh_log.write(line)
        fh_templog.close()
    except IOError:
        print("Nothing written to log file this replicate...")
    fh_log.write("likelihood = {}\n".format(rep_results[1]))
    fh_log.write("theta = {}\n".format(rep_results[4]))
    fh_log.write("Optimized parameters = {}\n".format(rep_results[5]))
    fh_log.close()

def Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True,
                         reps=None, maxiters=None, folds=None, in_params=None,
                         in_upper=None, in_lower=None, param_labels=None, optimizer="log_fmin"):
    """
    Main function for running dadi routine.

    Mandatory/Positional Arguments
    (1) fs:  spectrum object name
    (2) pts: grid size for extrapolation, list of three values
    (3) outfile:  prefix for output naming
    (4) model_name: a label to slap on the output files; ex. "no_mig"
    (5) func: access the model function from within 'moments_optimize.py' or from a separate python model script, ex. Models_2D.no_mig
    (6) rounds: number of optimization rounds to perform
    (7) param_number: number of parameters in the model selected (can count in params line for the model)
    (8) fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False). Default is True.

    Optional Arguments
    (9) reps: a list of integers controlling the number of replicates in each of three optimization rounds
    (10) maxiters: a list of integers controlling the maxiter argument in each of three optimization rounds
    (11) folds: a list of integers controlling the fold argument when perturbing input parameter values
    (12) in_params: a list of parameter values 
    (13) in_upper: a list of upper bound values
    (14) in_lower: a list of lower bound values
    (15) param_labels: a string, labels for parameters that will be written to the output file to keep track of their order
    (16) optimizer: a string, to select the optimizer. Choices include: log (BFGS method), 
                    log_lbfgsb (L-BFGS-B method), log_fmin (Nelder-Mead method, DEFAULT), and log_powell (Powell's method).
    """    

    #call function that determines if our params and bounds have been set or need to be generated for us
    params, upper_bound, lower_bound = parse_params(param_number, in_params, in_upper, in_lower)

    #call function that determines if our replicates, maxiter, and fold have been set or need to be generated for us
    reps_list, maxiters_list, folds_list = parse_opt_settings(rounds, reps, maxiters, folds)
    
    print("\n\n============================================================================"
              "\nModel {}\n============================================================================\n\n".format(model_name))

    #start keeping track of time it takes to complete optimizations for this model
    tbr = datetime.now()

    #optimizer dict
    optdict = {"log":"BFGS method", "log_lbfgsb":"L-BFGS-B method", "log_fmin":"Nelder-Mead method", "log_powell":"Powell's method"}
    
    # We need an output file that will store all summary info for each replicate, across rounds
    outname = "{0}.{1}.optimized.txt".format(outfile, model_name)
    with open(outname, 'a') as fh_out:
        if param_labels:
            fh_out.write("Model\tReplicate\tlog-likelihood\tAIC\tchi-squared\ttheta\toptimized_params({})\n".format(param_labels))
        else:
            fh_out.write("Model\tReplicate\tlog-likelihood\tAIC\tchi-squared\ttheta\toptimized_params\n")
        
    #Create list to store sublists of [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values] for every replicate
    results_list = []
    
    #for every round, execute the assigned number of replicates with other round-defined args (maxiter, fold, best_params)
    rounds = int(rounds)
    for r in range(rounds):
        print("\tBeginning Optimizations for Round {}:".format(r+1))
       
        #make sure first round params are assigned (either user input or auto generated)
        if r == int(0):
            best_params = params
        #and that all subsequent rounds use the params from a previous best scoring replicate
        else:
            best_params = results_list[0][5]

        #perform an optimization routine for each rep number in this round number
        for rep in range(1, (reps_list[r]+1) ):
            print("\n\t\tRound {0} Replicate {1} of {2}:".format(r+1, rep, (reps_list[r])))
                
            #keep track of start time for rep
            tb_rep = datetime.now()
            
            #create an extrapolating function 
            func_exec = dadi.Numerics.make_extrap_log_func(func)
            
            #perturb starting parameters
            params_perturbed = dadi.Misc.perturb_params(best_params, fold=folds_list[r],
                                                            upper_bound=upper_bound, lower_bound=lower_bound)
            
            if param_labels:
                print("\n\t\t\tModel parameters = {}".format(param_labels))
                print("\t\t\tStarting parameters = [{}]".format(", ".join([str(numpy.around(x, 6)) for x in params_perturbed])))
            else:
                print("\n\t\t\tStarting parameters = [{}]".format(", ".join([str(numpy.around(x, 6)) for x in params_perturbed])))
                
            #optimize from perturbed parameters
            if optimizer == "log_fmin":
                params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,
                                                                  lower_bound=lower_bound, upper_bound=upper_bound,
                                                                  verbose=1, maxiter=maxiters_list[r],
                                                                  output_file = "{}.log.txt".format(model_name))
            elif optimizer == "log":
                params_opt = dadi.Inference.optimize_log(params_perturbed, fs, func_exec, pts,
                                                                  lower_bound=lower_bound, upper_bound=upper_bound,
                                                                  verbose=1, maxiter=maxiters_list[r],
                                                                  output_file = "{}.log.txt".format(model_name))
            elif optimizer == "log_lbfgsb":
                params_opt = dadi.Inference.optimize_log_lbfgsb(params_perturbed, fs, func_exec, pts,
                                                                  lower_bound=lower_bound, upper_bound=upper_bound,
                                                                  verbose=1, maxiter=maxiters_list[r],
                                                                  output_file = "{}.log.txt".format(model_name))
            elif optimizer == "log_powell":
                params_opt = dadi.Inference.optimize_log_powell(params_perturbed, fs, func_exec, pts,
                                                                  lower_bound=lower_bound, upper_bound=upper_bound,
                                                                  verbose=1, maxiter=maxiters_list[r],
                                                                  output_file = "{}.log.txt".format(model_name))
            else:
                 raise ValueError("\n\nERROR: Unrecognized optimizer option: {}\nPlease select from: log, log_lbfgsb, log_fmin, or log_powell.\n\n".format(optimizer))
                 
            print("\t\t\tOptimized parameters =[{}]".format(", ".join([str(numpy.around(x, 6)) for x in params_opt])))
            print("\t\t\tOptimized using: {0} ({1})\n".format(optimizer, optdict[optimizer]))
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #collect results into a list using function above - [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values]
            roundrep = "Round_{0}_Replicate_{1}".format(r+1, rep)
            rep_results = collect_results(fs, sim_model, params_opt, roundrep, fs_folded)
            
            #reproduce replicate log to bigger log file, because constantly re-written
            write_log(outfile, model_name, rep_results, roundrep)
            
            #append results from this sim to larger list
            results_list.append(rep_results)
            
            #write all this info to our main results file
            with open(outname, 'a') as fh_out:
                #join the param values together with commas
                easy_p = ",".join([str(numpy.around(x, 4)) for x in rep_results[5]])
                fh_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(model_name, rep_results[0],
                                                                              rep_results[1], rep_results[2],
                                                                              rep_results[3], rep_results[4],
                                                                              easy_p))

            #calculate elapsed time for replicate
            tf_rep = datetime.now()
            print("\n\t\t\tReplicate time: {0} (H:M:S)\n".format(tf_rep - tb_rep))

        #Now that this round is over, sort results in order of likelihood score
        #we'll use the parameters from the best rep to start the next round as the loop continues
        results_list.sort(key=lambda x: float(x[1]), reverse=True)
        print("\n\t----------------------------------------------\n"
                  "\tBest replicate: {0}\n"
                  "\t\tLikelihood = {1:,}\n\t\tAIC = {2:,}\n"
                  "\t\tChi-Squared = {3:,}\n\t\tParams = [{4}]\n"
                  "\t----------------------------------------------\n\n".format(results_list[0][0],
                                                                              results_list[0][1],
                                                                              results_list[0][2],
                                                                              results_list[0][3],
                                                                              ", ".join([str(numpy.around(x, 4)) for x in results_list[0][5]])))

    #Now that all rounds are over, calculate elapsed time for the whole model
    tfr = datetime.now()
    print("\nAnalysis Time for Model '{0}': {1} (H:M:S)\n\n"
              "============================================================================".format(model_name, tfr - tbr))

    #cleanup file
    os.remove("{}.log.txt".format(model_name))
