import sys
import os
import numpy
import dadi
from datetime import datetime

def parse_params(param_number, in_params=None, in_upper=None, in_lower=None):
    #--------------------------------------------------------------------------------------
    # function to correctly deal with parameters and bounds, if none were provided, generate them automatically
    
    # Arguments
    # param_number: number of parameters in the model selected (can count in params line for the model)
    # in_params: a list of parameter values 
    # in_upper: a list of upper bound values
    # in_lower: a list of lower bound values
    #--------------------------------------------------------------------------------------
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
    #send back values
    return params, upper_bound, lower_bound

def parse_opt_settings(rounds, reps=None, maxiters=None, folds=None):
    #--------------------------------------------------------------------------------------
    # function to correctly deal with replicate numbers, maxiter and fold args
    
    # Arguments
    # rounds: number of optimization rounds to perform
    # reps: a list of integers controlling the number of replicates in each of three optimization rounds
    # maxiters: a list of integers controlling the maxiter argument in each of three optimization rounds
    # folds: a list of integers controlling the fold argument when perturbing input parameter values
    #--------------------------------------------------------------------------------------
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
    #send back values
    return reps_list, maxiters_list, folds_list

def collect_results(fs, sim_model, params_opt, roundrep, fs_folded):
    #--------------------------------------------------------------------------------------
    # gather up a bunch of results, return a list = [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values, sfs_sum] 
    
    # Arguments
    # fs: spectrum object name
    # sim_model: model fit with optimized parameters
    # params_opt: list of the optimized parameters
    # fs_folded: a Boolean (True, False) for whether empirical spectrum is folded or not
    #--------------------------------------------------------------------------------------

    #create empty list to store results
    temp_results = []

    #calculate likelihood
    ll = dadi.Inference.ll_multinom(sim_model, fs)
    ll = numpy.around(ll, 2)
    print "\t\t\tLikelihood = ", ll

    #calculate AIC 
    aic = ( -2*( float(ll))) + (2*len(params_opt))
    print "\t\t\tAIC = ", aic

    #calculate theta
    theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
    theta = numpy.around(theta, 2)
    print "\t\t\tTheta = ", theta

    if fs_folded is True:
        #calculate Chi^2 statistic for folded
        scaled_sim_model = sim_model*theta
        folded_sim_model = scaled_sim_model.fold()
        chi2 = numpy.sum((folded_sim_model - fs)**2/folded_sim_model)
        chi2 = numpy.around(chi2, 2)
        print "\t\t\tChi-Squared = ", chi2
        
    elif fs_folded is False:
        #calculate Chi^2 statistic for unfolded
        scaled_sim_model = sim_model*theta
        chi2 = numpy.sum((scaled_sim_model - fs)**2/scaled_sim_model)
        chi2 = numpy.around(chi2, 2)
        print "\t\t\tChi-Squared = ", chi2
        
    #get sum of sfs
    sfs_sum = numpy.around(fs.S(), 2)

    #store key results in temporary sublist, append to larger results list
    temp_results.append(roundrep)
    temp_results.append(ll)
    temp_results.append(aic)
    temp_results.append(chi2)
    temp_results.append(theta)
    temp_results.append(params_opt)
    temp_results.append(sfs_sum)

    #send list of results back
    return temp_results

def write_log(outfile, model_name, rep_results, roundrep):
    #--------------------------------------------------------------------------------------
    #reproduce replicate log to bigger log file, because constantly re-written
    
    # Arguments =
    # outfile: prefix for output naming
    # model_name: a label to slap on the output files; ex. "no_mig"
    # rep_results: the list returned by collect_results function: [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values]
    # roundrep: name of replicate (ex, "Round_1_Replicate_10")
    #--------------------------------------------------------------------------------------
    fh_log = open("{0}.{1}.log.txt".format(outfile, model_name), 'a')
    fh_log.write("\n{}\n".format(roundrep))
    templogname = "{}.log.txt".format(model_name)
    try:
        fh_templog = open(templogname, 'r')
        for line in fh_templog:
            fh_log.write(line)
        fh_templog.close()
    except IOError:
        print "Nothing written to log file this replicate..."
    fh_log.write("likelihood = {}\n".format(rep_results[1]))
    fh_log.write("theta = {}\n".format(rep_results[4]))
    fh_log.write("Optimized parameters = {}\n".format(rep_results[5]))
    fh_log.close()

def Optimize_Routine_GOF(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" "):
    #--------------------------------------------------------------------------------------
    # Mandatory Arguments =
    #(1) fs:  spectrum object name
    #(2) pts: grid size for extrapolation, list of three values
    #(3) outfile:  prefix for output naming
    #(4) model_name: a label to slap on the output files; ex. "no_mig"
    #(5) func: access the model function from within 'moments_optimize.py' or from a separate python model script, ex. Models_2D.no_mig
    #(6) rounds: number of optimization rounds to perform
    #(7) param_number: number of parameters in the model selected (can count in params line for the model)
    #(8) fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).

    # Optional Arguments =
    #(9) reps: a list of integers controlling the number of replicates in each of three optimization rounds
    #(10) maxiters: a list of integers controlling the maxiter argument in each of three optimization rounds
    #(11) folds: a list of integers controlling the fold argument when perturbing input parameter values
    #(12) in_params: a list of parameter values 
    #(13) in_upper: a list of upper bound values
    #(14) in_lower: a list of lower bound values
    #(15) param_labels: list of labels for parameters that will be written to the output file to keep track of their order
    #--------------------------------------------------------------------------------------

    #call function that determines if our params and bounds have been set or need to be generated for us
    params, upper_bound, lower_bound = parse_params(param_number, in_params, in_upper, in_lower)

    #call function that determines if our replicates, maxiter, and fold have been set or need to be generated for us
    reps_list, maxiters_list, folds_list = parse_opt_settings(rounds, reps, maxiters, folds)
    
    #start keeping track of time it takes to complete optimizations for this model
    tb_round = datetime.now()
    
    # We need an output file that will store all summary info for each replicate, across rounds
    outname = "{0}.{1}.optimized.txt".format(outfile,model_name)
    fh_out = open(outname, 'a')
    fh_out.write("Model"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"AIC"+'\t'+"chi-squared"+'\t'+"theta"+'\t'+"optimized_params({})".format(param_labels)+'\n')
    fh_out.close()
    
    #Create list to store sublists of [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values] for every replicate
    results_list = []
    
    #for every round, execute the assigned number of replicates with other round-defined args (maxiter, fold, best_params)
    rounds = int(rounds)
    for r in range(rounds):
        print "\tBeginning Optimizations for Round {}:".format(r+1)
       
        #make sure first round params are assigned (either user input or auto generated)
        if r == int(0):
            best_params = params
        #and that all subsequent rounds use the params from a previous best scoring replicate
        else:
            best_params = results_list[0][5]

        #perform an optimization routine for each rep number in this round number
        for rep in range(1, (reps_list[r]+1) ):
            print "\t\tRound {0} Replicate {1} of {2}:".format(r+1, rep, (reps_list[r]))
            
            #keep track of start time for rep
            tb_rep = datetime.now()
            
            #create an extrapolating function 
            func_exec = dadi.Numerics.make_extrap_log_func(func)

            
            #perturb starting parameters
            params_perturbed = dadi.Misc.perturb_params(best_params, fold=folds_list[r], upper_bound=upper_bound, lower_bound=lower_bound)

            #optimize from perturbed parameters
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts, lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=maxiters_list[r], output_file = "{}.log.txt".format(model_name))
            print "\t\t\tOptimized parameters = ", params_opt

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
            fh_out = open(outname, 'a')
            #join the param values together with commas
            easy_p = ",".join(str(numpy.around(x, 4)) for x in rep_results[5])
            fh_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(model_name, rep_results[0], rep_results[1], rep_results[2], rep_results[3], rep_results[4], easy_p))
            fh_out.close()

            #calculate elapsed time for replicate
            tf_rep = datetime.now()
            te_rep = tf_rep - tb_rep
            print "\n\t\t\tReplicate time: {0} (H:M:S)\n".format(te_rep)

        #Now that this round is over, sort results in order of likelihood score, we'll use the parameters from the best rep to start the next round as the loop continues
        results_list.sort(key=lambda x: float(x[1]), reverse=True)
        print "\tBest so far: {0}, ll = {1}\n\n".format(results_list[0][0], results_list[0][1])

    
    #Now that all rounds are over, calculate elapsed time for the whole model
    tf_round = datetime.now()
    te_round = tf_round - tb_round
    print "\n{0} Analysis Time: {1} (H:M:S)\n\n============================================================================".format(outfile, te_round)

    #cleanup file
    os.remove("{}.log.txt".format(model_name))

    #most important - sort the results list to find the top replicate for this simulation and return it to use
    results_list.sort(key=lambda x: float(x[1]), reverse=True)
    #remember the format of this list: [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values]
    return results_list[0]


def Optimize_Empirical(fs, pts, outfile, model_name, func, in_params, fs_folded=True):
    #--------------------------------------------------------------------------------------
    # Mandatory Arguments =
    #(1) fs:  spectrum object name
    #(2) pts: grid size for extrapolation, list of three values
    #(3) outfile:  prefix for output naming
    #(4) model_name: a label to slap on the output files; ex. "no_mig"
    #(5) func: access the model function from within script or from a separate python model script, ex. Models_2D.no_mig
    #(6) in_params: the previously optimized parameters to use
    #(7) fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False). Default is True.
    #--------------------------------------------------------------------------------------
    
    print "============================================================================\nFitting model '{}' to empirical data...\n============================================================================\n".format(model_name)
    print "Input parameters = {}\n".format(in_params)

    # We need an output file that will store empirical results
    outname = "{0}.{1}.optimized.txt".format(outfile,model_name)
    fh_out = open(outname, 'a')
    fh_out.write("Model"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"sfs_sum"+'\t'+"chi-squared"+'\n')
    
    #create an extrapolating function 
    func_exec = dadi.Numerics.make_extrap_log_func(func)

    #simulate the model with the optimized parameters
    sim_model = func_exec(in_params, fs.sample_sizes, pts)

    #calculate sum of sfs
    sfs_sum = numpy.around(fs.S(), 2)
    
    rep_results = collect_results(fs, sim_model, in_params, "1", fs_folded)
    fh_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(model_name, rep_results[0], rep_results[1], rep_results[4], rep_results[6], rep_results[3]))
    fh_out.close()
    
    theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
    scaled_sim_model = sim_model*theta

    if fs_folded is True:
        return_model = scaled_sim_model.fold()
    elif fs_folded is False:
        return_model = scaled_sim_model
    return return_model
    

def Perform_Sims(sim_number, model_fs, pts, model_name, func, rounds, param_number, fs_folded=True, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" "):
    #--------------------------------------------------------------------------------------
	# Mandatory Arguments =
		#(1) sim_number: the number of simulations to perform
		#(2) model_fs:  the scaled model spectrum object name
		#(3) pts: grid size for extrapolation, list of three values
		#(4) model_name: a label to help label on the output files; ex. "no_mig"
		#(5) func: access the model function from within this script
		#(6) rounds: number of optimization rounds to perform
		#(7) param_number: number of parameters in the model selected (can count in params line for the model)
        #(8) fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False). Default is True.

	# Optional Arguments =
		 # reps: a list of integers controlling the number of replicates in each of the optimization rounds
		 # maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
		 # folds: a list of integers controlling the fold argument when perturbing input parameter values
    #--------------------------------------------------------------------------------------

    #Define number of simulations to perform
    sims = int(sim_number)

    #Output file information
    sim_out = "Simulation_Results.txt"
    fh_out = open(sim_out, 'a')
    fh_out.write("Simulation"+'\t'+"Best_Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"sfs_sum"+'\t'+"chi-squared"+'\t'+"optimized_params"+'\n')
    fh_out.close()

    #Simulate data sets and optimize each using the general optimization routine
    for i in range(1,(sims+1)):
        fh_out = open(sim_out, 'a')

        #create the simulated data
        sim_fs = model_fs.sample()

        #prefix for output file naming
        outfile = "Simulation_{}".format(i)
        print "\n\n============================================================================\n{}\n============================================================================".format(outfile)

        #optimize the simulated SFS
        best_rep = Optimize_Routine_GOF(sim_fs, pts, outfile, model_name, func, rounds, param_number, fs_folded, reps, maxiters)

        easy_params = ",".join(str(numpy.around(x, 4)) for x in best_rep[5])

        #list is [roundnum_repnum, log-likelihood, AIC, chi^2 test stat, theta, parameter values, sfs_sum]
        #add result for this simulation to output file
        fh_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(i,best_rep[0],best_rep[1],best_rep[4],best_rep[6],best_rep[3],easy_params))
        fh_out.close()
