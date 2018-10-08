**Performing Simulations and Goodness of Fit Tests With dadi**
---------------------------------

**Purpose:**

Perform goodness of fit tests for demographic models using dadi.

This tool is designed to work with the Python package [dadi](https://bitbucket.org/gutenkunstlab/dadi) 
and assumes you already have the package installed. You'll need to be familiar with how dadi works, 
and some of the basic syntax for writing dadi scripts with python. A good resource for all dadi-related 
questions is the [user group](https://groups.google.com/forum/#!forum/dadi-user). Before attempting
to use these scripts, read over the user manual for dadi and try running the program with the 
example files.

**Overview:**

This is meant to be a general use script to run dadi to perform simulations and goodness of fit tests for any model on an afs/jsfs with one to three populations. To use this workflow, you'll need a SNPs input text file to create an allele frequency or joint site frequency spectrum object. Alternatively, you can import a frequency spectrum of your own creation, editing the script appropriately (see dadi manual). The user will have to edit information about their allele frequency spectrum, and a #************** marks lines in the *Simulate_and_Optimize.py* that will have to be edited. 
The frequency spectrum object can be unfolded or folded, which requires minimal script changes (see Caveats section).

The user provides a model and the previously optimized parameters for their empirical 
data. The model is fit using these parameters, and the resulting model SFS is used to
generate a user-selected number of Poisson-sampled SFS (ie simulated data). If the SNPs
are unlinked, these simulations represent parametric bootstraps. For each of 
the simulated SFS, the optimization routine is performed and the best scoring replicate
is saved. The important results for such replicates include the log-likelihood and 
Pearson's chi-squared test statistic, which are used to generate a distribution of the
simulated data values to compare the empirical values to. In addition, theta, the sum of 
the SFS, and the optimized parameters are also saved.


The *Simulate_and_Optimize.py* script and *Optimize_Functions_GOF.py* script must be in the same working directory to run properly.

**Empirical Data Optimization:**

Within the *Simulate_and_Optimize.py* script, let's assume you've supplied the correct information about your SNPs input file, population IDs, projection sizes, and are using the model in the script (sym_mig).

The model will first be fit to the empirical data using the following function:

***Optimize_Empirical(fs, pts, outfile, model_name, func, in_params, fs_folded)***
 
***Mandatory Arguments:***

+ **fs**:  spectrum object name
+ **pts**: grid size for extrapolation, list of three values
+ **outfile**:  prefix for output naming
+ **model_name**: a label help name the output files; ex. "sym_mig"
+ **func**: access the model function from within 'Simulate_and_Optimize.py' or from a separate model script
+ **in_params**: the previously optimized parameter values to use
+ **fs_folded**: A Boolean value indicating whether the empirical fs is folded (True) or not (False)

***Example:***

In the script you will need to define the extrapolation grid size and the parameter values. The 
number of parameter values must match the number in the model. 

    #Make sure to define your extrapolation grid size.
    pts = [50,60,70]
    
    #Provide best optimized parameter set for empirical data.
    #These will come from previous analyses you have already completed
	emp_params = [0.1487,0.1352,0.2477,0.1877]
    
    #Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
    fs_folded = True

    #Fit the model using these parameters and return the folded model SFS (scaled by theta).
	#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model function,
	#but everything else can stay as it is. See above for argument explanations.
	scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, pts, "Empirical", "sym_mig", sym_mig, emp_params, fs_folded=fs_folded)

	
**Performing Simulations and Optimizations:**

After the model is fit to the empirical data, the model SFS can be used to generate a user-selected number of Poisson-sampled SFS,
in other words, the simulated data. For each simulation, an optimization routine is performed that is similar in structure
to that described in the original model-fitting script [here](https://github.com/dportik/dadi_pipeline). The routine contains a 
user-defined number of rounds, each with a user-defined or default number of replicates. The starting parameters are initially random, 
but after each round is complete the parameters of the best scoring replicate from that round are used to generate perturbed starting 
parameters for the replicates of the subsequent round. The arguments controlling steps of the optimization algorithm (maxiter) and 
perturbation of starting parameters (fold) can be supplied by the user for more control across rounds. 

The simulations and optimizations are performed with the following function:

***Perform_Sims(sim_number, model_fs, pts, model_name, func, rounds, param_number, fs_folded, reps=None, maxiters=None, folds=None)***
 
***Mandatory Arguments:***

+ **sim_number**: the number of simulations to perform
+ **model_fs**:  the scaled model spectrum object name
+ **pts**: grid size for extrapolation, list of three values
+ **model_name**: a label to help label on the output files; ex. "sym_mig"
+ **func**: access the model function from within this script
+ **rounds**: number of optimization rounds to perform
+ **param_number**: number of parameters in the model to fit
+ **fs_folded**: A Boolean value indicating whether the empirical fs is folded (True) or not (False)

***Optional Arguments:***

+ **reps**: a list of integers controlling the number of replicates in each of the optimization rounds
+ **maxiters**: a list of integers controlling the maxiter argument in each of the optimization rounds
+ **folds**: a list of integers controlling the fold argument when perturbing input parameter values


***Example:***

The important arguments will need to be defined in the script. Below shows how to perform
100 simulations and define an optimization routine. 

    #Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
    fs_folded = True

    #Set the number of simulations to perform here. This should be ~100 or more.
    sims = 100
    
    #Enter the number of parameters found in the model to test.
    p_num = 4
    
    #Set the number of rounds here.
    rounds = 3
    
    #I strongly recommend defining the lists for optional arguments to control the settings 
    #of the optimization routine for all the simulated data.
    reps = [20,30,50]
    maxiters = [5,10,20]
    folds = [3,2,1]
    
    #Execute the optimization routine for each of the simulated SFS.
    #Here, you will want to change the "sym_mig" and sym_mig arguments to match your model 
    #function name, but everything else can stay as it is.
    Optimize_Functions_GOF.Perform_Sims(sims, scaled_fs, pts, "sym_mig", sym_mig, rounds, p_num, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds)

The optimization routine set here will have the following settings:

| Argument | Round 1 | Round 2  | Round 3|
| ------ |------:| -----:| -----:|
| ***reps***    | 20 | 30 | 50 |
| ***maxiter*** | 5 |  10  | 20 |
| ***fold*** |  3 |  2   | 1 |

If only the number of rounds is provided, but no additional optional arguments, the optimization
routine will use the default values for each round described [here](https://github.com/dportik/dadi_pipeline).

Because it may take some time to optimize each simulated SFS, the elapsed time is provided along
the way which can help provide an estimate of the total time necessary. You may choose to adjust
the optimization routine accordingly, or change the number of simulations.

**Outputs:**

The ***Optimize_Empirical*** function will produce an output file for the empirical fit, which will be in tab-delimited format:

     Model	Replicate	log-likelihood	theta	sfs_sum	chi-squared
     sym_mig	1	-591.21	619.83	1552.44	758.21

This is based on the parameter values supplied, as no optimization routine is performed. 

The ***Perform_Sims*** function will produce many output files.
For each simulation performed, a log file and optimization summary output file will be produced 
with a prefix matching the simulation number. The optimization summary output file will be in tab-delimited format:

     Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params( )
     sym_mig	Round_1_Replicate_1	-476.15	960.3	571.75	429.99	0.33,0.2753,0.2759,0.8528
     sym_mig	Round_1_Replicate_2	-7182.39	14372.78	40146874463.2	52.58	2.72,7.1758,4.8034,4.7965
     sym_mig	Round_1_Replicate_3	-2221.59	4451.18	58293.74	109.05	4.5185,0.6713,0.4644,3.7992
     sym_mig	Round_2_Replicate_1	-570.25	1148.5	893.16	401.4	0.2539,0.4217,0.0744,0.4792
     sym_mig	Round_2_Replicate_2	-1201.35	2410.7	3134.15	340.21	1.2311,0.1402,0.4655,0.3826
     sym_mig	Round_2_Replicate_3	-621.59	1251.18	1813.97	332.51	0.5275,0.3393,0.0769,1.1639
     sym_mig	Round_3_Replicate_1	-490.62	989.24	680.65	382.17	0.4025,0.3321,0.1819,0.8315
     sym_mig	Round_3_Replicate_2	-509.41	1026.82	705.95	339.42	0.5187,0.3782,0.1482,0.8757
     sym_mig	Round_3_Replicate_3	-467.93	943.86	474.19	513.94	0.2516,0.2106,0.4328,0.8059

After all simulations are complete, the main output file *Simulation_Results.txt* will be created.
This file contains the best scoring replicate for each simulation, and contains the 
log-likelihood, theta, sum of sfs, Pearson's chi-squared test statistic, and optimized parameter
values. It will also be in tab-delimited format:

     Simulation	Best_Replicate	log-likelihood	theta	sfs_sum	chi-squared	optimized_params
     1	Round_3_Replicate_3	-467.93	513.94	1497.0	474.19	0.2516,0.2106,0.4328,0.8059
     2	Round_3_Replicate_1	-907.83	250.22	1494.0	1757.27	1.6895,0.3219,0.0868,0.7076
     3	Round_3_Replicate_2	-458.62	315.17	1508.0	455.3	0.4111,0.3406,0.2915,3.5104
     4	Round_3_Replicate_3	-488.11	133.42	1568.0	688.36	1.175,1.0293,0.0753,5.6391
     5	Round_3_Replicate_3	-456.46	621.48	1522.0	397.07	0.1981,0.164,0.631,1.8222

**Plotting Results:**

The R-script *Plot_GOF.R* can be used to visualize the results. Essentially, the simulations will
be used to create a distribution of values to which the empirical value will be
compared. The script will create a histogram of the simulated data values, and a blue line will
be plotted showing the empirical value. This will be done for the log-likelihood scores and for 
the log-transformed Pearson's chi-squared test statistic. 

The script will need to be edited to indicate the path to the *Simulation_Results.txt* file and the 
output file for the empirical data. The script will also need to be tailored to create histograms that match the empirical distributions 
for your data set. This can be done by editing the 'seq' function to create an appropriate range of 
bin sizes and overall number of bins for the histograms. The 'seq' function takes three arguments: a 
minimum range value, maximum range value, and increment value. So, seq(0,20,1) will created bins ranging 
from zero to twenty by increments of one. You'll need to adjust the range to include the values
from the simulations (and the empirical data, if it falls outside this range!). In general,
if the empirical value falls outside the simulated value distribution in the direction of worse values, the goodness of fit
test is not considered passed. 

**Using Folded vs. Unfolded Spectra:**

 To change whether the frequency spectrum is folded vs. unfolded requires two changes in the script. The first is where the spectrum object is created, indicated by the *polarized* argument:
 
     #Convert this dictionary into folded AFS object
     #[polarized = False] creates folded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

The above code will create a folded spectrum. When calling the empirical optimization function for the model, this is indicated by the *fs_folded* argument:

     #Note that in the script fs_folded is assigned a variable and referred to in the empirical and simulation optimization functions:
     
     #**************
     #Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
     fs_folded = True
     
     scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, pts, "Empirical", "sym_mig", sym_mig, emp_params, fs_folded=fs_folded)
     Optimize_Functions_GOF.Perform_Sims(sims, scaled_fs, pts, "sym_mig", sym_mig, rounds, p_num, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds)
     
To create an unfolded spectrum, the *polarized* and *fs_folded*  arguments in the above lines need to be changed accordingly:

     #[polarized = True] creates an unfolded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
     

     #Change this variable to False to set the argument fs_folded in all the empirical and simulation optimizations
     fs_folded = False
     
     scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, pts, "Empirical", "sym_mig", sym_mig, emp_params, fs_folded=fs_folded)
     Optimize_Functions_GOF.Perform_Sims(sims, scaled_fs, pts, "sym_mig", sym_mig, rounds, p_num, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds)
     
It will be clear if either argument has been misspecified because the calculation of certain statistics will cause a crash with the following error:

     ValueError: Cannot operate with a folded Spectrum and an unfolded one.


**Caveats:**

 The data are simulated using the fs.sample() method in dadi, which is equivalent to a 
 parametric boostrap ONLY if SNPs are unlinked across loci. For ddRADseq data where a 
 single SNP is selected per locus, this is generally true, and this workflow is valid.
 

**Test Data Set:**

In the folder labeled *Example_Data* you will find a SNPs input file that will run with the *Simulate_and_Optimize.py* script.
You will only need to edit the path to the file in the script, and then the script should run normally. The 
outputs for five simulations (from truncated optimizations) are contained within the *Example_Data* folder, in a separate folder labeled *Example_Outputs*.
Running the *Simulate_and_Optimize.py* script as is will actually produce 100 simulations, rather than five. 
You may choose to test the script using these data to ensure everything is working properly before examining your own empirical data. 


**Citation Information:**

 If you use or modify these scripts for your work, please cite the following publications.
 
 For the general optimization routine:

+ Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. ***Molecular Ecology*** 26: 5245-5263. *https://doi.org/10.1111/mec.14266*

For the goodness of fit testing scripts:

+ Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., Portik, D.M., Streicher, J.W., and S.P. Loader. Vanishing refuge: testing the forest refuge hypothesis in coastal East Africa using genome-wide sequence data for five co-distributed amphibians. ***Molecular Ecology, Early Access.*** *https://doi.org/10.1111/mec.14862*


**Contact:**

Daniel Portik, PhD

Postdoctoral Researcher

University of Arizona

daniel.portik@gmail.com

