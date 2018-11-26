**Demographic Modeling With dadi**
---------------------------------

**Purpose:**

Perform demographic model optimizations and comparisons with this accessible and flexible dadi tool.

This tool is designed to work with the Python package [dadi](https://bitbucket.org/gutenkunstlab/dadi) 
and assumes you already have the package installed. You'll need to be familiar with how dadi works, 
and some of the basic syntax for writing dadi scripts with python. A good resource for all dadi-related 
questions is the [user group](https://groups.google.com/forum/#!forum/dadi-user). Before attempting
to use these scripts, read over the user manual for dadi and try running the program with the 
example files.

**Overview:**

In this main repository is a general use script that can be used to run dadi to fit any model on an allele frequency spectrum/joint-site frequency spectrum containing one to three populations. To use this workflow, you'll need a SNPs input text file to create an allele frequency or joint site frequency spectrum object. Alternatively, you can import a frequency spectrum of your own creation, editing the script appropriately (see dadi manual). The user will have to edit information about their allele frequency spectrum, and a #************** marks lines in the *dadi_Run_Optimizations.py* that will have to be edited. Any custom model can be used, and below are several examples of how to use various arguments to control the model optimizations. 

The *dadi_Run_Optimizations.py* script and *Optimize_Functions.py* script must be in the same working directory to run properly.

If you'd like to use the optimization routine of this script to analyze larger sets of published models through established pipelines, please look in the nested repositories ([Two_Population_Pipeline](https://github.com/dportik/dadi_pipeline/tree/master/Two_Population_Pipeline), [Three_Population_Pipeline](https://github.com/dportik/dadi_pipeline/tree/master/Three_Population_Pipeline)) to see how to import models from external model scripts for two population or three population comparisons.
There are a considerable number of 2D models that can be selected from (32), with a fewer number of 3D models (18).

If you'd like to assess the goodness of fit for your demographic model, please look in the [Goodness_of_Fit](https://github.com/dportik/dadi_pipeline/tree/master/Goodness_of_Fit) repository.

If you'd like to create a figure comparing the empirical SFS and model SFS for a demographic model (with residuals), please look in the [Plotting](https://github.com/dportik/dadi_pipeline/tree/master/Plotting) repository.


**Optimizations:**

The *dadi_Run_Optimizations.py* script will run the optimization routine, which contains a user-defined number of rounds, each with a user-defined or default number of replicates. The starting parameters are initially random, but after each round is complete the parameters of the best scoring replicate from that round are used to generate perturbed starting parameters for the replicates of the subsequent round. The arguments controlling steps of the optimization algorithm (maxiter) and perturbation of starting parameters (fold) can be supplied by the user for more control across rounds. The user can also supply their own set of initial parameters, or set custom bounds on the parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility should allow these scripts to be generally useful for fitting any model to any data set.

**Examples of Usage:**

Let's assume you've supplied the correct information about your SNPs input file, population IDs, projection sizes, and are using the model in the script (sym_mig).

I will show several ways to use the main function for model fitting to highlight different options. 

We will use always use the following function from the Optimize_Functions.py script, which requires some explanation:

***Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")***
 
***Mandatory Arguments:***

+ **fs**:  spectrum object name
+ **pts**: grid size for extrapolation, list of three values
+ **outfile**:  prefix for output naming
+ **model_name**: a label help name the output files; ex. "no_mig"
+ **func**: access the model function from within 'dadi_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
+ **rounds**: number of optimization rounds to perform
+ **param_number**: number of parameters in the model selected (can count in params line for the model)
+ **fs_folded**: A Boolean value indicating whether the empirical fs is folded (True) or not (False)

***Optional Arguments:***

+ **reps**: a list of integers controlling the number of replicates in each of the optimization rounds
+ **maxiters**: a list of integers controlling the maxiter argument in each of the optimization rounds
+ **folds**: a list of integers controlling the fold argument when perturbing input parameter values
+ **in_params**: a list of parameter values 
+ **in_upper**: a list of upper bound values
+ **in_lower**: a list of lower bound values
+ **param_labels**: list of labels for parameters that will be written to the output file to keep track of their order


***Example 1***

Let's use the function to run an optimization routine for our data and this model.
We need to specify the first eight required arguments in this function, but there are other options
we can also use if we wanted more control over the optimization scheme. We'll start with
the basic version here. The argument explanations are above. This would perform three
rounds of optimizations, using a default number of replicates for each round. At the end
of each round, the parameters of the best-scoring replicate are used to generate the starting
parameters for the replicates of the next round, and so on. This will help focus the parameter
search space as the rounds continue.

    #create a prefix to label the output files
    prefix = "V1"
    #make sure to define your extrapolation grid size
    pts = [50,60,70]
    
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True)

***Example 2***

It is a good idea to include the labels of the parameters so they can get written to the
output file, otherwise you'll have to go back to the model each time you wanted to see their
order.

    prefix = "V2"
    pts = [50,60,70]
    
    #here are the labels, given as a string
    p_labels = "nu1, nu2, m, T"
    
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, param_labels = p_labels)


***Example 3***

Here is the same example but also including your own custom parameter bounds. Notice
the optional arguments can be placed in any order following the mandatory arguments.

    prefix = "V3"
    pts = [50,60,70]
    p_labels = "nu1, nu2, m, T"

    #Here are the custom bounds
    upper = [20,20,10,15]
    lower = [0.01,0.01,0.01,0.1]

    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, param_labels = p_labels, in_upper = upper, in_lower = lower)

***Example 4***

You can also be very explicit about the optimization routine, controlling what happens
across each round. Let's keep the three rounds, but change the number of replicates,
the maxiter argument, and fold argument each time. We'll need to create a list of values
for each of these, that has three values within (to match three rounds).

    prefix = "V4"
    pts = [50,60,70]
    p_labels = "nu1, nu2, m, T"
    upper = [20,20,10,15]
    lower = [0.01,0.01,0.01,0.1]
    
    #Here are the optional arguments
    reps = [10,20,50]
    maxiters = [5,10,20]
    folds = [3,2,1]

    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

Using these arguments will cause round one to have 10 replicates, use 3-fold perturbed
starting parameters, and a maxiter of 5 for the optimization algorithm steps. Round two
will have 20 replicates, use 2-fold perturbed starting parameters, and a maxiter of 10
for the optimization algorithm steps, and etc. for round three. 


***Example 5***

It's also good run the optimization routine multiple times. Let's write a short
loop to do the above optimization routine five times. We will name the prefix based
on which point we are at, and include it within the looping. Note that when you use
the range argument in python it will go up to, but not include, the final number.
That's why I have written a range of 1-6 to perform this 5 times.

    pts = [50,60,70]
    p_labels = "nu1, nu2, m, T"
    upper = [20,20,10,15]
    lower = [0.01,0.01,0.01,0.1]
    reps = [10,20,50]
    maxiters = [5,10,20]
    folds = [3,2,1]

    for i in range(1,6):
        prefix = "V5_Number_{}".format(i)
        Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)


**Test Data Set:**

In the folder labeled *Example_Data* you will find a SNPs input file that will run with the *dadi_Run_Optimizations.py* script.
You will only need to edit the path to the file in the script, and then you will be able to run all five examples above. The 
outputs for these examples are also contained within the *Example_Data* folder, in a separate folder labeled *Example_Outputs*.
Please test the script using these data to ensure everything is working properly before examining your own empirical data. 


**Outputs:**

 For each model run, there will be a log file showing the optimization steps per replicate and a summary file that has all the important information. 
 
Here is an example of the output from a summary file, which will be in tab-delimited format:

     Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T) 
     sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
     sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
     sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
     sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
     sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688

**Using Folded vs. Unfolded Spectra:**

 To change whether the frequency spectrum is folded vs. unfolded requires two changes in the script. The first is where the spectrum object is created, indicated by the *polarized* argument:
 
     #Convert this dictionary into folded AFS object
     #[polarized = False] creates folded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

The above code will create a folded spectrum. When calling the optimization function, this must also be indicated in the *fs_folded* argument:

     #this is from the first example:
     Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True)
     
To create an unfolded spectrum, the *polarized* and *fs_folded*  arguments in the above lines need to be changed accordingly:

     #[polarized = True] creates an unfolded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
     
     #and the optimization routine function must also be changed:
     Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=False)
     
It will be clear if either argument has been misspecified because the calculation of certain statistics will cause a crash with the following error:

     ValueError: Cannot operate with a folded Spectrum and an unfolded one.

If you see this, check to make sure both relevant arguments actually agree on the spectrum being folded or unfolded.

**Default Settings:**

The optimization routine arguments offer a lot of flexibility, but the default settings can also be used. If only
the number of rounds is changed, here are the defaults for the optional arguments (reps, maxiters, folds)
based on the number of rounds selected:

**Three rounds (as in Examples 1-3):**

| Argument | Round 1 | Round 2  | Round 3|
| ------ |------:| -----:| -----:|
| ***reps***    | 10 | 10 | 20 |
| ***maxiter*** | 5 |  5  | 5 |
| ***fold*** |  3 |  2   | 1 |


**Two rounds:**

| Argument | Round 1 | Round 2  |
| ------ |------:| -----:|
| ***reps***    | 10 | 20 |
| ***maxiter*** | 5  | 5 |
| ***fold*** |  2   | 1 |

**X Rounds (>3):**

| Argument | Round 1 | Round 2  | Round 3| Round *X-1* | Round *X* |
| ------ |------:| -----:| -----:| -----:| -----:|
| ***reps***    | 10 | 10 | 10 | 10 | 20 |
| ***maxiter*** | 5 |  5  | 5 | 5 | 5 |
| ***fold*** |  3 |  3  | 3 | 2 | 1 |

In general, you should probably run multiple rounds and ensure the log-likelihoods are converging.

**Caveats:**

 The likelihood and AIC returned represent the true likelihood only if the SNPs are unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this is considered true, but if SNPs are linked across loci then the likelihood is actually a composite likelihood and using something like AIC is no longer appropriate for model comparisons. See the discussion group for more information on this subject. 


**Citation Information:**

***Using the pipeline.***
The scripts involved with this pipeline were originally published as part of the following work:

+ *Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology 26: 5245-5263. https://doi.org/10.1111/mec.14266*

If you use or modify these scripts for your own purposes, please cite the publication.


**Publications using this workflow:**

+ Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., Portik, D.M., Streicher, J.W., and S.P. Loader. Vanishing refuge: testing the forest refuge hypothesis in coastal East Africa using genome-wide sequence data for five co-distributed amphibians. ***Molecular Ecology, Early Access.*** *https://doi.org/10.1111/mec.14862*

+ Schield, D.R., Adams, R.H., Card, D.C., Corbin, A.C., Jezkova, T., Hales, N.R., Meik, J.M., Perry, B.W., Spencer, C.L., Smith, L.L., Garcia, G.C., Bouzid, N.M., Strickland, J.L., Parkinson, C.L., Borja, M., Castañeda-Gaytán, G., Bryson, R.W., Flores-Villela, O.A., Mackessy, S.P., and T.A. Castoe. 2018. Cryptic genetic diversity, population structure, and gene flow in the Mojave rattlesnake (Crotalus scutulatus). ***Molecular Phylogenetics and Evolution*** 127: 669-681. *https://doi.org/10.1016/j.ympev.2018.06.013*

+ Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K., Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. 2018. Sky, sea, and forest islands: diversification in the African leaf-folding frog Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). ***Journal of Biogeography*** 45: 1781-1794. *https://doi.org/10.1111/jbi.13365*

+ Schield, D.R., Adams, R.H., Card, D.C., Perry, B.W., Pasquesi, G.M., Jezkova, T., Portik, D.M., Andrew, A.L., Spencer, C.L., Sanchez, E.E., Fujita, M.K., Mackessy, S.P., and T.A. Castoe. 2017. Genomic patterns of divergence and admixture in a widely-distributed rattlesnake provide insight into speciation with gene flow. ***Ecology and Evolution*** 7: 3951–3966. *https://doi.org/10.1002/ece3.2996*

+ Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. ***Molecular Ecology*** 26: 5245-5263. *https://doi.org/10.1111/mec.14266*


**License:**

GNU Lesser General Public License v3.0

**Contact:**

Daniel Portik, PhD

Postdoctoral Researcher

University of Arizona

daniel.portik@gmail.com

