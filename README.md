# Perform Demographic Modeling with dadi using dadi_pipeline

---------------------------------

## Overview of Page Contents

+ [Purpose](#P)
+ [Version](#V)
+ [Dadi Pipeline Overview](#PO)
+ [Optimization Routine Overview](#OR)
+ [Examples of Usage](#EU)
+ [Test Data](#TD)
+ [Outputs](#O)
+ [Designating Folded vs. Unfolded Spectra](#FU)
+ [Default Optimization Routine Settings](#DOR)
+ [Why Perform Multiple Rounds of Optimizations?](#WMR)
+ [My Analysis Crashed! What Now?](#AC)
+ [Reporting Bugs/Errors](#RBE)
+ [Caveats](#C)
+ [Citation Information](#CI)
+ [Publications](#PUBS)
+ [License](#LIC)
+ [Contact](#CON)

## **Purpose** <a name="P"></a>

Perform demographic model optimizations and comparisons with this accessible and flexible  tool called `dadi_pipeline`.

This tool is designed to work with the Python package [dadi](https://bitbucket.org/gutenkunstlab/dadi) and assumes you already have the package installed. You'll need to be familiar with how dadi works, and some of the basic syntax for writing dadi scripts with python. A good resource for all dadi-related questions is the [dadi user group](https://groups.google.com/forum/#!forum/dadi-user). Before attempting to use these scripts, read over the user manual for dadi and try running the program with the example files.

The most current version of this pipeline is designed to run in Python 3 and requires the Python 3 version of `dadi` (tested for v2.0.5). The pipeline has been tested with Python 3.7. If you need to run the pipeline with Python 2, you could download and use the [v3.0 release here](https://github.com/dportik/dadi_pipeline/releases/tag/v3.0), but this also requires obtaining an older version of `dadi` (v1.x). For many reasons, this is not recommended. You should use the most current versions of `dadi` and `dadi_pipeline`, and run them using Python 3.

## **Version** <a name="V"></a>

**The current version of `dadi_pipeline` is v3.1.5, and it has been tested with `dadi` v2.0.5. Please ensure you have the latest versions of both `dadi_pipeline` and `dadi` installed.**

### What's new in version 3.1+?

Version 3.1.5:
+ Added 15 new 3D models that were developed for [Firneno et al. 2020]( https://doi.org/10.1111/mec.15496). Included is a new class of model involving an admixed (or "hybrid") origin for the third population. All new models are fully integrated in the [Three_Population_Pipeline](https://github.com/dportik/dadi_pipeline/tree/master/Three_Population_Pipeline).

Version 3.1.4:
+ Added option to select which optimizer to use. Includes the following functions from the `Inference.py` module in `dadi`: optimize_log (BFGS method), optimize_log_lbfgsb (L-BFGS-B method), optimize_log_fmin (Nelder-Mead method), and optimize_log_powell (Powell's method). For usage details please see section: [Default Optimization Routine Settings](#DOR).

Version 3.1.2:
+ Issue causing crash when referencing optimized params list is fixed. Compatible with `dadi` v2.0.5.

Version 3.1.1:
+ The 2D island model set has been revised and updated, fixing issues with parameter assignment. 

Version 3.1.0:
+ **All scripts have been upgraded to Python 3 (tested with 3.7), allowing compatibility with `dadi` v2.**
+ Prints parameter labels and perturbed starting parameters to screen for each replicate. This allows quick visual comparisons of the initial and final optimized parameters among replicates.


## **Dadi Pipeline Overview** <a name="PO"></a>

In this main repository of `dadi_pipeline` is a general use script (`dadi_Run_Optimizations.py`) that can be used to run dadi to fit any model on an allele frequency spectrum/joint-site frequency spectrum containing one to three populations. This script will perform a general optimization routine proposed by [Portik et al. (2017)](https://doi.org/10.1111/mec.14266) and will produce associated output files. To use this workflow, you'll need a SNPs input text file to create an allele frequency or joint site frequency spectrum object. Alternatively, you can import a frequency spectrum of your own creation, editing the script appropriately (see dadi manual). The user will have to edit information about their allele frequency spectrum, and a #************** marks lines in the `dadi_Run_Optimizations.py` that will have to be edited. Any custom model can be used, and below are several examples of how to use various arguments to control the model optimizations. 

The `dadi_Run_Optimizations.py` script and `Optimize_Functions.py` script must be in the same working directory to run properly.

If you'd like to use the optimization routine of this script to analyze larger sets of published 2D or 3D models, please look in the nested repositories ([Two_Population_Pipeline](https://github.com/dportik/dadi_pipeline/tree/master/Two_Population_Pipeline), [Three_Population_Pipeline](https://github.com/dportik/dadi_pipeline/tree/master/Three_Population_Pipeline)). These are essentially modified versions of the `dadi_Run_Optimizations.py` script that are designed to perform the optimization routine across the available 2D or 3D models.
There are a considerable number of 2D models that can be selected from, and many 3D models too. A visual depiction of these models can be found in the [Models_2D.pdf](https://github.com/dportik/dadi_pipeline/blob/master/Two_Population_Pipeline/Models_2D.pdf) and the [Models_3D.pdf](https://github.com/dportik/dadi_pipeline/blob/master/Three_Population_Pipeline/Models_3D.pdf) files.

If you'd like to assess the goodness of fit for your demographic model, please look in the [Goodness_of_Fit](https://github.com/dportik/dadi_pipeline/tree/master/Goodness_of_Fit) repository.

If you'd like to create a figure comparing the empirical SFS and model SFS for a demographic model (with residuals), please look in the [Plotting](https://github.com/dportik/dadi_pipeline/tree/master/Plotting) repository.

For information on how to cite `dadi_pipeline`, please see the Citation section at the bottom of this page.


## **Optimization Routine Overview** <a name="OR"></a>

The `dadi_Run_Optimizations.py` and associated 2D and 3D population pipelines are components of `dadi_pipeline` that each were designed to implement the optimization routine proposed by [Portik et al. (2017)](https://doi.org/10.1111/mec.14266). This optimization routine includes fitting the model using particular settings for a given number of replicates, then using the parameters from the best scoring replicate to seed a subsequent round of model fitting using updated settings. This process occurs across multiple rounds, which improves the log-likelihood scores and generally results in convergence in the final round.

In the `dadi_Run_Optimizations.py` script, the optimization routine contains a user-defined number of rounds, each with a user-defined or default number of replicates. The starting parameters are initially random, but after each round is complete the parameters of the best scoring replicate from that round are used to generate perturbed starting parameters for the replicates of the subsequent round. The arguments controlling steps of the optimization algorithm (maxiter) and perturbation of starting parameters (fold) can be supplied by the user for more control across rounds. The user can also supply their own set of initial parameters, or set custom bounds on the parameters (upper_bound and lower_bound) to meet specific model needs. As of v3.1.3, the user can also choose which optimizer to use. This flexibility should allow these scripts to be generally useful for fitting any model to any data set. 


## **Examples of Usage** <a name="EU"></a>

Let's assume you've supplied the correct information about your SNPs input file, population IDs, projection sizes, and are using the model in the script (sym_mig).

I will show several ways to use the main function for model fitting to highlight different options. 

We will use always use the following function from the `Optimize_Functions.py` script, which requires some explanation:

`Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=None, optimizer="log_fmin")`
 
***Mandatory Arguments:***

+ **fs**:  spectrum object name
+ **pts**: grid size for extrapolation, list of three values
+ **outfile**:  prefix for output naming
+ **model_name**: a name to help label the output files - ex. "no_mig"
+ **func**: access the model function from within `dadi_Run_Optimizations.py` or from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
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
+ **optimizer**: a string, to select the optimizer. Choices include: "log" (BFGS method), "log_lbfgsb" (L-BFGS-B method), "log_fmin" (Nelder-Mead method; the default), and "log_powell" (Powell's method).

The mandatory arguments must always be included when using the `Optimize_Routine` function, and the arguments must be provided in the exact order listed above (also known as positional arguments). The optional arguments can be included in any order after the required arguments, and are referred to by their name, followed by an equal sign, followed by a value (example: `reps = 4`). The usage is explained in the following examples.

### Example 1

Let's use the function to run an optimization routine for our data and this model.
We always need to specify the eight required arguments (in order) in this function, but there are other options
we can also use if we wanted more control over the optimization scheme. We'll start with
the basic version here. The argument explanations are above. This would perform three
rounds of optimizations, using a default number of replicates for each round. For details on these defaults, please see the [Default Optimization Routine Settings](#DOR) section. At the end of each round, the parameters of the best-scoring replicate are used to generate the starting parameters for the replicates of the next round, and so on. This will help focus the parameter search space as the rounds continue. 

    #create a prefix to label the output files
    prefix = "V1"
    #make sure to define your extrapolation grid size
    pts = [50,60,70]
    
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True)

### Example 2

It is a good idea to include the labels of the parameters so they can get written to the
output file, otherwise you'll have to go back to the model each time you wanted to see their order. Notice that the optional arguments require using the `=` sign to assign a variable or value to the argument.

    prefix = "V2"
    pts = [50,60,70]
    
    #here are the labels, given as a string
    p_labels = "nu1, nu2, m, T"
    
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, param_labels = p_labels)


### Example 3

Here is the same example but also including your own custom parameter bounds. Notice
the optional arguments can be placed in any order following the mandatory arguments.

    prefix = "V3"
    pts = [50,60,70]
    p_labels = "nu1, nu2, m, T"

    #Here are the custom bounds
    upper = [20,20,10,15]
    lower = [0.01,0.01,0.01,0.1]

    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, param_labels = p_labels, in_upper = upper, in_lower = lower)

### Example 4

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


### Example 5

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


## **Test Data** <a name="TD"></a>

In the folder labeled [Example_Data](https://github.com/dportik/dadi_pipeline/tree/master/Example_Data) you will find a SNPs input file that will run with the `dadi_Run_Optimizations.py` script.
You will only need to edit the path to the file in the script, and then you will be able to run all five examples above. The 
outputs for these examples are also contained within the [Example_Data](https://github.com/dportik/dadi_pipeline/tree/master/Example_Data) folder, in a separate folder labeled [Example_Outputs](https://github.com/dportik/dadi_pipeline/tree/master/Example_Data/Example_Outputs).
Please test the script using these data to ensure everything is working properly before examining your own empirical data. 


## **Outputs** <a name="O"></a>

 For each model run, there will be a log file showing the optimization steps per replicate and a summary file that has all the important information. 
 
Here is an example of the output from a summary file, which will be in tab-delimited format:

     Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T) 
     sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
     sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
     sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
     sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
     sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688

## **Designating Folded vs. Unfolded Spectra** <a name="FU"></a>

 To change whether the frequency spectrum is folded vs. unfolded requires two changes in the script. The first is where the spectrum object is created, indicated by the `polarized` argument:
 
     #Convert this dictionary into folded AFS object
     #[polarized = False] creates folded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

The above code will create a folded spectrum. When calling the optimization function, this must also be indicated in the `fs_folded` argument:

     #this is from the first example:
     Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True)
     
To create an unfolded spectrum, the `polarized` and `fs_folded`  arguments in the above lines need to be changed accordingly:

     #[polarized = True] creates an unfolded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
     
     #and the optimization routine function must also be changed:
     Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=False)
     
It will be clear if either argument has been misspecified because the calculation of certain statistics will cause a crash with the following error:

     ValueError: Cannot operate with a folded Spectrum and an unfolded one.

If you see this, check to make sure both relevant arguments actually agree on the spectrum being folded or unfolded.

## **Default Optimization Settings** <a name="DOR"></a>

The optimization routine arguments offer a lot of flexibility, but the default settings can also be used. If only
the number of rounds is changed, here are the defaults for the optional arguments (reps, maxiters, folds)
based on the number of rounds selected:

**Three rounds (as in Examples 1-3):**

| Argument | Round 1 | Round 2  | Round 3|
| ------ |------:| -----:| -----:|
| reps    | 10 | 10 | 20 |
| maxiter | 5 |  5  | 5 |
| fold |  3 |  2   | 1 |


**Two rounds:**

| Argument | Round 1 | Round 2  |
| ------ |------:| -----:|
| reps    | 10 | 20 |
| maxiter | 5  | 5 |
| fold |  2   | 1 |

**X Rounds (>3):**

| Argument | Round 1 | Round 2  | Round 3| Round *X-1* | Round *X* |
| ------ |------:| -----:| -----:| -----:| -----:|
| reps    | 10 | 10 | 10 | 10 | 20 |
| maxiter | 5 |  5  | 5 | 5 | 5 |
| fold |  3 |  3  | 3 | 2 | 1 |


The default optimizer used is the Nelder-Mead method (`Inference.py` function `optimize_log_fmin`). This can be changed by supplying the optional optimizer argument 
with any of the following choices:

+ log - Optimize log(params) to fit model to data using the BFGS method.
+ log_lbfgsb - Optimize log(params) to fit model to data using the L-BFGS-B method.
+ log_fmin - Optimize log(params) to fit model to data using Nelder-Mead. **This is the default option.**
+ log_powell - Optimize log(params) to fit model to data using Powell's method.

In Example 1 above, we could use the L-BFGS-B method instead by using the following command:

```
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, optimizer="log_lbfgsb")
```

Or we could use Powell's method with the following command:
```
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True, optimizer="log_powell")
```


## **Why Perform Multiple Rounds of Optimizations?** <a name="WMR"></a>

When fitting demographic models, it is important to perform multiple runs and ensure that final optimizations are converging on a similar log-likelihood score. In the 2D, 3D, and custom workflows of `dadi_pipeline`, the default starting parameters used for all replicates in first round are random. After each round is completed, the parameters of the best scoring replicate from the previous round are then used to generate perturbed starting parameters for the replicates of the subsequent round. This optimization strategy of focusing the parameter search space improves the log-likelihood scores and generally results in convergence in the final round. 

Below is a summary of the log-likelihood scores obtained using the default four-round optimization settings present in the 2D pipeline. This analysis was conducted for a particular model (nomig, the simplest 2D model) using the example data provided. You can clearly see the improvement in log-likelihood scores and decrease in variation among replicates as the optimization rounds progress. 

![Rounds](https://github.com/dportik/dadi_pipeline/blob/master/Example_Data/NoMig_Zoom.png)

**Please understand that it is possible for a single execution of the pipeline to get stuck on a local optima for any given model! This is why I strongly recommend running the pipeline multiple times for a given model.** If several independent runs for this model each converge on similar log-likelihood scores in their last optimization rounds, you can be mostly confident that analyses are not getting trapped on local optima, and that the true log-likelihood has been obtained.

## **My Analysis Crashed! What Now?** <a name="AC"></a>

For various reasons, sometimes an analysis can crash. In some cases, it is not desirable to re-start a model optimization routine from scratch. You can essentially pick up where you left off through a couple of simple actions. First, you will need to find the highest scoring replicate that occurred during the round that crashed. These parameter values will be used as input parameters. Second, the number of rounds and corresponding reps, maxiters, and folds arguments will need to be adjusted to start in the later rounds.

For example, let's say the program crashed fitting a model during round 2 (of 4). You can quickly sort the output file to find the best scoring replicate:

```
Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nuA, nu2, nu3, m1, m2, T1, T2, T3)
refugia_adj_1	Round_2_Replicate_1	-2227.89	4473.78	14271.56	180.26	0.7541,2.1299,12.2678,0.8252,0.6424,0.6868,0.3189,1.3291,0.9736
refugia_adj_1	Round_2_Replicate_7	-2283.95	4585.9	3131.07	182.59	0.6738,3.0342,5.7909,0.8692,0.3357,0.346,0.2572,2.2233,1.1109
refugia_adj_1	Round_2_Replicate_9	-2297.34	4612.68	4517.14	185.11	0.798,0.9479,6.9163,2.5229,0.3895,0.235,0.2362,1.2066,0.5539
```

In this example above, we want the parameters from Round_2_Replicate1: 
`0.7541,2.1299,12.2678,0.8252,0.6424,0.6868,0.3189,1.3291,0.9736`

In the script running this model, we will need to change the following arguments:
```
#**************
#Set the number of rounds here
rounds = 4

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]
```

These need to be changed to reflect the fact that we want to 'start' in round 3 and continue to round 4. We also want to use the best parameters from round 2 to seed round 3, so we will need to add a `params` variable. The arguments can be adjusted like so:
```
#**************
#Set the number of rounds here
rounds = 2

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [30,40]
maxiters = [10,15]
folds = [2,1]
params = [0.7541,2.1299,12.2678,0.8252,0.6424,0.6868,0.3189,1.3291,0.9736]
```

Finally, in the actual call for the model we will need to add the optional flag `in_params=params` to let the routine know we are supplying the starting parameters to seed the replicates.
 
For example, add the `in_params=params` argument to this:

`Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_adj_1", Models_3D.refugia_adj_1, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, m1, m2, T1, T2, T3")`

so that it looks like this:

`Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_adj_1", Models_3D.refugia_adj_1, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=params, param_labels = "nu1, nuA, nu2, nu3, m1, m2, T1, T2, T3")`

That will allow you to more or less pick up where you left off. Please note that if running multiple models in a given script, changing the rounds, reps, maxiters, and folds arguments will affect all of them. So, it is best to isolate a single model to jump-start a crashed analysis.

## **Reporting Bugs/Errors** <a name="RBE"></a>

If you encounter any issues while using `dadi_pipeline`, it could be the result of a dadi-specific problem or an error in `dadi_pipeline`. I strongly recommend looking through the [dadi user group](https://groups.google.com/forum/#!forum/dadi-user) to see if a similar issue has been raised previously. This community resource is extremely helpful for troubleshooting. Ryan Gutenkunst has always been extremely helpful and generous with his time when it comes to resolving issues with `dadi`. However, it is important not to bother him with questions related specifically to this pipeline, as it is an independent project. If you have questions about `dadi_pipeline`, or have found an issue you believe is specific to `dadi_pipeline`, please post on the [issues page](https://github.com/dportik/dadi_pipeline/issues) here and not on the dadi user group! 

## **Caveats** <a name="C"></a>

 The likelihood and AIC returned represent the true likelihood only if the SNPs are unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this is considered true, but if SNPs are linked across loci then the likelihood is actually a composite likelihood and using something like AIC is no longer appropriate for model comparisons. See the discussion group for more information on this subject. 


## **Citation Information** <a name="CI"></a>

### How to cite dadi_pipeline:

This demographic modeling pipeline was built with a novel multi-round optimization routine, it includes many original models, and it generates custom output files. Because of these important features, `dadi_pipeline` is best considered as an additional package. It was published as part of [Portik et al. (2017)](https://doi.org/10.1111/mec.14266). If you have used `dadi_pipeline` to run your analyses, please indicate so in your publication. Here is an example of how to cite this workflow:

> To explore alternative demographic models, we used the diffusion approximation method of dadi (Gutenkunst et al. 2009) to analyze joint site frequency spectra. We fit 15 demographic models using dadi_pipeline v3.1.5 (Portik et al. 2017).

The main motivation behind the creation of this workflow was to increase transparency and reproducibility in demographic modeling. In your publication you should report the key parameters of the optimization routine. The goal is to allow other researchers to plug your data into `dadi_pipeline` and run the same analyses. For example:

> For all models, we performed consecutive rounds of optimizations. For each round, we ran multiple replicates and used parameter estimates from the best scoring replicate (highest log-likelihood) to seed searches in the following round. We used the default settings in dadi_pipeline for each round (replicates = 10, 20, 30, 40; maxiter = 3, 5, 10, 15; fold = 3, 2, 2, 1), and optimized parameters using the Nelder-Mead method (optimize_log_fmin). Across all analyses, we used the optimized parameter sets of each replicate to simulate the 3D-JSFS, and the multinomial approach was used to estimate the log-likelihood of the 3D-JSFS given the model.

The above example explains all the parameters used to run the analyses. If you change any of the default options, you should report them here in your methods section. This can include changes to the number of rounds, replicates, maxiters, folds, or other optional features (such as supplying parameter values or changing the default parameter bounds).

There are additional features of `dadi_pipeline` that were developed for other publications. For example, several 2D and 3D models were published as part of [Charles et al. (2018)](https://doi.org/10.1111/jbi.13365)[Barratt et al. (2018)](https://doi.org/10.1111/mec.14862), and [Firneno et al. 2020]( https://doi.org/10.1111/mec.15496). The goodness of fit tests were published as part of [Barratt et al. (2018)](https://doi.org/10.1111/mec.14862). Depending on what you include in your own analyses, you may also choose to cite these papers.

Here is a list of the publications mentioned above, for easy reference:

+ *Gutenkunst, R.N., Hernandez, R.D., Williamson, S.H., and C.D. Bustamante. 2009. Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data. PLoS Genetics 5: e1000695.*

+ *Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology 26: 5245-5263. https://doi.org/10.1111/mec.14266*

+ *Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K., Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. 2018. Sky, sea, and forest islands: diversification in the African leaf-folding frog Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). Journal of Biogeography 45: 1781-1794. https://doi.org/10.1111/jbi.13365*

+ *Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., Portik, D.M., Streicher, J.W., and S.P. Loader. 2018. Vanishing refuge: testing the forest refuge hypothesis in coastal East Africa using genome-wide sequence data for five co-distributed amphibians. Molecular Ecology 27: 4289-4308. https://doi.org/10.1111/mec.14862*

+ *Firneno Jr., T.J., Emery, A.H., Gerstner, B.E., Portik, D.M., Townsend, J.H., and M.K. Fujita. 2020. Mito-nuclear discordance reveals cryptic genetic diversity, introgression, and an intricate demographic history in a problematic species complex of Mesoamerican toads. Molecular Ecology, 29: 3543–3559. https://doi.org/10.1111/mec.15496*

### **Publications** <a name="PUBS"></a>

The following is a running list of publications that have used the `dadi_pipeline` demographic modeling workflow:

#### Published Articles

+ Leder, E.H., Andre, C., Le Moan, A., Topel, M., Blomberg, A., Havenhand, J.N., Lindstrom, K., Volckaert, F.A.M., Kvarnemo, C., Johannesson, K., and O. Svensson. Post-glacial establishment of locally adapted fish populations over a steep salinity gradient. **2020**. ***Journal of Evolutionary Biology***, *Early View*.  *https://doi.org/10.1111/jeb.13668*

+ Liu, C., Sarnat, E.M., Friedman, N.R., Garcia, F.H., Darwell, C., Booher, D., Kubota, Y., Mikheyev, A.S., and E.P. Economo. **2020**. Colonize, radiate, decline: Unraveling the dynamics of island community assembly with Fijian trap-jaw ants. ***Evolution*** 74: 1082-1097. *https://doi.org/10.1111/evo.13983*

+ Firneno Jr., T.J., Emery, A.H., Gerstner, B.E., Portik, D.M., Townsend, J.H., and M.K. Fujita. **2020**. Mito-nuclear discordance reveals cryptic genetic diversity, introgression, and an intricate demographic history in a problematic species complex of Mesoamerican toads. ***Molecular Ecology*** 29: 3543–3559. *https://doi.org/10.1111/mec.15496*

+ Talla, V., Pierce, A.A., Adams, K.L., de Man, T.J.B., Nallu, S., Villablanca, F.X., Kronforst, M.R., and J.C. de Roode. **2020**. Genomic evidence for gene flow between monarchs with divergent migratory phenotypes and flight performance. ***Molecular Ecology*** 29: 2567-2582. *https://doi.org/10.1111/mec.15508*

+ Ebersbach, J., Posso-Terranova, A., Bogdanowicz, S., Gomez-Diaz, M., Garcia-Gonzalez, M.X., Bolivar-Garcia, W., and J. Andres. **2020**. Complex patterns of differentiation and gene flow underly the divergence of aposematic phenotypes in *Oophaga* poison frogs. ***Molecular Ecology*** 29: 1944-1956. *https://doi.org/10.1111/mec.15360*

+ Feng, L., Ruhsam, M., Wang, Y-.H., Li, Z-.H., and X-.M. Wang. **2020**. Using demographic model selection to untangle allopatric divergence and diversification mechanisms in the *Rheum palmatum* complex in the Eastern Asiatic Region. ***Molecular Ecology*** 29: 1791-1805. 

+ Darwell, C.T., Fischer, G., Sarnat, E.M., Friedman, N.R., Liu, C., Baiao, G., Mikheyev, A.S., and E.P. Economo. **2020**. Genomic and phenomic analysis of island ant community assembly. ***Molecular Ecology*** 29: 1611-1627. *https://doi.org/10.1111/mec.15326*

+ Mattingsdal, M., Jorde, P.E., Knutsen, H., Jentoft, S., Stenseth, N.C., Sodeland, M., Robalo, J.I., Hansen, M.M., Andre, C., and E.B. Gonzalez. **2020**. Demographic history has shaped the strongly differentiated corkwing wrasse populations in Northern Europe. ***Molecular Ecology*** 29: 160-171. *https://doi.org/10.1111/mec.15310*

+ Sun, Y., and X. Wen. **2020**. Demographic inference of divergence and gene exchange between *Castanopsis fabri* and *Castanopsis lamontii*. ***Frontiers in Plant Science*** 11: 198. *https://doi.org/10.3389/fpls.2020.00198*

+ Park, B., and J.M. Burke. **2020**. Phylogeography and the evolutionary history of sunflower (*Helianthus annuus** L.): wild diversity and the dynamics of domestication. ***Genes*** 11: 266. *https://doi.org/10.3390/genes11030266**

+ Kioukis, A., Michalopoulou, V.A., Briers, L., Pirintsos, S., Studholme, D.J., Pavlidis, P., and P.F. Sarris. **2020**. Intraspecific diversification of the crop wild relative *Brassica cretica* Lam. using demographic model selection. ***BMC Genomics*** 21: 48. *https://doi.org/10.1186/s12864-019-6439-x*

+ Choi, J.Y., Purugganan, M., and E.A. Stacy. **2020**. Divergent selection and primary gene flow shape incipient speciation of a riparian tree on Hawaii Island. ***Molecular Biology and Evolution*** 37: 695-710. *https://doi.org/10.1093/molbev/msz259*

+ Beuadry, F.E.G., Barrett, S.C.H., & S.I. Wright. **2020**. Ancestral and neo-sex chromosomes contribute to population divergence in a dioecious plant. ***Evolution*** 74: 256-269. *https://doi.org/10.1111/evo.13892*

+ Salle, G., Doyle, S.R., Cabaret, J., Berriman, M., Holroyd, N., and J.A. Cotton. **2019**. The global diversity of the major parasitic nematode *Haemonchus contortus* is shaped by human intervention and climate. ***Nature Communications*** 10: 4811. *https://doi.org/10.1038/s41467-019-12695-4*

+ Kreiner, J.M., Giacomini, D.A., Bemm, F., Waithaka, B., Regalado, J., Lanz, C., Hildebrandt, J., Sikkema, P.H., Tranel, P.J., Weigel, D., Stinchcombe, J.R., & S.I. Wright. **2019**. Multiple modes of convergent adaptation in the spread of glyphosate-resistant *Amaranthus tuberculatus*. ***Proceedings of the National Academy of Sciences, USA*** 2019 pii: 201900870. *https://doi.org/10.1073/pnas.1900870116*

+ Merondun, J., Murray, D.L., and A.B.A. Shafer. **2019**. Genome-scale sampling suggests cryptic epigenetic structuring and insular divergence in Canada lynx. ***Molecular Ecology*** 28: 3186-3196. *https://doi.org/10.1111/mec.15131*

+ Schield, D.R., Perry, B.W., Adams, R.H., Card, D.C., Jezkova, T., Pasquesi, G.I.M., Nikolakis, Z.L., Row, K., Meik, J.M., Smith, C.F., Mackessy, S.P., and T.A. Castoe. **2019**. Allopatric divergence and secondary contact with gene flow: a recurring theme in rattlesnake speciation. ***Biological Journal of the Linnaean Society*** 128: 149-169. *https://doi.org/10.1093/biolinnean/blz077*

+ Krohn, A.R., Diepeveen, E.T., Bi, K., and E.B. Rosenblum. **2019**. Local adaptation does not lead to genome-wide differentiation in lava flow lizards. ***Ecology and Evolution*** 9: 6810-6820. *https://doi.org/10.1002/ece3.5231*

+ Gray, L.N., Barley, A.J., Poe, S., Thomson, R.C., Nieto-Montes de Oca, A., and I.J. Wang. **2019**. Phylogeography of a widespread lizard complex reflects patterns of both geographic and ecological isolation. ***Molecular Ecology*** 28: 644-657. *https://doi.org/10.1111/mec.14970*

+ Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., Portik, D.M., Streicher, J.W., and S.P. Loader. **2018**. Vanishing refuge: testing the forest refuge hypothesis in coastal East Africa using genome-wide sequence data for five co-distributed amphibians. ***Molecular Ecology*** 27: 4289-4308. *https://doi.org/10.1111/mec.14862*

+ Schield, D.R., Adams, R.H., Card, D.C., Corbin, A.C., Jezkova, T., Hales, N.R., Meik, J.M., Perry, B.W., Spencer, C.L., Smith, L.L., Garcia, G.C., Bouzid, N.M., Strickland, J.L., Parkinson, C.L., Borja, M., Castañeda-Gaytán, G., Bryson, R.W., Flores-Villela, O.A., Mackessy, S.P., and T.A. Castoe. **2018**. Cryptic genetic diversity, population structure, and gene flow in the Mojave rattlesnake (Crotalus scutulatus). ***Molecular Phylogenetics and Evolution*** 127: 669-681. *https://doi.org/10.1016/j.ympev.2018.06.013*

+ Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K., Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. **2018**. Sky, sea, and forest islands: diversification in the African leaf-folding frog Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). ***Journal of Biogeography*** 45: 1781-1794. *https://doi.org/10.1111/jbi.13365*

+ Schield, D.R., Adams, R.H., Card, D.C., Perry, B.W., Pasquesi, G.M., Jezkova, T., Portik, D.M., Andrew, A.L., Spencer, C.L., Sanchez, E.E., Fujita, M.K., Mackessy, S.P., and T.A. Castoe. **2017**. Genomic patterns of divergence and admixture in a widely-distributed rattlesnake provide insight into speciation with gene flow. ***Ecology and Evolution*** 7: 3951–3966. *https://doi.org/10.1002/ece3.2996*

+ Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. **2017**. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. ***Molecular Ecology*** 26: 5245-5263. *https://doi.org/10.1111/mec.14266*

## License <a name="LIC"></a>

GNU Lesser General Public License v3.0


## Contact: <a name="CON"></a>

Daniel Portik, PhD

daniel.portik@gmail.com

