# Two Population Demographic Modeling Pipeline using dadi

---------------------------------

Explore demographic models capturing variation in migration rates, periods of isolation, and population size change associated with the divergence between two populations.

This is a modified version of the `dadi_Run_Optimizations.py` script in which we run optimizations for 2D comparisons for a large set of models that have been made available as part of published works. These models are stored in the `Models_2D.py` script, and will be called directly here. This 2D workflow is a component of the `dadi_pipeline` package.

## General Overview:

There are a variety of 2D models available that can be applied to your data set. The commands for using all of the available models are in the `dadi_Run_2D_Set.py `script, and these can be modified to only analyze a subset of the models. **In general, you should NOT try to run all of the models. Some models were created for specific projects and only make sense in particular biological contexts. Look at the available models in [Models_2D.pdf](https://github.com/dportik/dadi_pipeline/blob/master/Two_Population_Pipeline/Models_2D.pdf) and decide which ones make sense for your project and for your system.**

The optimization routine runs a user-defined number of rounds, each with a user-defined or default number of replicates. The starting parameters are initially random, but after each round is complete the parameters of the best scoring replicate from that round are used to generate perturbed starting parameters for the replicates of the subsequent round. The arguments controlling steps of the optimization algorithm (maxiter) and perturbation of starting parameters (fold) can be supplied by the user for more control across rounds. The user can also supply their own set of initial parameters, or set custom bounds on the parameters (upper_bound and lower_bound) to meet specific model needs. Because this script will generate many output files for all the models included to analyze, the `Summarize_Outputs.py` script can be used to find the best scoring replicate from each model, which will be written to a summary output file.

To use this workflow, you'll need a SNPs input text file to create the 2D joint site frequency spectrum object. Check the dadi website for instructions on the basic format for this file. This pipeline is written to create folded spectra (lacking outgroup information to polarize SNPs), but can easily be modified to created unfolded spectrum objects (see Caveats section).

The user will have to edit information about their allele frequency spectrum, and a #************** marks lines in the `dadi_Run_2D_Set.py` that will have to be edited. 

The `dadi_Run_2D_Set.py`, `Optimize_Functions.py` (from the [main](https://github.com/dportik/dadi_pipeline) repository), and `Models_2D.py` scripts must all be in the same working directory for `dadi_Run_2D_Set.py` to run properly.

## Available Two Population (2D) Models:

Here is a running list of the models currently available. The name of the model function in the `Models_2D.py` script is given, along with a brief description, and a corresponding visual representation of the model is provided in the file [Models_2D.pdf](https://github.com/dportik/dadi_pipeline/blob/master/Two_Population_Pipeline/Models_2D.pdf). If you use these models and scripts for your work, please provide proper citations (provided below).

**Diversification Model Set**

+ *no_mig*: Split into two populations, no migration.
+ *sym_mig*: Split into two populations, with continuous symmetric migration.
+ *asym_mig*: Split into two populations, with continuous asymmetric migration.
+ *anc_sym_mig*: Split with continuous symmetric migration, followed by isolation.
+ *anc_asym_mig*: Split with continuous asymmetric migration, followed by isolation.
+ *sec_contact_sym_mig*: Split with no gene flow, followed by period of continuous symmetrical gene flow.
+ *sec_contact_asym_mig*: Split with no gene flow, followed by period of continuous asymmetrical gene flow.
+ *no_mig_size*: Split with no migration, then instantaneous size change with no migration.
+ *sym_mig_size*: Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
+ *asym_mig_size*: Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
+ *anc_sym_mig_size*: Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.  
+ *anc_asym_mig_size*: Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
+ *sec_contact_sym_mig_size*: Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
+ *sec_contact_asym_mig_size*: Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
+ *sym_mig_twoepoch*: Split into two populations, with continuous symmetric migration, rate varying across two epochs.
+ *asym_mig_twoepoch*: Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
+ *sec_contact_sym_mig_three_epoch*: Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
+ *sec_contact_asym_mig_three_epoch*: Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
+ *sec_contact_sym_mig_size_three_epoch*: Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
+ *sec_contact_asym_mig_size_three_epoch*: Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.

**Island Model Set**

The island model set was recently revised as part of release v3.1.1.

For all the following models, the 'island' population is population 2 and the mainland population is population 1. Their sizes are based on the parameter s, where 0 < s < 1. Here, s represents the fraction of Nref that enters the island population, with pop2 (island) = s, and pop1 (mainland) = 1-s. When these models are called in the pipeline the upper bound on s is set to 0.5, such that the island population (pop2) must be composed of <50% of the ancestral population. This is to reflect the idea that the size of the island population should be smaller (at least initially). 

The 'vicariance' models involve no population size change, whereas the 'founder' event models enforce exponential growth in the island population. For vicariance models, the size parameters of pop1 (mainland) and pop2 (island) are 1-s and s, respectively. For the founder effect models, there is exponential growth in the island population (s), such that it achieves the final size of nu2. Therefore, for the founder event models the size parameters of pop1 (mainland) and pop2 (island) are 1-s and nu2, respectively. Please see the `Models_2D.py` script for additional details, and see the [Models_2D.pdf](https://github.com/dportik/dadi_pipeline/blob/master/Two_Population_Pipeline/Models_2D.pdf) file for a visual representation.


+ *vic_no_mig*: Vicariance with no migration.
+ *vic_anc_sym_mig*: Vicariance with ancient symmetric migration.
+ *vic_anc_asym_mig*: Vicariance with ancient asymmetric migration.
+ *vic_sec_contact_sym_mig*: Vicariance with no migration, secondary contact with symmetric migration.
+ *vic_sec_contact_asym_mig*: Vicariance with no migration, secondary contact with asymmetric migration.
+ *founder_nomig*: Founder event with no migration.
+ *founder_sym*: Founder event with symmetric migration.
+ *founder_asym*: Founder event with asymmetric migration.
+ *vic_no_mig_admix_early*: Vicariance, early unidirectional discrete admixture event (before drift).
+ *vic_no_mig_admix_late*: Vicariance, late unidirectional discrete admixture event (after drift).
+ *vic_two_epoch_admix*: Vicariance, middle unidirectional discrete admixture event (between two drift events).
+ *founder_nomig_admix_early*: Founder event with no migration, early unidirectional discrete admixture event (before drift).
+ *founder_nomig_admix_late*: Founder event with no migration, late unidirectional discrete admixture event (after drift).
+ *founder_nomig_admix_two_epoch*: Founder event, middle unidirectional discrete admixture event (between two drift events).


## Optimization Routine Settings:

To control the optimization routine is relatively easy, and the arguments are located in the script before all the
commands to call the models:

    #Set the number of rounds here
    rounds = 4

    #define the lists for optional arguments
    #you can change these to alter the settings of the optimization routine
    reps = [10,20,30,40]
    maxiters = [3,5,10,15]
    folds = [3,2,2,1]

The settings I've written here will provide four rounds of increasingly focused optimizations, and the values
of the arguments across rounds are summarized in the table below. 


| Argument | Round 1 | Round 2  | Round 3 | Round 4 |
| ------ |------:| -----:| -----:| -----:|
| reps    | 10 | 20 | 30 | 40 | 
| maxiter | 3 |  5  | 10 | 15 |
| fold |  3 |  2   | 2 | 1 |


If you change the number of rounds, you have to change the list length of the reps, maxiters, and folds arguments to match.

The default optimizer used is the Nelder-Mead method (`Inference.py` function `optimize_log_fmin`). This can be changed by supplying the optional optimizer argument 
with any of the following choices:

+ log - Optimize log(params) to fit model to data using the BFGS method.
+ log_lbfgsb - Optimize log(params) to fit model to data using the L-BFGS-B method.
+ log_fmin - Optimize log(params) to fit model to data using Nelder-Mead. **This is the default option.**
+ log_powell - Optimize log(params) to fit model to data using Powell's method.


## Why Perform Multiple Rounds of Optimizations?

When fitting demographic models, it is important to perform multiple runs and ensure that final optimizations are converging on a similar log-likelihood score. In this workflow, the starting parameters used for all replicates in first round are random. After each round is completed, the parameters of the best scoring replicate from the previous round are then used to generate perturbed starting parameters for the replicates of the subsequent round. This optimization strategy of focusing the parameter search space improves the log-likelihood scores and generally results in convergence in the final round. 

Below is a summary of the log-likelihood scores obtained using the default four-round optimization settings present in the 2D pipeline. This analysis was conducted for a particular model (nomig, the simplest 2D model) using the example data provided. You can clearly see the improvement in log-likelihood scores and decrease in variation among replicates as the optimization rounds progress. 

![Rounds](https://github.com/dportik/dadi_pipeline/blob/master/Example_Data/NoMig_Zoom.png)

## Performing Multiple Independent Optimization Routines:

It is also a good idea to optimize from multiple starting points, that is to run the complete optimization routine multiple times. If the true log-likelihood is being reached, the same score (or close) should be present in the final round across multiple independent executions of the optimization routine.

Executing the optimization routine multiple times can be accomplished by writing loops in the script, or by running the script multiple times. Here is an example of a custom loop:

    for i in range(1,6):
        Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")

The above loop will run the optimization routine to completion five separate times. 
Note that when you use the range argument in python it will go up to, but not include, the final number.
That's why I have written a range of 1-6 to perform this 5 times. Also note the indentation required when
writing loops in python, the proper indentation required to create a loop is four spaces!

As an alternative to custom loops, you can also just execute the script multiple times, and the same output files will simply be added to similar to what occurs with the loop. In both cases,
the outputs across models can be easily summarized as described below. 


## Modifying the Model Set to Analyze:

The model set can easily be reduced by either blocking out or deleting relevant sections. 
Let's say you no longer wish to include the first two models in the example below. There are three easy options for doing this.

    # Split into two populations, no migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")

    # Split into two populations, with continuous symmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")

    # Split into two populations, with continuous asymmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

    # Split with continuous symmetric migration, followed by isolation.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")

**Option 1:** You can hash out the Optimize_Functions.Optimize_Routine function for those models, using the # character.

    # Split into two populations, no migration.
    #Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")

    # Split into two populations, with continuous symmetric migration.
    #Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")

    # Split into two populations, with continuous asymmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

    # Split with continuous symmetric migration, followed by isolation.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")

**Option 2:** You can block out this section of the Optimize_Functions.Optimize_Routine function for those models using triple quotes. Anything contained within the set of ''' will be ignored. 

    '''
    # Split into two populations, no migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")

    # Split into two populations, with continuous symmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")
    '''

    # Split into two populations, with continuous asymmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

    # Split with continuous symmetric migration, followed by isolation.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")

**Option 3:** You can simply delete the lines!

    # Split into two populations, with continuous asymmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

    # Split with continuous symmetric migration, followed by isolation.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")


The model set can be added to by inserting your model in the *Models_2D.py* script, then adding an appropriate
call for the Optimize_Functions.Optimize_Routine function, similar to the other models. However,
the simplest and easiest way to analyze a new or custom model is to use the flexible *dadi_Run_Optimizations.py* script,
changing the optional arguments to match the settings used here for the four rounds. 


## Outputs:

 For each model run, there will be a log file showing the optimization steps per replicate and a summary file that has all the important information. 
 
Here is an example of the output from a summary file, which will be in tab-delimited format:

     Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T) 
     sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
     sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
     sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
     sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
     sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688

## Summarizing Outputs Across Models:

 The information contained in each model summary file can be quickly sorted and extracted using the `Summarize_Outputs.py` script. 
 
 The usage of this script is straightforward, and only requires placing the path to the directory containing the outputs on the command line. Here is an example of the usage, where
 the path to a folder called *My_Output_files* is specified:
 
     python Summarize_Outputs.py /Users/dan/dadi_pipeline/TwoPopulationComparisons/My_Output_files
 
 Here, the information for the best-scoring replicate for each model will be compiled and written to a tab-delimited output file called `Results_Summary_Short.txt`. Here is an example of the contents:
 

    Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params
    asym_mig_size	Round_4_Replicate_25	-840.1	1696.2	589.73	254.86	2.9196,6.337,1.5024,0.103,0.0673,0.0419,5.0356,0.0374	
    anc_asym_mig_size	Round_4_Replicate_22	-842.02	1700.04	589.99	570.99	1.0406,4.3882,1.6501,0.0587,0.2556,0.0109,1.5419,0.0249	
    sym_mig_size	Round_4_Replicate_4	-845.24	1704.48	589.66	294.33	2.2229,9.7653,2.6256,0.178,0.082,4.1593,0.0897	
	sec_contact_sym_mig_size	Round_4_Replicate_39	-874.78	1763.56	714.24	546.93	1.62,4.4826,0.9601,0.0897,0.2029,1.4814,0.0389	

This file can be sorted for the purpose of performing model comparisons. However, it is a good idea to inspect not only the top replicate, but several of the highest scoring replicates to check for consistency. For this reason, a second output file is created that is called `Results_Summary_Extended.txt`. Rather than simply containing the top-scoring replicate, it will contain the top five replicates for each model. Here is an example of the contents:

    Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params
    anc_asym_mig	Round_4_Replicate_19	-1131.52	2275.04	1500.39	267.16	2.8073,1.4011,0.0808,0.607,8.8704,0.7333	
    anc_asym_mig	Round_4_Replicate_29	-1132.42	2276.84	1458.8	173.39	4.2429,2.3167,0.0661,0.2715,11.0258,0.9718	
    anc_asym_mig	Round_4_Replicate_4	-1133.08	2278.16	1467.05	362.0	2.0616,1.1283,0.1343,0.5381,4.2699,0.4544	
    anc_asym_mig	Round_4_Replicate_6	-1133.57	2279.14	1469.86	337.43	2.2279,1.0477,0.0629,0.8164,17.85,0.5144	
    anc_asym_mig	Round_4_Replicate_17	-1134.03	2280.06	1483.91	189.57	3.9022,2.2126,0.077,0.2782,9.057,0.934	
    anc_asym_mig_size	Round_4_Replicate_22	-842.02	1700.04	589.99	570.99	1.0406,4.3882,1.6501,0.0587,0.2556,0.0109,1.5419,0.0249	
    anc_asym_mig_size	Round_4_Replicate_12	-853.1	1722.2	620.43	443.74	1.2273,6.594,4.1622,0.0728,0.2343,0.012,2.1966,0.032	
    anc_asym_mig_size	Round_4_Replicate_45	-865.39	1746.78	658.24	667.7	0.923,3.0538,1.3797,0.0971,0.3135,0.0125,1.1266,0.0404	
    anc_asym_mig_size	Round_4_Replicate_3	-867.69	1751.38	689.86	519.07	1.3274,4.7944,1.556,0.1002,0.1578,0.0134,1.5759,0.0386	
    anc_asym_mig_size	Round_4_Replicate_19	-869.06	1754.12	632.47	508.51	0.8137,6.8294,2.0211,0.1258,0.4475,0.0172,1.9328,0.0614	

This information can be used for model comparisons using AIC, etc., but see below for an important note about comparisons.

## Using Folded vs. Unfolded Spectra:

 To change whether the frequency spectrum is folded vs. unfolded requires two changes in the script. The first is where the spectrum object is created, indicated by the *polarized* argument:
 
     #Convert this dictionary into folded AFS object
     #[polarized = False] creates folded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

The above code will create a folded spectrum. When calling the optimization function for each model, this is indicated by the *fs_folded* argument:

     #Note that in the script fs_folded is assigned a variable and referred to in the optimization functions:
     
     #**************
     #Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
     fs_folded = True
     Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")
     
To create an unfolded spectrum, the *polarized* and *fs_folded*  arguments in the above lines need to be changed accordingly:

     #[polarized = True] creates an unfolded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
     
     #and the optimization routine function must also be changed:
     #Change this variable to False to set the argument fs_folded in all the model optimizations
     fs_folded = False
     Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")
     
It will be clear if either argument has been misspecified because the calculation of certain statistics will cause a crash with the following error:

     ValueError: Cannot operate with a folded Spectrum and an unfolded one.

If you see this, check to make sure both relevant arguments actually agree on the spectrum being folded or unfolded.

## Caveats:

 The likelihood and AIC returned represent the true likelihood only if the SNPs are unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this is true, but if SNPs are linked across loci then the likelihood is actually a composite likelihood and using something like AIC is no longer appropriate for model comparisons. See the discussion group for more information on this subject. 


## Contributing to the 2D Model Set:
 
 If you would like to contribute your custom 2D models to the model set, please email me at daniel.portik@gmail.com. If your models have been or will be published, I will provide explicit citation information for the contributed models (as below). 
 
## Citation Information:

### How to cite `dadi_pipeline`:

This demographic modeling pipeline was built with a novel multi-round optimization routine, it includes many original models, and it generates custom output files. Because of these important features, `dadi_pipeline` is best considered as an additional package. It was published as part of [Portik et al. (2017)](https://doi.org/10.1111/mec.14266). If you have used `dadi_pipeline` to run your analyses, please indicate so in your publication. Here is an example of how to cite this workflow:

> To explore alternative demographic models, we used the diffusion approximation method of dadi (Gutenkunst et al. 2009) to analyze joint site frequency spectra. We fit 15 demographic models using dadi_pipeline v3.1.4 (Portik et al. 2017).

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


## Contact:

Daniel Portik, PhD

Postdoctoral Researcher

California Academy of Sciences

daniel.portik@gmail.com

