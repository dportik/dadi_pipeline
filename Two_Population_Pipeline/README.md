**Two Populations dadi Pipeline**
---------------------------------

Explore demographic models capturing variation in migration rates, periods of isolation, and population size change associated with the divergence between two populations.

This is a modified version of the *dadi_Run_Optimizations.py* script in which we run optimizations for 2D comparisons for a large set of models that have been made available as part of published works. These models are stored in the *Models_2D.py* script, and will be called directly here. 

**General Overview:**

There are many 2D models available that can be applied to your data set. The commands for using all of the available models are in the *dadi_Run_2D_Set.py* script, and these can be modified to only analyze a subset of the models. The optimization routine runs a user-defined number of rounds, each with a user-defined or default number of replicates. The starting parameters are initially random, but after each round is complete the parameters of the best scoring replicate from that round are used to generate perturbed starting parameters for the replicates of the subsequent round. The arguments controlling steps of the optimization algorithm (maxiter) and perturbation of starting parameters (fold) can be supplied by the user for more control across rounds. The user can also supply their own set of initial parameters, or set custom bounds on the parameters (upper_bound and lower_bound) to meet specific model needs. Because this script will generate many output files for all the models included to analyze, the *Summarize_Outputs.py* script can be used to find the best scoring replicate from each model, which will be written to a summary output file.

To use this workflow, you'll need a SNPs input text file to create the 2D joint site frequency spectrum object. Check the dadi website for instructions on the basic format for this file. This pipeline is written to create folded spectra (lacking outgroup information to polarize SNPs), but can easily be modified to created unfolded spectrum objects (see Caveats section).

The user will have to edit information about their allele frequency spectrum, and a #************** marks lines in the *dadi_Run_2D_Set.py* that will have to be edited. 

The *dadi_Run_2D_Set.py*, *Optimize_Functions.py* (from the [main](https://github.com/dportik/dadi_pipeline) repository), and *Models_2D.py* scripts must all be in the same working directory for *dadi_Run_2D_Set.py* to run properly.

**Available 2D Models**

Here is a running list of the models currently available. The name of the model function in the *Models_2D.py* script is given, along with a brief description, and a corresponding visual representation of the model is provided in the file *Models_2D.pdf*. If you use these models and scripts for your work, please provide proper citations (provided below).

***Diversification Model Set:***

1. *no_mig*: Split into two populations, no migration.
2. *sym_mig*: Split into two populations, with continuous symmetric migration.
3. *asym_mig*: Split into two populations, with continuous asymmetric migration.
4. *anc_sym_mig*: Split with continuous symmetric migration, followed by isolation.
5. *anc_asym_mig*: Split with continuous asymmetric migration, followed by isolation.
6. *sec_contact_sym_mig*: Split with no gene flow, followed by period of continuous symmetrical gene flow.
7. *sec_contact_asym_mig*: Split with no gene flow, followed by period of continuous asymmetrical gene flow.
8. *no_mig_size*: Split with no migration, then instantaneous size change with no migration.
9. *sym_mig_size*: Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
10. *asym_mig_size*: Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
11. *anc_sym_mig_size*: Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.  
12. *anc_asym_mig_size*: Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
13. *sec_contact_sym_mig_size*: Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
14. *sec_contact_asym_mig_size*: Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
15. *sym_mig_twoepoch*: Split into two populations, with continuous symmetric migration, rate varying across two epochs.
16. *asym_mig_twoepoch*: Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
17. *sec_contact_sym_mig_three_epoch*: Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
18. *sec_contact_asym_mig_three_epoch*: Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
19. *sec_contact_sym_mig_size_three_epoch*: Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
20. *sec_contact_asym_mig_size_three_epoch*: Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.

***Island Model Set:***

For all the following models, pop2 is assumed to be the 'island' population, and pop2=nuA(s), pop1=nuA(1-s), where NuA = ancestral population and 's' is a fraction. Vicariance events involve no population size change, whereas founder event models always enforce exponential growth in pop2 (the island population). Discrete admixture events are included as a comparison to intervals of continuous migration, and are represented by parameter 'f', in which a fraction f of the mainland population is instantaneously present in the post-admixture island population. Because of the way these models are constructed, you should not analyze these in combination with the Diversification Model Set for the purpose of performing model selection.

21. *vic_no_mig*: Vicariance with no migration.
22. *vic_anc_asym_mig*: Vicariance with with ancient continuous asymmetric migration.
23. *vic_sec_contact_asym_mig*: Vicariance with no migration, secondary contact with continuous asymmetric migration
24. *founder_nomig*: Founder event with no migration.
25. *founder_sym*: Founder event with continuous symmetric migration.
26. *founder_asym*: Founder event with continuous asymmetric migration.
27. *vic_no_mig_admix_early*: Vicariance, early unidirectional discrete admixture event (before any drift).
28. *vic_no_mig_admix_late*: Vicariance, late unidirectional discrete admixture event (after any drift).
29. *vic_two_epoch_admix*: Vicariance, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
30. *founder_nomig_admix_early*: Founder event, early unidirectional discrete admixture event (before any drift).
31. *founder_nomig_admix_late*: Founder event, late unidirectional discrete admixture event (after any drift).
32. *founder_nomig_admix_two_epoch*: Founder event, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.


**Optimization Settings:**

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
| ***reps***    | 10 | 20 | 30 | 40 | 
| ***maxiter*** | 3 |  5  | 10 | 15 |
| ***fold*** |  3 |  2   | 2 | 1 |


If you change the number of rounds, you have to change the list length of the reps, maxiters, and folds arguments to match.

It is also a good idea to optimize from multiple starting points, that is to run the above configuration multiple times.
This can be accomplished by writing loops or by running the main script multiple times. Here is an example of a custom loop:

    for i in range(1,6):
        Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")

The above loop will run the optimization routine to completion five separate times. 
Note that when you use the range argument in python it will go up to, but not include, the final number.
That's why I have written a range of 1-6 to perform this 5 times. Also note the indentation required when
writing loops in python.

As an alternative to custom loops, you can also just execute the script multiple times, and the same output files will simply be added to similar to what occurs with the loop. In both cases,
the outputs across models can be easily summarized as described below. 

**Modifying the Model Set to Analyze:**

The model set can easily be reduced by either blocking out or deleting relevant sections. 
Let's say you no longer wish to include the first two models in the example below:

    # Split into two populations, no migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")

    # Split into two populations, with continuous symmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")

    # Split into two populations, with continuous asymmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

    # Split with continuous symmetric migration, followed by isolation.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")

You can hash out the Optimize_Functions.Optimize_Routine function for those models as here:

    # Split into two populations, no migration.
    #Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")

    # Split into two populations, with continuous symmetric migration.
    #Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")

    # Split into two populations, with continuous asymmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

    # Split with continuous symmetric migration, followed by isolation.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")

Or you can block out this section of the Optimize_Functions.Optimize_Routine function for those models using triple quotes:


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

Anything contained within the set of ''' will be ignored. 

Finally, you can simply delete these lines:

    # Split into two populations, with continuous asymmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

    # Split with continuous symmetric migration, followed by isolation.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")


The model set can be added to by inserting your model in the *Models_2D.py* script, then adding an appropriate
call for the Optimize_Functions.Optimize_Routine function, similar to the other models. However,
the simplest and easiest way to analyze a custom model is to use the flexible *dadi_Run_Optimizations.py* script,
changing the optional arguments to match the settings used here. 


**Outputs:**

 For each model run, there will be a log file showing the optimization steps per replicate and a summary file that has all the important information. 
 
Here is an example of the output from a summary file, which will be in tab-delimited format:

     Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T) 
     sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
     sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
     sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
     sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
     sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688

**Summarizing Outputs Across Models:**

 The information contained in each model summary file can be quickly sorted and extracted using the *Summarize_Outputs.py* script. 
 
 The usage of this script is straightforward, and only requires placing the path to the directory containing the outputs on the command line. Here is an example of the usage, where
 the path to a folder called *My_Output_files* is specified:
 
     python Summarize_Outputs.py /Users/dan/dadi_pipeline/TwoPopulationComparisons/My_Output_files
 
 Here, the information for the best-scoring replicate for each model will be compiled and written to a tab-delimited output file called *Results_Summary_Short.txt*. Here is an example of the contents:
 

    Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params
    asym_mig_size	Round_4_Replicate_25	-840.1	1696.2	589.73	254.86	2.9196,6.337,1.5024,0.103,0.0673,0.0419,5.0356,0.0374	
    anc_asym_mig_size	Round_4_Replicate_22	-842.02	1700.04	589.99	570.99	1.0406,4.3882,1.6501,0.0587,0.2556,0.0109,1.5419,0.0249	
    sym_mig_size	Round_4_Replicate_4	-845.24	1704.48	589.66	294.33	2.2229,9.7653,2.6256,0.178,0.082,4.1593,0.0897	
	sec_contact_sym_mig_size	Round_4_Replicate_39	-874.78	1763.56	714.24	546.93	1.62,4.4826,0.9601,0.0897,0.2029,1.4814,0.0389	

This file can be sorted for the purpose of performing model comparisons. However, it is a good idea to inspect not only the top replicate, but several of the highest scoring replicates to check for consistency. For this reason, a second output file is created that is called *Results_Summary_Extended.txt*. Rather than simply containing the top-scoring replicate, it will contain the top five replicates for each model. Here is an example of the contents:

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

**Using Folded vs. Unfolded Spectra:**

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

**Caveats:**

 The likelihood and AIC returned represent the true likelihood only if the SNPs are unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this is true, but if SNPs are linked across loci then the likelihood is actually a composite likelihood and using something like AIC is no longer appropriate for model comparisons. See the discussion group for more information on this subject. 


**Contributing to the 2D Model Set:**
 
 If you would like to contribute your custom 2D models to the model set, please email me at daniel.portik@gmail.com. If your models have been or will be published, I will provide explicit citation information for the contributed models (see below). 
 
**Citation Information:**

***Using the pipeline.***
The scripts involved with this pipeline were originally published as part of the following work:

+ Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. ***Molecular Ecology*** 26: 5245-5263. *https://doi.org/10.1111/mec.14266*


***Using specific 2D models.***
The model sets were developed through several publications, and we are making our models available for replication and future usage.

Models 1-14 were written for:

+ Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. ***Molecular Ecology*** 26: 5245-5263. *https://doi.org/10.1111/mec.14266*

Models 15-20 and 21-32 were written for:

+ Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K., Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. 2018. Sky, sea, and forest islands: diversification in the African leaf-folding frog Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). ***Journal of Biogeography*** 45: 1781-1794. *https://doi.org/10.1111/jbi.13365*


**Contact**

Daniel Portik, PhD

Postdoctoral Researcher

University of Arizona

daniel.portik@gmail.com

