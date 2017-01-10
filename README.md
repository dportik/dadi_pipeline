# dadi_pipeline

This workflow is designed to work with the Python package dadi (https://bitbucket.org/gutenkunstlab/dadi) and assumes you already 
have the package installed. You'll need to be familiar with how dadi works, and some of the basic syntax for writing your own 
dadi scripts. A good resource for all dadi-related questions is the user group: https://groups.google.com/forum/#!forum/dadi-user.

To begin, you will need a SNPs input text file. Mine is generated as part of a STACKs pipeline I have written, which converts the 
haplotypes tsv file into the correct format. Depending on the program you've used to generate your SNP data file, you may or may not
have easy options to get this formatted correctly. See the example file I've provided here and the explanation on the dadi website for 
instructions on the basic format.

Any number of demographic models can be created and tested by the user, I have written my pipeline to explore 14 models of joint 
demographic history for two-population comparisons. These models are written in the TwoD_models.py script, which must be in the same
working directory when running the other scripts (they import these models) I assume outgroup data are lacking and create folded 
2D-spectra. If you have outgroup data in your SNP input file, you can easily change this option in the script. None of the scripts 
are fully automated, as they are inherently specific to the input file being used. You'll need to hand edit sections of the script to match your data set.

The general overview is:
1. run dadi_2D_initial_optimization.py
  -Initial optimizations are performed by generating 20 sets of randomly perturbed parameters and optimizing each parameter set 
  using the Nelder-Mead method (optimize_log_fmin), running each optimization algorithm for a maximum of 100 iterations. Each 
  optimized parameter set is used to simulate the 2D-JSFS, and the multinomial approach is used to estimate the log-likelihood 
  of the 2D-JSFS given the model. The output of this script: {outfile}_INITIAL_model_results.txt, is used to find the best rep 
  for each model and these parameters are manually inserted into the next script for final optimizations.

2. run dadi_2D_final_optimization.py
  -Final optimizations should use the best param sets from each model resulting from the previous script. These parameters are also
  perturbed to start optimizations, just less so than in the previous step. For each replicate, a perturbed parameter set is then optimized
  using the Nelder-Mead method (optimize_log_fmin), running each optimization algorithm for a maximum of 60 iterations. Each 
  optimized parameter set is used to simulate the 2D-JSFS, and the multinomial approach is used to estimate the log-likelihood 
  of the 2D-JSFS given the model. The output of this script: {outfile}_OPTIMZED_model_results.txt, is tab-delimited and can be opened in 
  excel. It reports rep number, ln-lik, AIC, params, etc. You can select the best replicate for each model from here.
  
3. run dadi_2D_plotting.py
  -Manually insert the best param set for each model from the previous step. There is NO OPTIMIZATION here, so the values will not change
  in the final model. The idea is simply to run each model and then plot the results, a spectrum of the data and model, plus residuals.
  Output pdf files will be created for each plot run.
  
Sections of the scripts which require hand editing are flagged with a #**************
You'll need to figure out the best projections for your data set, often requiring a down-projection. The grid choice will 
depend on the projection values, and also needs to be changed. You can choose the number of reps in model searching, and 
will need to change output file naming schemes.

The point of this script was to automate functions of dadi and produce nice looking output files. You can create your own set of models
and perform similar analyses using the basic structure of these scripts, or use the set of 14 models here. The model set here can be reduced
by deleting relevant sections of the loops, but adding models will require adding substantial blocks of code within the loops I designed.

I noticed that optimizations sometimes fail, crashing the script. If this happened repeatedly, I blocked out the code for this model and 
skipped it - sometimes they were too complex for a given data set. Sometimes the optimizations don't completely fail, but display 
optimization warnings. These aren't translated to the output files, so pay attention to the program output and make sure this isn't 
a problem - the likelihoods from these are bizarre and should not be used. 
