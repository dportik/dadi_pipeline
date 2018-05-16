# older 2D dadi pipeline

**This represents my initial attempt to streamline model testing in dadi. Please use the more recent pipeline, rather than this one here. For the sake of continuity, I am keeping these scripts stored here.**

This workflow is designed to work with the Python package dadi (https://bitbucket.org/gutenkunstlab/dadi) and assumes you already 
have the package installed. You'll need to be familiar with how dadi works, and some of the basic syntax for writing your own 
dadi scripts. A good resource for all dadi-related questions is the user group: https://groups.google.com/forum/#!forum/dadi-user.
Before using these scripts, read over the user manual for dadi and try running the program with their example files.

To begin, you will need a SNPs input text file. Mine is generated as part of a STACKs pipeline I have written, which converts the 
haplotypes tsv file into the correct format. Depending on the program you've used to generate your SNP data file, you may or may not
have easy options to get this formatted correctly. See the example file I've provided here and the explanation on the dadi website for 
instructions on the basic format.

Any number of demographic models can be created and tested by the user, I have written my pipeline to explore 14 models of joint 
demographic history for two-population comparisons. These models are written in the TwoD_models.py script, which must be in the same
working directory when running the other scripts (they import these models). A visual summary of these models is included as a png file.

The general overview is

1. run dadi_2D_initial_optimization.py
  -Initial optimizations are performed by generating 20 sets of randomly perturbed parameters and optimizing each parameter set 
  using the Nelder-Mead method (optimize_log_fmin), running each optimization algorithm for a maximum of 100 iterations. Each 
  optimized parameter set is used to simulate the 2D-JSFS, and the multinomial approach is used to estimate the log-likelihood 
  of the 2D-JSFS given the model. The output of this script: {outfile}_INITIAL_model_results.txt, is used to find the best rep 
  for each model and these parameters are manually inserted into the next script for final optimizations.

2. run dadi_2D_final_optimization.py
  -Final optimizations should use the best param sets from each model resulting from the previous script. These parameters are also
  perturbed to start optimizations, just less so than in the previous step. For each replicate, a perturbed parameter set is then
  optimized using the Nelder-Mead method (optimize_log_fmin), running each optimization algorithm for a maximum of 60 iterations. Each 
  optimized parameter set is used to simulate the 2D-JSFS, and the multinomial approach is used to estimate the log-likelihood 
  of the 2D-JSFS given the model. The output of this script: {outfile}_OPTIMIZED_model_results.txt, is tab-delimited and can be opened in 
  excel. It reports rep number, ln-lik, AIC, params, etc. You can select the best replicate for each model from here.
  
3. run dadi_2D_plotting.py
  -Manually insert the best param set for each model from the previous step. There is NO OPTIMIZATION here, so the values will not change
  in the final model. The idea is simply to run each model and then plot the results, a spectrum of the data and model, plus residuals.
  Output pdf files will be created for each plot run.

I assume outgroup data are lacking and create folded  2D-spectra. If you have outgroup data in your SNP input file, you can easily change
this option in the script. None of the scripts are fully automated, as they are inherently specific to the input file being used.  You'll
need to figure out the best projections for your data set, often requiring a down-projection. The grid choice will depend on the
projection values, and also needs to be changed. You can choose the number of reps in model searching, and will need to change output 
file naming schemes. You'll need to hand edit sections of the script to match your data set.

Sections of the scripts which require hand editing are flagged with a #**************

The point of this workflow was to automate functions of dadi and produce output files that summarize a large number of replicates run for 
many different models. You can create your own set of models and perform similar analyses using the basic structure of these scripts, or
use the set of 14 models here. The model set here can be reduced by either blocking out or deleting relevant sections of the loops, but
adding models will require adding code within the loops I designed. This should be straightforward to imitate the style of the code, and
change relevant names, output text, parameter numbers, and AIC penalty depending on the model. I'm sure a more generic function could
be used to load in models and parameters and run optimizations, but I found this code intuitive and easy to follow so I could troubleshoot
the many issues I had along the way.

I noticed that optimizations sometimes fail, crashing the script. If this happened repeatedly, I blocked out the code for this model and 
skipped it - sometimes they were too complex for a given data set, particularly with a low number of segregating sites or low number of 
individuals. Sometimes the optimizations don't completely fail, but display severe warnings (that extrapolation may have failed). These
aren't translated to the output files, so pay attention to the program output and make sure this isn't a problem - the likelihoods from
these can be bizarre and should not be used in your model comparisons! 

The scripts here are all tailored to run the example file I've provided. You should be able to perform each step by only editing the 
path to the example file (and assuming you have the dadi python module already installed!).

Contact: daniel.portik@uta.edu

# Daniel Portik

Postdoctoral Researcher

University of Texas at Arlington

December 2016
