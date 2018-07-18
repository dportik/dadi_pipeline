# Two Populations Dadi Pipeline v3

**This represents my initial attempts to streamline model testing in dadi. Please use the more recent pipeline, rather than this one here. For the sake of continuity, I am keeping these scripts stored here.**

Explore demographic models capturing variation in migration rates, periods of isolation, 
and population size change associated with the divergence between two populations.

To use this workflow, you'll need a SNPs input text file to create the 2D joint site frequency
spectrum object. There are many options for converting files to the correct format, and if 
you've used STACKS to process your ddRADseq data I have written a conversion script which
can be found in the Stacks_pipeline repository. Check the dadi website for instructions on 
the basic format for this file. This pipeline is written to create folded spectra (lacking
outgroup information to polarize SNPs), but can easily be modified to created unfolded
spectrum objects.

None of the scripts are fully automated, as they are inherently specific to the input file 
being used. Sections of the scripts which require hand editing are flagged with a #**************.
Look through each script before attempting to run it to understand what needs to be modified.
At the top of each script, there are also explanations and instructions.

**Version 3 updates** (Oct 2017) 
Moved all functions out of main scripts and placed in *Optimize_Functions.py*
script for clarity. The workflow is nearly identical to V2, with less clutter. New models 
added to the model set, here is a complete list of models (new models marked with +):

1. Standard neutral model, populations never diverge.
2. Split into two populations, no migration.
3. Split into two populations, with continuous symmetric migration.
4. Split into two populations, with continuous asymmetric migration.
5. Split with continuous symmetric migration, followed by isolation.
6. Split with continuous asymmetric migration, followed by isolation.
7. Split with no gene flow, followed by period of continuous symmetrical gene flow.
8. Split with no gene flow, followed by period of continuous asymmetrical gene flow.
9. Split with no migration, then instantaneous size change with no migration.
10. Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
11. Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
12. Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.  
13. Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
14. Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
15. Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
16. +Split into two populations, with continuous symmetric migration, rate varying across two epochs.
17. +Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
18. +Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
19. +Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
20. +Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
21. +Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.

'Island' models:
For all the following models, pop2 is assumed to be the 'island' population, and pop2=nuA(s), 
pop1=nuA(1-s), where NuA = ancestral population and 's' is a fraction. Vicariance events 
involve no population size change, whereas founder event models always enforce exponential growth 
in pop2 (the island population).

22. +Island: Vicariance with no migration.
23. +Island: Vicariance with with ancient continuous asymmetric migration.
24. +Island: Vicariance with no migration, secondary contact with continuous asymmetric migration
25. +Island: Founder event with no migration.
26. +Island: Founder event with continuous symmetric migration.
27. +Island: Founder event with continuous asymmetric migration.
28. +Island: Vicariance, early unidirectional discrete admixture event (before any drift).
29. +Island: Vicariance, late unidirectional discrete admixture event (after any drift).
30. +Island: Vicariance, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.
31. +Island: Founder event, early unidirectional discrete admixture event (before any drift).
32. +Island: Founder event, late unidirectional discrete admixture event (after any drift).
33. +Island: Founder event, two epochs with unidirectional discrete admixture event occurring at beginning of the second epoch.


A new graphics file with the visual representation has been provided. The names below models 
match their label in the *Models_2D.py* allowing them to be easily matched. The example file and
all outputs have been updated to include results from several of the new models (but not all).


**Suggested Workflow**


The basic workflow involves running the scripts sequentially and in order:

**1. run dadi_2D_00_projections.py**

Determine the best down projection numbers for your data set, to maximize the number of 
segregating sites. The projection numbers represent the number of alleles, not individuals. 
The full path to the SNPs input file must be provided, along with specific population labels 
within this file, and the projection numbers. The syntax to create the 2D-JSFS will be the 
same in all subsequent scripts.

Usage:

`python dadi_2D_00_projections.py`


**2. run dadi_2D_01_first_optimizations.py**

This script will perform optimizations from multiple starting points using a 3-fold perturbed 
set of random starting values for parameters. As before, you'll need to provide the full path 
to the SNPs input file, the specific population labels, and the projection numbers. You'll
also need to provide numbers for the grid size to use. The default is to run 50 replicates
for each of the 15 models, but this can also be edited at the bottom of the script. This 
script requires the *Models_2D.py* and *Optimize_Functions.py* scripts to be in same working 
directory. This is where all the population model functions are stored. 


The output for each model is a tab-delimited  text file which can be opened and sorted to 
find the best scoring replicate. The parameter values from this run should be used as the 
starting values in the next script, *dadi_2D_02_second_optimizations.py*. The order of 
output parameters matches the input order, so they can be copied and pasted from the 
appropriate output file into the next script within the corresponding parameter list.

The outputs from here will be labeled according to model and a user selected prefix name, 
and are written to the working directory as:

`Round1_[prefix]_[model_name]_optimized.txt`

Example output file from one model:
```
Model	param_set	Replicate	log-likelihood	theta	AIC	optimized_params
Divergence with no migration	parameter set = [nu1, nu2, T]	1	-443.41	317.15	892.82	0.6795	0.7315	0.141	
Divergence with no migration	parameter set = [nu1, nu2, T]	2	-179.8	141.7	365.6	4.3798	1.6121	0.7995	
Divergence with no migration	parameter set = [nu1, nu2, T]	3	-241.28	116.99	488.56	3.1119	17.7282	1.1496	
Divergence with no migration	parameter set = [nu1, nu2, T]	4	-1763.04	54.2	3532.08	9.0388	0.2684	5.3911	
Divergence with no migration	parameter set = [nu1, nu2, T]	5	-177.49	150.28	360.98	2.7819	2.6884	0.8431	
```

Usage:

`python dadi_2D_01_first_optimizations.py`

There are potentially a lot of outputs to sort through, depending on how many models you are running.
An easy way to summarize them is to move all the *Round1_something_optimized.txt* files to a new
directory, move in to the directory, and concatenate them:

`cat *.txt > Round1_concatenated.txt`

You can open this tab-delimited summary file and sort by model name, then by AIC (smallest to
largest), to find the best scoring replicate per model. You'll use the corresponding parameter
values as the starting base parameters in the next script.

  
**3. run dadi_2D_02_second_optimizations.py**

This script will perform optimizations from multiple starting points using a 2-fold perturbed 
set of USER SELECTED starting values for parameters. These values should come from the best
replicate per model of the previous script. As before, you'll need to provide the full path 
to the SNPs input file, the specific population labels, and the projection numbers. You'll
also need to provide numbers for the grid size to use. The default is to run 50 replicates
for each of the 15 models based on the user selected starting parameters. This  script also 
requires the *Models_2D.py* and *Optimize_Functions.py* scripts to be in same working directory.

The output for each model is an identical style tab-delimited text file. Similar to the 
previous step, the parameter values from this run should be used as the starting values 
for the next and final optimization script, *dadi_2D_03_third_optimizations.py*.

The outputs from here will be labeled according to model and a user selected prefix name, 
and are written to the working directory as:

`Round2_[prefix]_[model_name]_optimized.txt`

Usage:

`python dadi_2D_02_second_optimizations.py`


**4. run dadi_2D_03_third_optimizations.py**

This script will perform optimizations from multiple starting points using a 1-fold perturbed 
set of USER SELECTED starting values for parameters. The default is to run 100 replicates
for each of the 15 models. As before, you'll need to provide the full path to the SNPs input 
file, the specific population labels, the projection numbers, and grid size. This  script 
also requires the *Models_2D.py* and *Optimize_Functions.py* scripts to be in same working 
directory. The output for each model is the same style tab-delimited text file which can be 
opened and sorted to find the best scoring replicate and do AIC comparisons, delta AIC, or 
weights for the model set. The optimized parameter values for each model should be used for 
the subsequent plotting script. 

The outputs from here will be labeled according to model and a user selected prefix name, 
and are written to the working directory as:

`Round3_[prefix]_[model_name]_optimized.txt`

Usage:

`python dadi_2D_03_third_optimizations.py`


**5. run dadi_2D_04_plotting_functions.py**

This script will simulate the model using a fixed set of input parameters. That is, no 
optimizations will occur here, and the assumption is you've searched thoroughly for the 
parameter set with the highest likelihood value. The goal of this script is to produce 
plots of the 2D-JSFS for the data and model, plus the residuals. A pdf file for each of 
the models will be written for each model to the current working directory. This script 
requires the *Models_2D.py* and *Optimize_Functions.py* scripts to be in same working 
directory.

Usage:

`python dadi_2D_04_plotting_functions.py`


**Summary**

The goal of this workflow was to automate functions of dadi to produce useful output files
resulting from many replicate runs of the model set, and to perform increasingly focused 
optimizations. The first optimization round involves 50 replicates, with the *optimize_log_fmin*
function (maxiter=20), using 3-fold perturbed default starting parameters. The second
optimization round involves 50 replicates, with the *optimize_log_fmin* function (maxiter=30), 
using 2-fold perturbed starting parameters (from user input values). The third optimization 
round involves 100 replicates, with the *optimize_log_fmin* function (maxiter=50), using 1-fold 
perturbed starting parameters (from user input values). This is done for each of the 15
models included in the Models_2D.py script. These values worked well in my analyses, but 
you may wish to change them if optimizations are not performing well or you want to be even
more thorough.

You can create your own set of models and perform similar analyses using the 
basic structure of these scripts, or use the relatively simple model set here. Adding models
will require adding additional functions to the *Models_2D.py* and *Optimize_Functions.py* 
scripts. Although this is more advanced, the model set can easily be reduced by either 
blocking out or deleting relevant sections of the first, second, and third optimization 
scripts. As each model is called at the bottom of the optimization scripts, eliminating 
them is fairly simple. For example, the partial block of code controlling the models examined at 
the bottom of the first optimization script:

```
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_divergence")
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_mig")
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig")
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig")
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_sym_mig")
``` 

can be modified to only include the first model this way:
```
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_divergence")
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_mig")
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig")
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig")
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_sym_mig")
``` 

or by deleting the other calls to the models:
```
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_divergence")
``` 

and executing the script. That is, you don't need to edit anything within the *Optimize_Round1*
function found in the *Optimize_Functions* script or eliminate anything from the 
*Models_2D.py* script to reduce the model set, just change the bottom section as above. 
Similar changes will need to be made to the following optimization and plotting scripts, 
but this will be straightforward and nearly identical.

Also, if you feel more strongly about a particular optimization routine, you can change the 
function from *optimize_log_fmin* to any of the other options available, inside the 
*Optimize_Round1*, *Optimize_Round2*, and *Optimize_Round3* functions in the 
*Optimize_Functions* script. There are many available and some may be more appropriate or 
perhaps faster for your data sets. 


**Cautions**

Optimizations sometimes fail and crash the script. The commands to call each model are at 
the bottom of the script, and the default set up is to call each model sequentially. You
can try running this as is, or create multiple verisions of this script, block out some 
of the models (use hashes, block annotation, or delete), and execute one script version 
for every core you have available. This will speed things up, and reduce potential 
crashes from optimization fails.

Be aware of other potential issues when running models, including data being masked, extrapolation
failures, etc. To troubleshoot these problems, see the manual or head to the user group: 
https://groups.google.com/forum/#!forum/dadi-user. They may not always cause the scripts 
to crash, but if the extrapolation fails the likelihoods can be entirely misleading. Make
sure you inspect the output when the replicates are running to make sure this isn't the case.

I encountered some specific issues when plotting, which terminated in errors and the plots
were not generated. Some of this required editing or duplicating the dadi code to fix. These
specific issues and workarounds are detailed at the top of the plotting script. The most
common error is very easy to fix, and this is explained in the next section.


**Plotting Errors**

Sometimes you might see this error when plotting:

`ValueError: Data has no positive values, and therefore can not be log-scaled.`

This can be fixed. This stems from an error in the *plot_2d_comp_multinom* function located
in the dadi *Plotting.py* script. You will need to change the vmin value in this function
from *None* to something between 0 and 1 (generally a super low decimal value) where it
is called in the plotting script I've written. I have changed vmin to 0.005-0.01 with good results.

To implement this change open the *dadi_2D_04_plotting_functions.py* and scroll to the 
bottom section where the code reads:
```
def plot_all(sim_model, data, outfile, model_name):
    print '{0}_{1}.pdf'.format(outfile,model_name), '\n'
    outname = '{0}_{1}.pdf'.format(outfile,model_name)
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(sim_model, data, resid_range = 3)
    fig.savefig(outname)
```

You'll want to change:

`dadi.Plotting.plot_2d_comp_multinom(model, data, resid_range = 3)`

to:

`dadi.Plotting.plot_2d_comp_multinom(model, data, resid_range = 3, vmin = 0.005)`


This should eliminate the ValueError and allow the plotting to continue. You can play around
with the vmin values to see the resulting effects on the plots. 


**Example Data**

I've provided an example SNPS input file along with the results of each optimization round
and the final spectrum plots. The scripts above will work with this input file if you
provide the correct path to its location, as the population IDs, projections, grid size,
and optimized param lists are all specific to this data set.


**Requirements**

These scripts are written for Python 2.7, and to use them requires the most up to date 
versions of the following Python modules (besides dadi):

-Numpy

-Scipy

-Matplotlib




**Citation**

If you decide to use these scripts or modify the code for your purposes, please cite:

*Portik, D.M., Leaché, A.D., Rivera, D., Blackburn, D.C., Rödel, M.-O., Barej, M.F., 
Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification 
in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology, 
Early View, doi: 10.1111/mec.14266*



**Daniel Portik**

Contact: daniel.portik@uta.edu -> danielportik@email.arizona.edu

Postdoctoral Researcher

University of Arizona

October 2017



 