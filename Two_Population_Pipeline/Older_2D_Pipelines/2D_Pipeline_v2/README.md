# Two_Populations V2

**This represents my initial attempts to streamline model testing in dadi. Please use the more recent pipeline, rather than this one here. For the sake of continuity, I am keeping these scripts stored here.**

Explore 15 demographic models capturing variation in migration rates, periods of isolation, 
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
script requires the *Models_2D.py* script to be in same working directory. This is where
all the population model functions are stored. 


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


  
**3. run dadi_2D_02_second_optimizations.py**

This script will perform optimizations from multiple starting points using a 2-fold perturbed 
set of USER SELECTED starting values for parameters. These values should come from the best
replicate per model of the previous script. As before, you'll need to provide the full path 
to the SNPs input file, the specific population labels, and the projection numbers. You'll
also need to provide numbers for the grid size to use. The default is to run 50 replicates
for each of the 15 models based on the user selected starting parameters. This  script also 
requires the *Models_2D.py* script to be in same working directory.

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
also requires the *Models_2D.py* script to be in same working directory. The output for 
each model is the same style tab-delimited text file which can be opened and sorted to find 
the best scoring replicate and do AIC comparisons, delta AIC, or weights for the model set. 
The optimized parameter values for each model should be used for the subsequent plotting script. 

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
the 15 models will be written for each model to the current working directory. 

Usage:

`python dadi_2D_04_plotting_functions.py`


**Summary**

The goal of this workflow was to automate functions of dadi to produce useful output files
resulting from many replicate runs of the model set, and to perform increasingly focused 
optimizations. The first optimization round involves 50 replicates, with the *optimize_log_fmin*
function (maxiter=20), using 3-fold perturbed starting parameters (from values=1). The second
optimization round involves 50 replicates, with the *optimize_log_fmin* function (maxiter=30), 
using 2-fold perturbed starting parameters (from user input values). The third optimization 
round involves 100 replicates, with the *optimize_log_fmin* function (maxiter=50), using 1-fold 
perturbed starting parameters (from user input values). This is done for each of the 15
models included in the Models_2D.py script. 

You can create your own set of models and perform similar analyses using the 
basic structure of these scripts, or use the relatively simple model set here. Adding models
will require adding additional functions in the same style to all the scripts. The model 
set can be reduced by either blocking out or deleting relevant sections of the code. As each
model is called at the bottom of the optimization scripts, eliminating them is fairly simple.
For example the block of code controlling the models examined at the bottom of the first
optimization script:

```
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_divergence")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sym_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "asym_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_sym_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_asym_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig")

Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_mig_size")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sym_mig_size")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "asym_mig_size")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_sym_mig_size")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_asym_mig_size")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size")
``` 

can be modified to only include the first four models this way:
```
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_divergence")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sym_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "asym_mig")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_sym_mig")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_asym_mig")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig")

#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_mig_size")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sym_mig_size")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "asym_mig_size")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_sym_mig_size")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_asym_mig_size")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size")
#Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size")
``` 

or simply by deleting the other calls to the models:
```
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_divergence")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sym_mig")
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "asym_mig")
``` 

and executing the script. That is, you don't need to edit anything within the *Two_Pop_Model*
function or eliminate anything from the *Models_2D.py* script to reduce the model set, just
change the bottom section as above. Similar changes will need to be made to the following
optimization and plotting scripts, but this will be straightforward.

Also, if you feel more strongly about a particular optimization routine, you can change the 
function from *optimize_log_fmin* to any of the other options available, throughout the code.
There are many available and some may be more appropriate or faster. 


**Cautions**

Optimizations sometimes fail and crash the script. The function to call each model is at 
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
common error is very easy to fix, and this is explained in the *Plotting_Error_Fix* section.

**Example Data**

I've provided an example SNPS input file along with the results of each optimization round
and the final spectrum plots. The scripts above will work with this input file if you
provide the correct path to its location, as the population IDs, projections, grid size,
and optimized param lists are all specific to this data set. Note the standard plotting 
routine will fail for this data set, and the instructions for fixing this common issue are
available in the *Plotting_Error_Fix* section.


**Requirements**

These scripts are written for Python 2.7, and to use them requires the most up to date 
versions of the following Python modules (besides dadi):

-Numpy

-Scipy

-Matplotlib




**Citation**

If you decide to use these scripts or modify the code for your purposes, please cite:

*Portik, D.M., Leaché, A.D., Rivera, D., Blackburn, D.C., Rödel, M.-O., Barej, M.F., 
Hirschfeld, M., Burger, M., and M.K. Fujita. Evaluating mechanisms of diversification 
in a Guineo-Congolian forest frog using demographic model selection. 
In Review, Molecular Ecology.*



**Daniel Portik**

Contact: daniel.portik@uta.edu

Postdoctoral Researcher

University of Texas at Arlington

April 2017



 