**Create SFS Comparison Plots**
---------------------------------

**Purpose:**

Create figures comparing the data and model sfs for 1D, 2D, or 3D spectra.

This tool is designed to work with the Python package [dadi](https://bitbucket.org/gutenkunstlab/dadi) 
and assumes you already have the package installed. You'll need to be familiar with how dadi works, 
and some of the basic syntax for writing dadi scripts with python. A good resource for all dadi-related 
questions is the [user group](https://groups.google.com/forum/#!forum/dadi-user). Before attempting
to use these scripts, read over the user manual for dadi and try running the program with the 
example files.

**Overview:**

This is meant to be a general use script to run dadi to fit any model on an afs/jsfs with one to three populations, then create figures comparing the data and model sfs. To use this workflow, you'll need a SNPs input text file to create an allele frequency or joint site frequency spectrum object. Alternatively, you can import a frequency spectrum of your own creation, editing the script appropriately (see dadi manual). The user will have to edit information about their allele frequency spectrum, and a #************** marks lines in the *Make_Plots.py* that will have to be edited. 

The user provides a model and the previously optimized parameters for their empirical 
data. The model is fit using these parameters, and the resulting model SFS is used to
create a figure comparing the data and model SFS, including the residuals.


The *Make_Plots.py* script and *Plotting_Functions.py* script must be in the same working directory to run properly.

**Empirical Data Optimization:**

Within the *Make_Plots.py* script, let's assume you've supplied the correct information about your SNPs input file, population IDs, projection sizes, and are using the model in the script (sym_mig).

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

    #create a prefix based on the population names to label the output files
    #ex. Pop1_Pop2
    prefix = "_".join(pop_ids)

    #Make sure to define your extrapolation grid size.
    pts = [50,60,70]
    
    #Provide best optimized parameter set for empirical data.
    #These will come from previous analyses you have already completed
	emp_params = [0.1487,0.1352,0.2477,0.1877]
     
    #Fit the model using these parameters and return the model SFS.
	#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model function,
	#but everything else can stay as it is. See above for argument explanations.
	model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "sym_mig", sym_mig, emp_params, fs_folded=True)

	
**Plotting Functions:**

After the model is fit to the empirical data, the model SFS can be compared to the data SFS to create plots.
Depending on whether the SFS is 1D, 2D, or 3D, a different plotting function will be used. They all follow
the basic argument structure. I will show the 2D example below.

The 2D plotting is performed with the following function:

***Plot_2D(fs, model_fit, outfile, model_name, vmin_val=None)***
 
***Mandatory Arguments:***

+ **fs**:  spectrum object name
+ **pts**: grid size for extrapolation, list of three values
+ **outfile**:  prefix for output naming
+ **model_name**: a label help name the output files; ex. "sym_mig"

***Optional Arguments:***

+ **vmin_val**: Minimum values plotted for sfs. The default is 0.05, and to fix the common plotting error this should be changed to something between 0 and 0.05.

***Example:***

The important arguments will need to be defined in the script. Below shows how to perform
the basic plotting and also change the vmin in the second plot. 

    #Now we simply call the function with the correct arguments (notice many are the same from the
    #empirical fit).
    
    Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig")
    
    
    #Although the above function does not produce an error, let's pretend it did and
    #change the vmin value. You can see the effect on the colors in the plot. We will
    #edit the "model_name" string so the output file will be called something different.
      
    vmin_val = float(0.01)
    Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig_vmin", vmin_val = vmin_val)
    
Notice that running the ***Plot_2D*** function creates a pop-up window with the plot in it. To move along to the second plot,
simply close this window. The plot will be saved to a PDF file in the working directory automatically.

The functions are nearly identical for 1D and 3D plotting:

***Plot_3D(fs, model_fit, outfile, model_name, vmin_val=None)***

***Plot_1D(fs, model_fit, outfile, model_name)***

Note that there is no vmin_val optional argument for 1D plots, but there is for the 3D plotting.

You would use these in the same fashion:

    Plotting_Functions.Plot_1D(fs, model_fit, prefix, "something")

    Plotting_Functions.Plot_3D(fs, model_fit, prefix, "something")
    Plotting_Functions.Plot_3D(fs, model_fit, prefix, "something", vmin_val = vmin_val)


**Outputs:**

The ***Optimize_Empirical*** function will produce an output file for the empirical fit, which will be in tab-delimited format:

     Model	Replicate	log-likelihood	theta	sfs_sum	chi-squared
     sym_mig	1	-591.21	619.83	1552.44	758.21

This is based on the parameter values supplied, as no optimization routine is performed. 

Each of the plotting functions (1D, 2D, or 3D) will produce a PDF output file each
time the plotting routine is called.


**Using Folded vs. Unfolded Spectra:**

 To change whether the frequency spectrum is folded vs. unfolded requires two changes in the script. The first is where the spectrum object is created, indicated by the *polarized* argument:
 
     #Convert this dictionary into folded AFS object
     #[polarized = False] creates folded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

The above code will create a folded spectrum. When calling the empirical optimization function, this must also be indicated in the *fs_folded* argument:

     #this is from the first example:
     model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "sym_mig", sym_mig, emp_params, fs_folded=True)
     
To create an unfolded spectrum, the *polarized* and *fs_folded*  arguments in the above lines need to be changed accordingly:

     #[polarized = True] creates an unfolded spectrum object
     fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
     
     #and the optimization routine function must also be changed:
     model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "sym_mig", sym_mig, emp_params, fs_folded=False)
     
It will be clear if either argument has been misspecified because the calculation of certain statistics will cause a crash with the following error:

     ValueError: Cannot operate with a folded Spectrum and an unfolded one.

If you see this, check to make sure both relevant arguments actually agree on the spectrum being folded or unfolded.

**Helpful Notes:**

 Sometimes you may see the following error when plotting 2D or 3D, after the script crashes:
 
     "ValueError: Data has no positive values, and therefore can not be log-scaled."
 
 To deal with this, you can change the optional argument vmin_val to a value
 between 0 and 0.05 (0.05 is the default). I have used values from 0.005-0.01
 with good visual results.
 

**Test Data Set:**

In the folder labeled *Example_Data* you will find a SNPs input file that will run with the *Make_Plots.py* script.
You will only need to edit the path to the file in the script, and then the script should run normally. The 
output file of the model fit and both 2D plots from the above example are contained within the *Example_Data* folder, in a separate folder labeled *Example_Outputs*.
You should test the script using these data to ensure everything is working properly before examining your own empirical data. 


**Citation Information:**

 If you use or modify these scripts for your work, please cite the following publication.
 
+ Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. ***Molecular Ecology*** 26: 5245-5263. *https://doi.org/10.1111/mec.14266*


**Contact:**

Daniel Portik, PhD

Postdoctoral Researcher

University of Arizona

daniel.portik@gmail.com

