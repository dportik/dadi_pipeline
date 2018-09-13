import sys
import os
import numpy
import dadi
import Plotting_Functions

'''
Usage: python Make_Plots.py

This is meant to be a general use script to fit a demographic model and create
plots of the data and model sfs. The user will have to edit information about
their allele frequency spectrum and provide a custom model. The instructions
are annotated below, with a #************** marking sections that will have to
be edited.

This script must be in the same working directory as Plotting_Functions.py, which
contains all the functions necessary for generating simulations and optimizing the model.

General workflow:
 The user provides a model and the previously optimized parameters for their empirical 
 data. The model is fit using these parameters, and the resulting model SFS is used to
 create the plot comparing the data to the model. 
 
Outputs:
 A summary file of the model fit will be created, along with a pdf of the plot.
 
Notes/Caveats:
 Sometimes you may see the following error when plotting 2D or 3D:
 
 "ValueError: Data has no positive values, and therefore can not be log-scaled."
 
 To deal with this, you can change the optional argument vmin_val to a value
 between 0 and 0.05 (0.05 is the default). I have used values from 0.005-0.01
 with good visual results.

Citations:
 If you use these scripts for your work, please cite the following publications.
 
 For the general optimization routine:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 52455263.
    doi: 10.1111/mec.14266
    
-------------------------
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated September 2018
'''

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
snps = "/Users/dan/Dropbox/dadi_inputs/General_Script/dadi_2pops_North_South_snps.txt"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["North", "South"]

#**************
#projection sizes, in ALLELES not individuals
proj = [16,32]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'



#================================================================================
# Fit the empirical data based on prior optimization results, obtain model SFS
#================================================================================
''' 
 We will use a function from the Plotting_Functions.py script:
 	Fit_Empirical(fs, pts, outfile, model_name, func, in_params, fs_folded)

Mandatory Arguments =
 	fs:  spectrum object name
 	pts: grid size for extrapolation, list of three values
 	outfile: prefix for output naming
 	model_name: a label to help name the output files; ex. "sym_mig"
 	func: access the model function from within script
 	in_params: the previously optimized parameters to use
    fs_folded: A Boolean value indicating whether the empirical fs is folded (True) or not (False).
'''


#Let's start by defining our model
def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
prefix = "_".join(pop_ids)

#**************
#Make sure to define your extrapolation grid size.
pts = [50,60,70]

#**************
#Provide best optimized parameter set for empirical data.
#These will come from previous analyses you have already completed.
emp_params = [0.1487,0.1352,0.2477,0.1877]

#Fit the model using these parameters and return the model SFS.
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model, but
#everything else can stay as it is. See above for argument explanations.
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "sym_mig", sym_mig, emp_params, fs_folded=True)



#================================================================================
# Plotting 2D
#================================================================================
'''
 We will use a function from the Plotting_Functions.py script to plot:
 	Plot_2D(fs, model_fit, outfile, model_name, vmin_val=None)

Mandatory Arguments =
	fs:  spectrum object name
    model_fit:  the model spectrum object name
 	outfile: prefix for output naming
 	model_name: a label to help name the output files; ex. "sym_mig"

Optional Arguments =
     vmin_val: Minimum values plotted for sfs, default is 0.05 and to fix the common error this should
               be changed to something between 0 and 0.05.
'''

#Now we simply call the function with the correct arguments (notice many are the same from the
#empirical fit).
Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig")


#Although the above function does not produce an error, let's pretend it did and
#change the vmin value. You can see the effect on the colors in the plot. We will
#edit the "model_name" string so the output file will be called something different.  
vmin_val = float(0.01)
Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig_vmin", vmin_val = vmin_val)



#================================================================================
# Plotting 1D or 3D spectra
#================================================================================

#Although we've used the script to plot a 2D jsfs, it can also be used to create
#plots for 1D and 3D spectra in much the same way.
 
'''
For a 1D sfs:
 We will use a function from the Plotting_Functions.py script to plot:
 	Plot_1D(fs, model_fit, outfile, model_name)

 Mandatory Arguments =
	fs:  spectrum object name
    model_fit:  the model spectrum object name
 	outfile: prefix for output naming
 	model_name: a label to help name the output files; ex. "sym_mig"

Example 1D plot call (will not work with this 2D sfs!):

Plotting_Functions.Plot_1D(fs, model_fit, prefix, "something")



For a 3D jsfs:
 We will use a function from the Plotting_Functions.py script to plot:
 	Plot_3D(fs, model_fit, outfile, model_name, vmin_val=None)

 Mandatory Arguments =
	fs:  spectrum object name
    model_fit:  the model spectrum object name
 	outfile: prefix for output naming
 	model_name: a label to help name the output files; ex. "sym_mig"

 Optional Arguments =
     vmin_val: Minimum values plotted for sfs, default is 0.05 and to fix the common error this should
               be changed to something between 0 and 0.05.

#Example 3D plot call (will not work with this 2D sfs!):

Plotting_Functions.Plot_3D(fs, model_fit, prefix, "something")
'''

