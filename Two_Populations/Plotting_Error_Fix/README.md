# Plotting Error


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



 