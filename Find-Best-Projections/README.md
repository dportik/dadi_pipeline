# Test Projection Sizes for JSFS with dadi

---------------------------------

## Purpose:

Quickly generate combinations of projection sizes to determine the resulting numbers of segregating sites. This workflow is a component of the `dadi_pipeline` package.

This tool is designed to work with the Python package [dadi](https://bitbucket.org/gutenkunstlab/dadi) 
and assumes you already have the package installed. You'll need to be familiar with how dadi works, 
and some of the basic syntax for writing dadi scripts with python. A good resource for all dadi-related 
questions is the [user group](https://groups.google.com/forum/#!forum/dadi-user). Before attempting
to use these scripts, read over the user manual for dadi and try running the program with the 
example files.

## Overview:

This is meant to be a general use script to generate combinations of projection sizes to determine the resulting numbers of segregating sites. To use this workflow, you'll need a SNPs input text file to create an allele frequency or joint site frequency spectrum object. Alternatively, you can import a frequency spectrum of your own creation, editing the script appropriately (see dadi manual). The user will have to edit information about their allele frequency spectrum, and a #************** marks lines in the script that will have to be edited. 

The user provides an input SNPs file along with the relevant information (population IDs, folded vs. unfolded). Finally, a maximum projection size is set and a minimum fraction required for the projection sizes. The script will then generate all permutations of projection size combinations and write the results to stdout. 

The `dadi-test-projections.py` script can be run independently and does not require any other scripts from the `dadi_pipeline`.

## What to Edit:

Within the script, you will need to edit information about your input file and spectrum characteristics, similar to all other scripts in `dadi_pipeline`. In addition, you will need to set the minimum fraction for the projection size (default fraction is 0.5, or 50%). Using this default value of 0.5, if your max projection sizes are 40 and 50, then this would set a lower limit of 20 and 25 for the projection size tests, respectively. 


## Outputs:

The combinations of projections sizes and resulting numbers of segregating sites will be written to stdout. The output information is provided in two ways. First, the projection size combinations are shown in descending size order. Second, the projection size combinations are shown by the number of segregating sites (largest to smallest). This should help users to decide on the best projection sizes for their dataset, which reflects a trade-off between the number of alleles (e.g., projection size) vs. the number of segregating sites.

Using the 2D example file will produce the following:

```
--------------------------------------------------------------------------------
Running test projections for 2D JSFS containing: North, South

Maximum projection sizes = [22, 46]
Minimum projection sizes = [11, 23]


Found 72 projection combinations...


Results sorted by combination order:

	243 segregating sites with projection sizes of (22, 46)
	400 segregating sites with projection sizes of (22, 44)
	495 segregating sites with projection sizes of (22, 42)
	557 segregating sites with projection sizes of (22, 40)
	592 segregating sites with projection sizes of (22, 38)
	628 segregating sites with projection sizes of (22, 36)
	656 segregating sites with projection sizes of (22, 34)
	679 segregating sites with projection sizes of (22, 32)
	697 segregating sites with projection sizes of (22, 30)
	710 segregating sites with projection sizes of (22, 28)
	730 segregating sites with projection sizes of (22, 26)
	742 segregating sites with projection sizes of (22, 24)


	312 segregating sites with projection sizes of (20, 46)
	539 segregating sites with projection sizes of (20, 44)
	696 segregating sites with projection sizes of (20, 42)
	796 segregating sites with projection sizes of (20, 40)
	877 segregating sites with projection sizes of (20, 38)
	935 segregating sites with projection sizes of (20, 36)
	996 segregating sites with projection sizes of (20, 34)
	1,044 segregating sites with projection sizes of (20, 32)
	1,086 segregating sites with projection sizes of (20, 30)
	1,116 segregating sites with projection sizes of (20, 28)
	1,151 segregating sites with projection sizes of (20, 26)
	1,173 segregating sites with projection sizes of (20, 24)


	343 segregating sites with projection sizes of (18, 46)
	609 segregating sites with projection sizes of (18, 44)
	817 segregating sites with projection sizes of (18, 42)
	952 segregating sites with projection sizes of (18, 40)
	1,064 segregating sites with projection sizes of (18, 38)
	1,160 segregating sites with projection sizes of (18, 36)
	1,248 segregating sites with projection sizes of (18, 34)
	1,320 segregating sites with projection sizes of (18, 32)
	1,388 segregating sites with projection sizes of (18, 30)
	1,446 segregating sites with projection sizes of (18, 28)
	1,498 segregating sites with projection sizes of (18, 26)
	1,533 segregating sites with projection sizes of (18, 24)


	378 segregating sites with projection sizes of (16, 46)
	672 segregating sites with projection sizes of (16, 44)
	923 segregating sites with projection sizes of (16, 42)
	1,081 segregating sites with projection sizes of (16, 40)
	1,221 segregating sites with projection sizes of (16, 38)
	1,344 segregating sites with projection sizes of (16, 36)
	1,454 segregating sites with projection sizes of (16, 34)
	1,552 segregating sites with projection sizes of (16, 32)
	1,643 segregating sites with projection sizes of (16, 30)
	1,732 segregating sites with projection sizes of (16, 28)
	1,812 segregating sites with projection sizes of (16, 26)
	1,862 segregating sites with projection sizes of (16, 24)


	398 segregating sites with projection sizes of (14, 46)
	722 segregating sites with projection sizes of (14, 44)
	997 segregating sites with projection sizes of (14, 42)
	1,179 segregating sites with projection sizes of (14, 40)
	1,331 segregating sites with projection sizes of (14, 38)
	1,473 segregating sites with projection sizes of (14, 36)
	1,601 segregating sites with projection sizes of (14, 34)
	1,725 segregating sites with projection sizes of (14, 32)
	1,833 segregating sites with projection sizes of (14, 30)
	1,946 segregating sites with projection sizes of (14, 28)
	2,043 segregating sites with projection sizes of (14, 26)
	2,104 segregating sites with projection sizes of (14, 24)


	425 segregating sites with projection sizes of (12, 46)
	765 segregating sites with projection sizes of (12, 44)
	1,055 segregating sites with projection sizes of (12, 42)
	1,247 segregating sites with projection sizes of (12, 40)
	1,406 segregating sites with projection sizes of (12, 38)
	1,559 segregating sites with projection sizes of (12, 36)
	1,697 segregating sites with projection sizes of (12, 34)
	1,831 segregating sites with projection sizes of (12, 32)
	1,955 segregating sites with projection sizes of (12, 30)
	2,089 segregating sites with projection sizes of (12, 28)
	2,192 segregating sites with projection sizes of (12, 26)
	2,260 segregating sites with projection sizes of (12, 24)


Results sorted by highest number of segregating sites:

	2,260 segregating sites with projection sizes of (12, 24)
	2,192 segregating sites with projection sizes of (12, 26)
	2,104 segregating sites with projection sizes of (14, 24)
	2,089 segregating sites with projection sizes of (12, 28)
	2,043 segregating sites with projection sizes of (14, 26)
	1,955 segregating sites with projection sizes of (12, 30)
	1,946 segregating sites with projection sizes of (14, 28)
	1,862 segregating sites with projection sizes of (16, 24)
	1,833 segregating sites with projection sizes of (14, 30)
	1,831 segregating sites with projection sizes of (12, 32)
	1,812 segregating sites with projection sizes of (16, 26)
	1,732 segregating sites with projection sizes of (16, 28)
	1,725 segregating sites with projection sizes of (14, 32)
	1,697 segregating sites with projection sizes of (12, 34)
	1,643 segregating sites with projection sizes of (16, 30)
	1,601 segregating sites with projection sizes of (14, 34)
	1,559 segregating sites with projection sizes of (12, 36)
	1,552 segregating sites with projection sizes of (16, 32)
	1,533 segregating sites with projection sizes of (18, 24)
	1,498 segregating sites with projection sizes of (18, 26)
	1,473 segregating sites with projection sizes of (14, 36)
	1,454 segregating sites with projection sizes of (16, 34)
	1,446 segregating sites with projection sizes of (18, 28)
	1,406 segregating sites with projection sizes of (12, 38)
	1,388 segregating sites with projection sizes of (18, 30)
	1,344 segregating sites with projection sizes of (16, 36)
	1,331 segregating sites with projection sizes of (14, 38)
	1,320 segregating sites with projection sizes of (18, 32)
	1,248 segregating sites with projection sizes of (18, 34)
	1,247 segregating sites with projection sizes of (12, 40)
	1,221 segregating sites with projection sizes of (16, 38)
	1,179 segregating sites with projection sizes of (14, 40)
	1,173 segregating sites with projection sizes of (20, 24)
	1,160 segregating sites with projection sizes of (18, 36)
	1,151 segregating sites with projection sizes of (20, 26)
	1,116 segregating sites with projection sizes of (20, 28)
	1,086 segregating sites with projection sizes of (20, 30)
	1,081 segregating sites with projection sizes of (16, 40)
	1,064 segregating sites with projection sizes of (18, 38)
	1,055 segregating sites with projection sizes of (12, 42)
	1,044 segregating sites with projection sizes of (20, 32)
	997 segregating sites with projection sizes of (14, 42)
	996 segregating sites with projection sizes of (20, 34)
	952 segregating sites with projection sizes of (18, 40)
	935 segregating sites with projection sizes of (20, 36)
	923 segregating sites with projection sizes of (16, 42)
	877 segregating sites with projection sizes of (20, 38)
	817 segregating sites with projection sizes of (18, 42)
	796 segregating sites with projection sizes of (20, 40)
	765 segregating sites with projection sizes of (12, 44)
	742 segregating sites with projection sizes of (22, 24)
	730 segregating sites with projection sizes of (22, 26)
	722 segregating sites with projection sizes of (14, 44)
	710 segregating sites with projection sizes of (22, 28)
	697 segregating sites with projection sizes of (22, 30)
	696 segregating sites with projection sizes of (20, 42)
	679 segregating sites with projection sizes of (22, 32)
	672 segregating sites with projection sizes of (16, 44)
	656 segregating sites with projection sizes of (22, 34)
	628 segregating sites with projection sizes of (22, 36)
	609 segregating sites with projection sizes of (18, 44)
	592 segregating sites with projection sizes of (22, 38)
	557 segregating sites with projection sizes of (22, 40)
	539 segregating sites with projection sizes of (20, 44)
	495 segregating sites with projection sizes of (22, 42)
	425 segregating sites with projection sizes of (12, 46)
	400 segregating sites with projection sizes of (22, 44)
	398 segregating sites with projection sizes of (14, 46)
	378 segregating sites with projection sizes of (16, 46)
	343 segregating sites with projection sizes of (18, 46)
	312 segregating sites with projection sizes of (20, 46)
	243 segregating sites with projection sizes of (22, 46)

--------------------------------------------------------------------------------

```



## Example Data Set:

In the folder labeled *Example_Data* you will find a 2D and 3D SNPs input file that will run with the `dadi-test-projections.py` script.
You will only need to edit the path to these files in the script, and then the script should run normally. You should test the script using these data to ensure everything is working properly before examining your own empirical data. 


## Citation Information:

The `dadi_pipeline` was originally published as part of the following work:

+ *Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology 26: 5245-5263. https://doi.org/10.1111/mec.14266*

If you use these scripts for your own purposes, please cite this publication.


**Contact:**

Daniel Portik, PhD

daniel.portik@gmail.com

