'''
Usage: dadi-test-projections.py

The purpose of this script is to quickly generate combinations of projection sizes 
to determine the resulting numbers of segregating sites. Information is entered 
for the JSFS, including the maximum projection sizes (in alleles) and the minimum fraction 
of alleles to retain for projections. The resulting segregating sites for each projection 
combination is shown on screen. The results are first sorted by descending sizes, and then 
by the maximum number of segregating sites.

The sections with #************** must be edited using information for your dataset. 

This information should help users decide on the best projection sizes, reflecting 
a trade-off between the number of alleles vs. segregating sites.

The script should work for any 2D or 3D snps input format. Example files are provided 
for both 2D and 3D files, and will run with the settings below (but edit the path to 
the input files based on their downloaded locations!).

-------------------------
Written for Python 3.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated October 2020
'''

import dadi
from itertools import product
from datetime import datetime

#===========================================================================
# Function to generate and test all projection combinations, prints to screen
#===========================================================================

def run_projections(dd, pop_ids, polarized, maxproj, min_frac):
    # get minprojections for all maxproj items using fraction
    minproj = [int(x * min_frac) for x in maxproj]
    # initiate empty list to store projection ranges
    sizes = []
    # iterate over each item in the maxprojections list
    for i in range(0, len(maxproj)):
        # create a list of values, max to min by decreasing increments of 2
        sizes.append(list(range(maxproj[i], minproj[i], -2)))
    # create all possible combinations of projections using itertools product module
    sizecombinations = list(product(*sizes))
    # print info for this jsfs and projections
    print("\n\n\n{}\nRunning test projections for {}D JSFS containing: {}\n".format("-"*80, len(pop_ids), ", ".join(pop_ids)))
    print("Maximum projection sizes = {}\nMinimum projection sizes = {}".format(maxproj, minproj))
    print("\n\nFound {} projection combinations...".format(len(sizecombinations)))

    # list to store projection sizes and segregating sites
    results = []
    # iterate over all projection combinations
    for combo in sizecombinations:
        # create a spectrum using projection sizes
        fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = combo, polarized = polarized)
        #print("projection: {},\tsites: {}".format(combo, int(fs.S())))
        # add sizes and segregating sites in a sub-list to results list
        results.append([combo, int(fs.S())])
        
    # show results by descending order of projections
    print("\n\nResults sorted by combination order:\n")
    for i in range(0, len(results)):
        if i != 0:
            if results[i][0][0] != results[i-1][0][0]:
                print("\n")
        print("\t{1:,} segregating sites with projection sizes of {0}".format(results[i][0], results[i][1]))

    # show results by highest to lowest number segregating sites
    results.sort(key=lambda x: x[1], reverse=True)
    print("\n\nResults sorted by highest number of segregating sites:\n")
    for r in results:
        print("\t{1:,} segregating sites with projection sizes of {0}".format(r[0], r[1]))
    print("\n{}\n\n".format("-"*80))


#===========================================================================
# Example for a 2D joint-site frequency spectrum
#===========================================================================

#**************
snps = "/Users/portik/Documents/GitHub/dadi_pipeline/Find-Best-Projections/Example_data/dadi_2pops_North_South_snps.txt"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["North", "South"]

#**************
#Indicate whether your frequency spectrum object is polarized/unfolded (set to True) or unpolarized/folded (set to False)
polarized = False

#**************
# Maximum projection sizes. This is the maximum possible ALLELES (not individuals) that can be
# found in each of the populations. 
maxproj = [22, 46]

# Choose the minimum fraction of alleles required for both populations.
# For example, if you want to include a minimum of 50% of the alleles, use min_frac = 0.5.
# For a minimum of 75% of alleles, use min_frac = 0.75, etc.
min_frac = 0.5

run_projections(dd, pop_ids, polarized, maxproj, min_frac)





#===========================================================================
# Example for a 3D joint-site frequency spectrum
#===========================================================================

#**************
snps = "/Users/portik/Documents/GitHub/dadi_pipeline/Find-Best-Projections/Example_data/dadi_3pops_CVLS_CVLN_Cross_snps.txt"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["CVLS", "CVLN", "Cross"]

#**************
#Indicate whether your frequency spectrum object is polarized/unfolded (set to True) or unpolarized/folded (set to False)
polarized = False

#**************
# Maximum projection sizes. This is the maximum possible ALLELES (not individuals) that can be
# found in each of the populations. 
maxproj = [26, 46, 20]

# Choose the minimum fraction of alleles required for both populations.
# For example, if you want to include a minimum of 50% of the alleles, use min_frac = 0.5.
# For a minimum of 75% of alleles, use min_frac = 0.75, etc.
min_frac = 0.5

run_projections(dd, pop_ids, polarized, maxproj, min_frac)
