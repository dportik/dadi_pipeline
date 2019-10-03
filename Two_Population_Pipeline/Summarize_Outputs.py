'''
usage: python Summarize_Outputs.py [full path to directory with results files]

example: python Summarize_Outputs.py users/dan/moments_analyses/results/

The purpose of this script is to make sense of the results files
that are produced after optimizations have been run on multiple
models, as there will be an output file for every model included
that is filled with many replicates across several rounds.

Two summary files will be produced:
1. Results_Summary_Extended.txt
    This contains the top five replicates for each results file, with all the information
    including: "Model"	"Replicate"	"log-likelihood"	"AIC"	"chi-squared"	"theta"	"optimized_params"

2. Results_Summary_Short.txt
    This is essentially a simplified version of the above file, and only contains the
    top scoring replicate per model and it is already sorted in order of AIC.

You should probably inspect that the top-scoring replicates were ERROR-FREE. 
Often errors are logged on screen and log-likelihoods are still produced.

-------------------------
Written for Python 2.7 and 3.7
No dependencies
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated September 2019
'''

import sys
import os

#===========================================================================
file_dir = sys.argv[1]
os.chdir(file_dir)

#check if summary files already exist in output directory
results = sorted([os.path.abspath(f) for f in os.listdir(file_dir) if f.startswith("Results_Summary")])
if results:
    raise ValueError("\n\n\nWARNING: Summary files are already located in the directory specified."
                         " Please remove them before running this script.\n\n")
    

#initiate empty lists that we will fill with summary information
summary_list = []
simple_list = []

#list comprehension to find output files
flist = sorted([os.path.abspath(f) for f in os.listdir(file_dir) if f.endswith("optimized.txt")])

#do tasks depending on whether files located or not
if flist:
    print("\n\nFound {} output files to summarize:".format(len(flist)))
    for f in flist:
        print("\t{}".format(f.split('/')[-1]))
else:
    print("\n\nFound {} output files to summarize!\n\n".format(len(flist)))
    print("Please check to make sure the output files end with '.optimized.txt' and are"
              "located in the directory specified:\n\t{}".format(file_dir))

#iterate over output files
for f in flist:
    print("\nExtracting contents from: {}".format(f.split('/')[-1]))
    #content list items will have order: "Model"	"Replicate"	"log-likelihood"	"AIC"	"chi-squared"	"theta"	"optimized_params(xxx)"
    with open(f, 'r') as fh:
        lines = [line.strip().split('\t') for line in fh if not line.startswith("Model")]
    #strict filtering: sublists must have 7 elements and not contain any "nan" entries
    content = [l for l in lines if len(l) == 7 and "nan" not in l]
    if len(lines) == len(content):
        print("\tFound {} total replicates.".format(len(lines)))
    else:
        print("\tFound {} row entries.".format(len(lines)))
        print("\tRemoved {} row entries due to presence of 'nan' values or incomplete data.".format(len(lines)-len(content)))
        
    #let's sort all the rows by AIC, lowest to highest
    content.sort(key=lambda x: float(x[3]))
        
    #add top five entries to summary list
    for i in content[:5]:
        summary_list.append(i)
        
    #add top entry to easy list
    simple_list.append(content[0])
   
#make sure results are actually in lists
if simple_list and summary_list:
    
    #sort the list containing only the top entry for each model by order of AIC
    simple_list.sort(key=lambda x: float(x[3]))

    #create output file 1
    out1 = "Results_Summary_Extended.txt"
    with open(out1, 'a') as fh:
        #write tab-delimited file with extended results
        fh.write("Model\tReplicate\tlog-likelihood\tAIC\tchi-squared\ttheta\toptimized_params\n")
        for row in summary_list:
            for val in row:
                fh.write("{}\t".format(val))
            fh.write("\n")

    #create output file 2
    out2 = "Results_Summary_Short.txt"
    with open(out2, 'a') as fh:
        #write tab-delimited file with simplified results
        fh.write("Model\tReplicate\tlog-likelihood\tAIC\tchi-squared\ttheta\toptimized_params\n")
        for row in simple_list:
            for val in row:
                fh.write("{}\t".format(val))
            fh.write("\n")

    print("\n\nSummary files '{0}' and '{1}' have been written to: \n\t{2}\n\n".format(out1, out2, file_dir))

else:
    print("\n\nNo results were written!\n\n")
