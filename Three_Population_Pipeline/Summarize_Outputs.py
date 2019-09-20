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

#initiate empty lists that we will fill with summary information
summary_list = []
simple_list = []

#list comprehension to find output files
flist = sorted([f for f in os.listdir('.') if f.endswith(".optimized.txt")])
print("\n\nFound {} output files to summarize.\n".format(len(flist)))

#iterate over output files
for f in flist:
    content = []
    print("\tExtracting contents from: {}".format(f))
    #open file, skip first line
    with open(f, 'r') as fh:
        next(fh)
        #split lines and add contents to list
        for line in fh:
            content.append(line.strip().split('\t'))
    #content list items will have order: "Model"	"Replicate"	"log-likelihood"	"AIC"	"chi-squared"	"theta"	"optimized_params(xxx)"
    #let's sort all the rows by AIC, lowest to highest
    content.sort(key=lambda x: float(x[3]))
        
    #add top five entries to summary list
    for i in content[:5]:
        summary_list.append(i)
        
    #add top entry to easy list
    simple_list.append(content[0])
   
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
