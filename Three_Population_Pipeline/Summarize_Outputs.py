import sys
import os
import shutil
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
Written for Python 2.7
-------------------------

Daniel Portik
danielportik@email.arizona.edu
https://github.com/dportik
May 2018
'''

#===========================================================================
file_dir = sys.argv[1]
os.chdir(file_dir)

#initiate empty lists that we will fill with summary information
summary_list = []
simple_list = []


#search working directory
for filetype in os.listdir('.'):
    #find output files by extension name
    if filetype.endswith(".optimized.txt"):
        #now to get contents into an empty list
        content = []
        #open file
        fh_temp = open(filetype, 'r')
        print "Examining file {}".format(filetype)
        #read lines into memory
        lines = fh_temp.readlines()
        #iterate over lines but skip header line
        for line in lines[1:]:
            line = line.strip()
            line_items = line.split('\t')
            content.append(line_items)
        #close file
        fh_temp.close()
        
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
fh_out1 = open(out1, 'a')

#write tab-delimited file with extended results
fh_out1.write("Model\tReplicate\tlog-likelihood\tAIC\tchi-squared\ttheta\toptimized_params\n")
for row in summary_list:
    for val in row:
        fh_out1.write("{}\t".format(val))
    fh_out1.write("\n")

#create output file 2
out2 = "Results_Summary_Short.txt"
fh_out2 = open(out2, 'a')

#write tab-delimited file with simplified results
fh_out2.write("Model\tReplicate\tlog-likelihood\tAIC\tchi-squared\ttheta\toptimized_params\n")
for row in simple_list:
    for val in row:
        fh_out2.write("{}\t".format(val))
    fh_out2.write("\n")

#close output files
fh_out1.close()
fh_out2.close()

print "\n\nSummary files '{0}' and '{1}' have been written to {2}\n\n".format(out1, out2, file_dir)

#===========================================================================
