#!/usr/bin/python3

################################################################################
# repeat_summary.py
# Author: Linnéa Smeds
# Date: 21-August-2024
# Description: Takes a bed file from stdin with names in 4th column, and sums up
# the total length for each type. Prints a list (names and lengths) to stdout.
################################################################################
##### Import libraries
import fileinput

################################################################################
##### Main code

# SAVE REPEATS AND LENGTHS IN DICT
r_dict={}
with fileinput.input() as f_input:
    for line in f_input:
        tabs = line.rstrip().split("\t")
        #print(tabs)
        if tabs[3] not in r_dict:
            r_dict[tabs[3]]=int(tabs[2])-int(tabs[1])
        else:
            r_dict[tabs[3]]=r_dict[tabs[3]]+int(tabs[2])-int(tabs[1])

# PRINT THE DICT (SORTED)
for rep in sorted(r_dict.keys()):
    print(rep+"\t"+str(r_dict[rep]))
