#!/usr/bin/python3

################################################################################
# upset_summary.py
# Author: Linnea Smeds
# Date: 13-September-2024
# Updated: 11-July-2025, changed GQ to G4 and included an optional chr parameter
# Updated: 03-January-2025, included TRI
# Description: Goes through bedfiles of each non-B motif type (listed below)
# adding each base to a joint matrix with the types as columns and positions 
# as rows.
################################################################################
##### Import libraries
import argparse
import re
from datetime import datetime


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "Finds overlaps of non-B motifs \
from multiple bedfiles")
parser.add_argument('-b', '--bedprefix', help = 'Prefix of the bed files', type = str, required = True)
parser.add_argument('-c', '--chr', help = "Chromosome (if given, only this chr \ will be used)", type = str, required = False)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()

################################################################################
##### Main code

nonB=["APR", "DR", "G4", "IR", "TRI", "STR", "Z"]
#nonB=["APR", "DR", "DRfilt", "G4", "G4quadron", "IR", "IRall","MRall", "TRI", "STR", "Z", "Zgfa", "Zseeker", "ZDNAm1" ]

# prepare dict with positions 
p_dict={}

# Go through the bed files 
for i in range(len(nonB)): 
    file=args.bedprefix+nonB[i]+".bed"
    print("looking at the file: "+file)
    with open(file, 'r') as infile:
        for line in infile:
            tabs = line.rstrip().split("\t")
            if args.chr is not None and tabs[0] == args.chr:
                for j in range(int(tabs[1]),int(tabs[2])):
                    pos=tabs[0]+":"+str(j)
                    if pos in p_dict:
                        p_dict[pos][i]="1"
                    else:
                        p_dict[pos]=["0"]*len(nonB)
                        p_dict[pos][i]="1"
                

# Go through dict and count the different types
c_dict={}
for k in p_dict:
    type="-".join(p_dict[k])
    if type in c_dict:
        c_dict[type]=c_dict[type]+1
    else:
        c_dict[type]=1

# Delete original dict and print summary
p_dict={}
with open(args.output, 'w') as outfile:
    for k in c_dict:
        name=""
        new_list=k.split("-")
        for i in range(len(new_list)):
            if new_list[i] == "1":
                if name == "":
                    name=nonB[i]
                else:
                    name=name+"-"+nonB[i]

        print(name+"\t"+str(c_dict[k]), file=outfile)