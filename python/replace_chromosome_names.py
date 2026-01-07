#!/usr/bin/python3

################################################################################
#replace_chromosome_names.py
# Author: Linnea Smeds
# Date: 21-April-2025
# Description: Takes a translation table with new and old names, and replaces
# the old names in a gtf/gff file with the new names. The translation table is
# expected to be a two-column file with the old names in the first column and
# the new names in the second column. 
# Usage: python replace_chromosome_names.py -g <gtffile> -t <table> -o <output>
# where -t is the translation table.
################################################################################
##### Import libraries
import argparse
import re
from datetime import datetime


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script replaces chromosome \
                    names in a file based on a translation table.")
parser.add_argument('-g', '--gtffile', help = 'gtf/gff file name', type = str, required = True)
parser.add_argument('-t', '--transl', help = 'Translation table with old name \
                     in first column, new in second', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()

################################################################################
##### MAIN CODE 

# SAVE NAMES IN DICT
t_dict={}
with open(args.transl, 'r') as infile:
    for line in infile:
        tabs = line.rstrip().split("\t")
        if tabs[0] not in t_dict:
            t_dict[tabs[0]]=tabs[1]
        else:
            print("ERROR: Problem with translation table: duplicate name", tabs[0])


# GO THROUGH THE GTF/GFF
with open(args.output, 'w') as outfile:
    with open(args.gtffile, 'r') as infile:
        for line in infile:
            s=re.match("^#", line)
            if s:
                h=re.match("^##sequence-region", line)
                if h:
                    tabs = line.rstrip().split(" ")
                    if tabs[1] in t_dict:
                        tabs[1] = t_dict[tabs[1]]
                        print(" ".join(tabs), file=outfile)
                    else:
                        print("Error: Cannot find translation of", tabs[1], "in translation table")
                        print(line.rstrip(), file=outfile)
                else:
                    print(line.rstrip(), file=outfile)
            else:
                tabs = line.rstrip().split("\t")
                if tabs[0] in t_dict:
                    tabs[0] = t_dict[tabs[0]]
                    print("\t".join(tabs), file=outfile)
                else:
                    print("Error: Cannot find translation of", tabs[0], "in translation table")
                    print(line.rstrip(), file=outfile)

