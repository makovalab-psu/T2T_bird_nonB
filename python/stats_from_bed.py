#!/usr/bin/python3

################################################################################
# stats_from_bed.py
# Author: Linnea Smeds
# Originally written 23-January-2025 for the T2T primates project.
# Description: Takes a bed file and a desired output (median, mean, max, min)
################################################################################
##### Import libraries
import argparse
import re
import numpy as np
from datetime import datetime


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script matches bed entries \
with a list of desired elements")
parser.add_argument('-b', '--bedfile', help = 'Bed input file name', type = str, required = True)
parser.add_argument('-p', '--op', help = 'Operation to perform', type = str, required = True)
args = parser.parse_args()

################################################################################
##### Main code


# SAVE LENGTHS IN LIST
llist=[]
n=0
with open(args.bedfile, 'r') as infile:
    for line in infile:
        tabs = line.rstrip().split("\t")
        llist.append(int(tabs[2])-int(tabs[1]))
        n=n+1

llist=np.array(llist)

if args.op == "median":
    print(np.median(llist))
elif args.op == "mean":
    print(np.mean(llist))
elif args.op == "min":
    print(np.min(llist))
elif args.op == "max":
    print(np.max(llist))
elif args.op == "all":
    print(str(n)+" "+str(np.median(llist))+" "+str(np.mean(llist))+" "+str(np.min(llist))+" "+str(np.max(llist)))
else:
    print("ERROR: operation "+args.op+" is not supported.")
