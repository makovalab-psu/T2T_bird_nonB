#!/usr/bin/python3

################################################################################
# repam_zdna.py
# Written by Linnéa Smeds using the help of generative AI
# Date: 21-August-2024
# Description: Takes full genome output from ZDNA Hunter and converts the 
# positions to chromsome coordinates. 
################################################################################
##### Import libraries
import argparse
import bisect

################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "Converts Z-DNA Hunter genome output \
to bed format with proper chromosome coordinates")
parser.add_argument('-b', '--bedgraph', help = 'Z-DNA Hunter bedgraph file for the full genome', type = str, required = True)
parser.add_argument('-f', '--fai', help = "Genome fasta index file (Need to be ordered \
the same as the Z-DNA Hunter input genome!)", type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()


################################################################################
##### Define functions 


def read_fai(fai_path):
    """
    Reads a FASTA index (.fai) and returns:
    - chroms: list of chromosome names in FASTA order
    - starts: list of cumulative start positions
    """
    chroms = []
    starts = []
    cum = 0

    with open(fai_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.split()
            chrom = fields[0]
            length = int(fields[1])

            chroms.append(chrom)
            starts.append(cum)
            cum += length

    return chroms, starts


def remap_zdna(fai_path, zdna_path, out_path):
    chroms, starts = read_fai(fai_path)

    with open(zdna_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if not line.strip():
                continue

            fields = line.split()

            # Skip header or track lines
            if len(fields) < 4:
                continue
            if not fields[1].isdigit():
                continue

            gstart = int(fields[1])
            gend   = int(fields[2])
            score  = fields[3]

            # Find chromosome via binary search
            idx = bisect.bisect_right(starts, gstart) - 1

            chrom = chroms[idx]
            chr_start = gstart - starts[idx]
            chr_end   = gend   - starts[idx]

            fout.write(f"{chrom}\t{chr_start}\t{chr_end}\t{score}\n")

################################################################################
##### Main code


if __name__ == "__main__":

    remap_zdna(args.fai, args.bedgraph, args.output)

