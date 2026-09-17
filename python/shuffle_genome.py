#!/usr/bin/env python3
"""
shuffle_genome.py

Generate a randomized-sequence control chromosome that preserves sequence length,
exact nucleotide composition (and exact positions of N-runs but this doesn't matter
# for the T2T assembly).

Usage:
    python shuffle_genome.py input.fa chrom_name output.fa [seed]

"""


import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def shuffle_preserving_N(seq_str, rng):
    seq_arr = np.frombuffer(seq_str.encode(), dtype='S1').copy()
    is_N = (seq_arr == b'N') | (seq_arr == b'n')

    # Indices of real bases (non-N) get shuffled among themselves;
    # N positions stay exactly where they are.
    real_idx = np.where(~is_N)[0]
    real_bases = seq_arr[real_idx].copy()
    rng.shuffle(real_bases)
    seq_arr[real_idx] = real_bases

    return seq_arr.tobytes().decode()

def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: shuffle_genome.py input.fa chrom_name output.fa [seed]")

    infile, chrom_name, outfile = sys.argv[1:4]
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else None
    rng = np.random.default_rng(seed)

    record = None
    for rec in SeqIO.parse(infile, "fasta"):
        if rec.id == chrom_name:
            record = rec
            break
    if record is None:
        sys.exit(f"Chromosome {chrom_name} not found in {infile}")

    shuffled_seq = shuffle_preserving_N(str(record.seq).upper(), rng)

    new_record = SeqRecord(
        Seq(shuffled_seq),
        id=chrom_name,
        description="shuffled_control"
    )

    with open(outfile, "w") as fh:
        SeqIO.write(new_record, fh, "fasta")

    print(f"Wrote shuffled {chrom_name} ({len(shuffled_seq)} bp) to {outfile}")

if __name__ == "__main__":
    main()
