#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Seq
# Author: Yu Fu

# 2017-06-08: Only keep the output if the position is 56 57 58 59 60
# 58 is the canonical cut: Miwi cleaves between the 10th and 11th position
# This should decrease the file size as RNAfold is quite slow

# 2017-05-15: Also report the g1 nt in the target
# (that should pair with the 1st nt in piRNA, though it is actually paired)
# to help me determine seed type

# In this version, I will proceed with the minimal assumptions:
# Use 2-7 as the seed. As long as 2-7 are matches, output gt pairs

# SEED_START and SEED_END are 0-based half open half closed ranges
# SEED_START = 1 and SEED_END = 12 mean 2-12 as the seed
SEED_START = 1
SEED_END = 7
SEED_LEN = SEED_END - SEED_START

## g8-21 as the non-seed region
NONSEED_START = SEED_END
NONSEED_END = 21
NONSEED_LEN = NONSEED_END - NONSEED_START

# piRNA species with read counts >= READ_COUNT_CUTOFF are kept
# 0.99 is equivalent to keeping reads with RPM>=1
READ_COUNT_CUTOFF = 0.99

def print_seed_pirna_count(seed2pirna, pirna2count):
    for seed in seed2pirna:
        print(seed)
        for pirna in seed2pirna[seed]:
            print((pirna, pirna2count[pirna]))
        print


def n_matches(a, b):
    '''Total number of matches between two strings
'''
    ret = 0
    # a = a.replace("U", "T")
    # b = b.replace("U", "T")
    for i in range(0, min(len(a), len(b))):
        if a[i] == b[i]:
            ret += 1
    return ret


def n_contig_matches(a, b):
    '''Total number of contiguous matches from the begining
'''
    ret = 0
    for i in range(0, min(len(a), len(b))):
        if a[i] != b[i]:
            break
        else:
            ret += 1
    return ret

# Test the two functions
# one = "abcdefg"
# two = "abczzfg"
# print n_matches(one, two)
# print n_contig_matches(one, two)
        
parser = argparse.ArgumentParser(description="Find piRNA targets. piRNAs \
are used as guide to search for targets in given sequences")
parser.add_argument("-p", "--pirna-insert", type=str, 
                    help="piRNA insert file", required=True)
parser.add_argument("-s", "--sequence", type=str, help="Transcript sequences",
                    required=True)
parser.add_argument("--min-match", type=int, default=6,
                    help="Minimum number of non-seed matches required in g8-21. \
                    This is to prevent the output from being huge")
parser.add_argument("--distance-from-canonical-cut-site", type=int, default=2,
                    help="With the default value of 2, the script will output \
                    cleavage that happens at positions 56, 57, 58, 59, 60. \
                    The canonical cleavage site should be exactly 58, which \
                    is due to the fact that I use 50 nt upstream and \
                    downstream of degradome cleavage sites and that PIWI \
                    proteins cleave between the \
                    10th and 11th position. \
                    This is to prevent the output from being huge")

# parser.add_argument("-w", "--allow-gu-wobble", type=bool, help="Allow GU wobble for seed matching")

args = parser.parse_args()
min_match = args.min_match
distance_from_canonical = args.distance_from_canonical_cut_site

# print args.pirna_insert
# print args.sequence

# Seed to sequence and piRNA to copy number
seed2pirna = {}
seed2pirna_gu = {}
pirna2count = {}

with open(args.pirna_insert) as fh:
    for line in fh:
        col = line.strip().split()
        pirna, count = col
        # This is to make count more general:
        # I may need to consider multimapper someday
        count = float(count)
        seed = pirna[SEED_START: SEED_END]
        if count >= READ_COUNT_CUTOFF:
            if seed in seed2pirna:
                seed2pirna[seed].append(pirna)
            else:
                seed2pirna[seed] = [pirna, ]

            assert pirna not in pirna2count
            pirna2count[pirna] = count

# print_seed_pirna_count(seed2pirna, pirna2count)

seqs = {}
with open(args.sequence) as fh:
    for record in SeqIO.parse(fh, "fasta"):
        seqs[str(record.id)] = record.seq



print("#transcript_id\tzero_based_pos_matching_first_nt_in_seed\tpirna_seed\tpirna\tpirna_count\tpirna_nonseed\ttarget_non_seed_rc\tt1")
for i in seqs:
    seq = seqs[i]
    # Note: use len(seq)-SEED_LEN instead of len(seq)-SEED_LEN+1
    # so that I always have g1 nt available
    for j in range(0, len(seq)-SEED_LEN):
        # Potential target
        pot_t_seed = seq[j:j+SEED_LEN]
        pot_t_seed_rc = pot_t_seed.reverse_complement()
        ## This is used to determine the seed type
        pot_t1 = '*'
        if j+SEED_LEN < len(seq):
            # pot_t1 = seq[j+SEED_LEN]
            pot_t1 = seq[j+SEED_LEN]
        if pot_t_seed_rc in seed2pirna:
            for pirna in seed2pirna[pot_t_seed_rc]:
                pirna_nonseed = pirna[SEED_END:]
                pirna_nonseed_len = len(pirna_nonseed)
                start = j-pirna_nonseed_len
                if start < 0:
                    start = 0
                pot_t_nonseed = seq[start: j]
                pot_t_nonseed_rc = pot_t_nonseed.reverse_complement()
                # Check out how many matched are there in the nonseed
                # gt_n_matches = n_matches(pirna_nonseed, pot_t_nonseed_rc)
                # This is the range of the limited nonseed region
                g_nonseed = pirna_nonseed[:NONSEED_LEN]
                t_nonseed = pot_t_nonseed_rc[:NONSEED_LEN]
                gt_n_matches_nonseed = n_matches(g_nonseed, t_nonseed)

                ## Number of contiguous matches (not needed right now)
                # gt_n_contig_matches = n_contig_matches(pirna_nonseed, pot_t_nonseed_rc)
                # Just output the most relevant information. Number of
                # matches and number of wobble pairs will be processed using other scripts
                if len(pot_t_nonseed_rc) == 0:
                    # Placeholder for those boundary cases
                    pot_t_nonseed_rc = "N"
                if gt_n_matches_nonseed >= min_match and abs(j+SEED_LEN-1-58) <= distance_from_canonical:
                    print("%s\t%d\t%s\t%s\t%f\t%s\t%s\t%s" % (i, j+SEED_LEN-1, pot_t_seed_rc, pirna, pirna2count[pirna], pirna_nonseed, pot_t_nonseed_rc, pot_t1))

