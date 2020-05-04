#!/usr/bin/env python

# Author: Yu Fu
# Compare the nonseed region of guides and targets
# Input is a two-column file containing the guide nonseed and target nonseed
# Just select the 2 columns from output of pirna_targete_finder.v1.py
# and use that as the input

import argparse
import sys

def longest_cont_mat(a, b):
    '''It turns out that this is not a good rubic
    Given two strings, it returns the starting position of the
    longest match and the length of the longest match. Gaps disallowed.
'''
    s = -1
    l = 0
    s_tmp = -1
    l_tmp = 0
    prev_matched = False
    for i in range(0, min(len(a), len(b))):
        ai = a[i]
        bi = b[i]
        if ai == bi:
            if prev_matched:
                l_tmp += 1
            else:
                l_tmp = 1
                s_tmp = i
            prev_matched = True
        else:
            # print a[i], b[i]
            if l_tmp > l:
                l = l_tmp
                s = s_tmp
            prev_matched = False
    # A possible case: the longest match goes all the way to the end.
    if l_tmp > l and prev_matched:
        l = l_tmp
        s = s_tmp
    return (s, l)
            
# print longest_cont_mat("AABBCCDDDDDDD", "AXBBCYDDDDDDD")
# print longest_cont_mat("X", "Y")

def n_mat_and_pos(a, b):
    '''Total number of matches between two strings and the matching positions
'''
    mat_pos = []
    ret = 0
    # a = a.replace("U", "T")
    # b = b.replace("U", "T")
    for i in range(0, min(len(a), len(b))):
        if a[i] == b[i]:
            mat_pos.append(i)
            ret += 1
    return (ret, mat_pos)

def n_gu_wobble_and_pos(a, b):
    """The same with n_matches(), except that this function considers 
    GU as a match. Nota that a and b are already in the same direction,
    so GU wobble is actually G-A or T-C
    
    Before:
    5'-GAA-3'
       |||
    3'-UTT-5'
    After rc:
    5'-GAA-3'
    5'-AAA-3'

    Before:
    5'-UAA-3'
    3'-GTT-5'
    After:
    5'-UAA-3'
    5'-CAA-3'
    
"""
    ret = 0
    gu_pos = []
    a = a.replace("U", "T")
    b = b.replace("U", "T")
    for i in range(0, min(len(a), len(b))):
        if a[i] == "G" and b[i] == "A":
            ret += 1
            gu_pos.append(i)
        elif a[i] == "T" and b[i] == "C":
            ret += 1
            gu_pos.append(i)
    return (ret, gu_pos)

# Test
# print n_gu_wobble("GGGAAT", "AGGAAC")
# print n_matches("GGGAAT", "AGGAAC")


parser = argparse.ArgumentParser(description="Compare the nonseed region. \
Input is the nonseed regions of the potential guide and target pairs.")
parser.add_argument("-i", "--input", type=str,
                    help="A two-column file containing the GT pairs")


args = parser.parse_args()
input_fn = args.input

print "#nonseed_guide\tnonseed_target_rc\tn_matches\tn_gu\tmatching_positions\tgu_positions"
with open(input_fn) as fh:
    for line in fh:
        col = line.strip().split()
        g = col[0]
        t_rc = col[1]
        n_gt_gu, gu_pos = n_gu_wobble_and_pos(g, t_rc)
        n_gt_matches, mat_pos = n_mat_and_pos(g, t_rc)
        # pos_longest_cont_mat, n_longest_cont_mat = longest_cont_mat(g, t_rc)
        # pos_longest_cont_mat += 1
        # Use 1-based index for reporting purposes
        gu_pos = [str(i+1) for i in gu_pos]
        mat_pos = [str(i+1) for i in mat_pos]
        gu_pos_str = ",".join(gu_pos) + ","
        mat_pos_str = ",".join(mat_pos) + ","
        print "%s\t%s\t%d\t%d\t%s\t%s" % (g, t_rc, n_gt_matches, n_gt_gu,
                                          mat_pos_str, gu_pos_str)




