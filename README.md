# GTBuster
A set of scripts for looking piRNA guide-target pairs by integrating small RNA-seq and degradome-seq data

## Usage

`bash pirna_target_finder.wrapper.v3.sh insert output_prefix ref/p`

`insert` is a tab-delimited file containing two columns: piRNA sequence and its normalized abundance (e.g., reads per million).

```
TGATATGTAGGATCCCTCTCAGATGTCAAG  28.1522
TGTTAAGTTCTAAAGGTTCTCTCTGAAGATC 25.9157
TGATATGTAGGATCCCTCTCAGATGTCAAGT 18.283
TGTTAAGTTCTAAAGGTTCTCTCTGAAGAT  10.8278
TGATATGTAGGATCCCTCTCAGATGTCAAGA 10.5083
TAAAAAAATGTAAGGAATTGTTAAGGACAGC 9.72725
```

`ref/p` is the prefix of many fasta files, such as `ref/p0`, `ref/p1`, etc. Collectively, these files represent the cleavage sites of degradome-seq reads. Below are a few example lines:

```
>chr1|4142665|4142666|-|1.96241|2.08546(-)
GCTGTTGCACTGTGAGGCCAGTGACAAATCCCAGTACTGGTACTGTGAGA
AGATCATCGTCAAAGATCCAGGCAGCTCCTCAGAATCCATCTTCACCTGT
G
>chr1|4163909|4163910|-|1.28572|0.678986(-)
GACTTTGCATTTCATCCAGGTTGTACTAAGAGGCATTGGAGAGATATACA
AGATCAGAATAGGACATGATGGAACTGGTGGGCAACCAGAGTGGACCTTA
C
>chr1|4170396|4170397|-|1.96241|0.339493(-)
TGTACCTGATCCGCCCCTTATGTTCTTCTCTTACTCTGTACAGAAGAAAG
CAAGCATCTGTGCAGGATTGAGACACTTCTGTTAGAAGAAAAATGGAAGG
T
```

## Output format

The output is a tab-delimited files describing the piRNA guides and their corresponding targets:

```
chr1|78665204|78665205|+|1.47849|0.670895(+)    19      CTAGCC  TCTAGCCTCAGACTGAGACTTGTACAAATGT 1.232850        TCAGACTGAGACTTGTACAAATGT
        TCATAATAGGACTG  A       TCAGACTGAGACTT  TCATAATAGGACTG  9       1       1,2,3,5,7,10,11,12,13,  8,
chr1|78665204|78665205|+|1.47849|0.670895(+)    20      TCTAGC  TTCTAGCCTGTTCTTACTTCAAAGATGCCT  3.857610        CTGTTCTTACTTCAAAGATGCCT CTCATAATAGGACTG G       CTGTTCTTACTTCA  CTCATAATAGGACT  6       0       1,2,5,8,9,13,   ,
chr1|78665204|78665205|+|1.47849|0.670895(+)    20      TCTAGC  TTCTAGCCTGTTCTTACTTCAAAGATGCCTGC        1.511230        CTGTTCTTACTTCAAAGATGCCTGC       CTCATAATAGGACTG G       CTGTTCTTACTTCA  CTCATAATAGGACT  6       0       1,2,5,8,9,13,   ,
chr1|78665204|78665205|+|1.47849|0.670895(+)    20      TCTAGC  TTCTAGCCTGTTCTTACTTCAAAGATGCCTG 1.352150        CTGTTCTTACTTCAAAGATGCCTG
        CTCATAATAGGACTG G       CTGTTCTTACTTCA  CTCATAATAGGACT  6       0       1,2,5,8,9,13,   ,
chr1|78665204|78665205|+|1.47849|0.670895(+)    20      TCTAGC  TTCTAGCCTGTTCTTACTTCAAAGATGCCTGT        1.232850        CTGTTCTTACTTCAAAGATGCCTGT       CTCATAATAGGACTG G       CTGTTCTTACTTCA  CTCATAATAGGACT  6       0       1,2,5,8,9,13,   ,
chr1|78665204|78665205|+|1.47849|0.670895(+)    21      CTCTAG  TCTCTAGCTTTGTTCTGGAAGGTGATCTGT  1.988460        CTTTGTTCTGGAAGGTGATCTGT CCTCATAATAGGACTG        C       CTTTGTTCTGGAAG  CCTCATAATAGGAC  6       4       1,3,6,9,11,13,  2,4,5,10,
```

Column 1: description of the cleavage site from degradome-seq: chr|start position|end position|strand|normalized degradome read counts under condition 1 (WT)|normalized degradome read counts under condition 2 (mut). The extra strand in the end is not useful, just a byproduct of some function for sanity check.

Column 2: the position of the cleavage (0-based). The reference is 50nt upstream + degradome cleavage site (1nt) + 50nt downstream. 58 means the cut happens with the exactly 10nt overlap of piRNA and target sites. I believe only rows with 58 in this column are kept in these files.

Column 3: 2–7 nt of the piRNA, i.e., the seed portion (or, at least, this is what I use for initial identification of potential target pairs)

Column 4: piRNA sequence

Column 5: normalized piRNA abundance (RPM)

Column 6: the remaining portion of the piRNA (column 3 subtracted from column 4), i.e., so-called non-seed portion of piRNA (assuming 2–7 nt as the seed region)

Column 7: the part of the target that corresponds to the non-seed part of piRNA. This is reverse-complemented to make it easier for a human to read.

Column 8: the part of the target corresponding to the first nt of piRNA. Pairing of this position is not required, just like in Wei's paper

Columns 9&10: please ignore these, as these are just short versions of columns 7&8

Column 11–14: # matches, # GU wobble pairs, matching positions, GU wobble positions


