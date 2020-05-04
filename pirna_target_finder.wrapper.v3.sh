#!/usr/bin/env bash

# RNAfold constraint: g2-21: seed must be matched; the artifical loop must not be matched
# <<<<<<..............xxxxxxx..............>>>>>>

## This script calculates the number of nonseed matches, gu wobbles and MFE

R=${RANDOM}${RANDOM}
## RNAfold does not help much: it cannot rescue targets with 9 or fewer non-seed matches
## Thus, RNAfold process is optional here
# RUN_RNAFOLD=true
RUN_RNAFOLD=false
cpu=56

# Usage: this.sh input.insert output.prefix deg_ref.prefix

# input=WT.sptd.pi6.plus.all.piRNA.20rpm.insert
input=$1
prefix=$2

# prefix=${input%.insert}
tmp_dir=tmp.${prefix}.${R}
# deg_prefix=deg/partitions_multi_only/p
deg_prefix=$3

## Just a hack. Rememebr to make sure to split the fasta file into fewer than 100 parts
for i in {00..99}; do
    mkdir -p ${tmp_dir}
    if [ -f ${deg_prefix}${i}.fa ]; then
	para_prefix=${tmp_dir}/${prefix}
	## With RNAfold
	if [ ${RUN_RNAFOLD} = true ]; then
	    echo "python pirna_target_finder.v3.py -p ${input} -s ${deg_prefix}${i}.fa > ${para_prefix}.gt.${i} && \
grep -v '#' ${para_prefix}.gt.${i} | awk '{ g2to21=substr(\$4, 2, 20); tmp=\$8\$3\$7; if(length(tmp)<21) {tmp=\"NNNNNNNNNNNNNNNNNNNNN\"}; t2to21rc=substr(tmp, 2, 20); print g2to21 > \"/dev/stderr\"; print t2to21rc; }' 2>${para_prefix}.g2to21.${i} | tr \"[ATGCatgc]\" \"[TACGtacg]\" | rev > ${para_prefix}.t2to21.${i} && \
paste ${para_prefix}.t2to21.${i} ${para_prefix}.g2to21.${i} | awk -v constraint='..............((((((xxxxxxx))))))..............' '{ print \$1 \"NNNNNNN\" \$2; print constraint  }' | RNAfold -C | awk 'NR%2==0' > ${para_prefix}.gt2to21.mfe.${i}
"
	else
	## Without RNAfold
	    ## echo "python pirna_target_finder.v3.py -p ${input} -s ${deg_prefix}${i}.fa > ${para_prefix}.gt.${i}"
	    echo "python pirna_target_finder.v3.py -p ${input} -s ${deg_prefix}${i}.fa --distance-from-canonical-cut-site 40 > ${para_prefix}.gt.${i}"
	fi
    fi
done | parallel --progress -j ${cpu}

cat ${tmp_dir}/${prefix}.gt.* | grep -v '#' > ${prefix}.gt

if [ ${RUN_RNAFOLD} = true ]; then
    cat ${tmp_dir}/${prefix}.gt2to21.mfe.* | grep -v '#' > ${prefix}.gt2to21.mfe
fi

## Limit the search in 12-21nt (2-7 are already perfectly matched) (in v1, seed is 2-11)
## Notice that if the nonseed is too short (the predicted target sites are too close to the end of transcripts, they will have X\tY as placeholder
grep -v '#' ${prefix}.gt | awk -v OFS="\t" -v nonseed_len=14 '{ a=substr($6, 1, nonseed_len); b=substr($7, 1, nonseed_len); if(length(a)==length(b) && length(a)==nonseed_len) { print a "\t" b } else { print "X\tY" } }' > ${prefix}.gt.nonseed_limited
python gt_nonseed_comparator.py -i ${prefix}.gt.nonseed_limited > ${prefix}.gt.nonseed_limited.match

## Columns: #transcript_id zero_based_pos_matching_first_nt_in_seed pirna_seed pirna pirna_count pirna_nonseed target_non_seed_rc t1
## pirna_g8-21 target_t8-21 n_matches n_gu match_pos gu_pos sec_struct mfe

## The sed command here is to make sure the last column of RNAduplex output is good for R parising
if [ ${RUN_RNAFOLD} = true ]; then
    paste <(grep -v '#' ${prefix}.gt) <(grep -v '#' ${prefix}.gt.nonseed_limited.match) <(cat ${prefix}.gt2to21.mfe | sed 's/ /\t/') > ${prefix}.gt2
else
    paste <(grep -v '#' ${prefix}.gt) <(grep -v '#' ${prefix}.gt.nonseed_limited.match) > ${prefix}.gt2
fi

## Comment these two commands for degbugging purposes
rm ${prefix}.gt.nonseed_limited ${prefix}.gt.nonseed_limited.match
rm ${prefix}.gt

