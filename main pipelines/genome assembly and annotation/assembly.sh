#!/bin/bash

both=
left=
~/spades/assembler/./spades.py -1 $left -2 $right -o spades_primary -t 32 -careful
~/./platanus assemble -f $left $right -t 56
~/./platanus scaffold -o scaf -c out_contig.fa -b out_contigBubble.fa -IP1 $left $right -t 56 2>gapclose.log
~/./platanus gap_close -o gap -c scaf_scaffold.fa -IP1 $left $right -t 56 2 > gapclose.log
velveth velvet_primary/ 31 -shortPaired -fastq $left $right
velvetg velvet_primary/ -exp_cov auto -cov_cutoff auto -ins_length 350 -min_contig_lgth 1000