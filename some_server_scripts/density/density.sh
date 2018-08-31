#!/bin/bash
echo "input your seq!"
read seq
touch dens.txt
n_contig=$(cat ~/genomes/dmel_test/refgenome/dmel-all-chromosome-r6.19.fasta | grep ">" | wc -l)
echo "Contig count is $n_contig"
#cat ~/genomes/dmel_test/refgenome/dmel-all-chromosome-r6.19.fasta | grep -aob $seq | grep -a [0-9] | sed 's/[^0-9]*//g' >> dens.txt
#Rscript density.R
#rm dens.txt
