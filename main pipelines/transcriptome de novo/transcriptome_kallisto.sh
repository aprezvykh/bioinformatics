kallisto index -i index/transcripts.idx index/Drosophila_melanogaster.BDGP6.dna.toplevel.fa
kallisto quant -t 32 --single -l 200 -s 20 -i index/transcripts.idx -o Kallisto_output/K5/Quants -b 64 K_5.fastq
kallisto quant -t 32 --single -l 200 -s 20 -i index/transcripts.idx -o Kallisto_output/K6/Quants -b 64 K_6.fastq
kallisto quant -t 32 --single -l 200 -s 20 -i index/transcripts.idx -o Kallisto_output/F3/Quants -b 64 F_3.fastq
kallisto quant -t 32 --single -l 200 -s 20 -i index/transcripts.idx -o Kallisto_output/F4/Quants -b 64 F_4.fastq