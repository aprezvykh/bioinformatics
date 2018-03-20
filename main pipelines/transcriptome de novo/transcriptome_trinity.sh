#!/bin/bash

~/sofware/trinityrnaseq/./Trinity --seqType fq --single SRR1533748_1.fastq --CPU 24 --max_memory 100G
~/sofware/trinityrnaseq/util/./TrinityStats.pl trinity_out_dir_primary_assembly/Trinity.fasta
bowtie2-build
bowtie2 -x b2_index_spades/index_spades SRR1533748_1.fastq -p 16 --no-unal -S spades_align.sam > spades.self-align.simmary.txt
blastx -query trinity_out_dir_primary_assembly/Trinity.fasta -db ~/uniprot_db/uniprot_sprot.fasta -out blastx.outfmt6 -num_threads 32 -evalue 1e-20 -max_target_seqs 1 -outfmt 6
./analyze_blastPlus_topHit_coverage.pl spades.transdecoder.cds.blastx.outfmt6 transdecoder_out/transcripts.fasta.transdecoder.cds ~/uniprot_db/uniprot_sprot.fasta | column -t


~/TransDecoder/./TransDecoder.LongOrfs -t spades_assembly/transcripts.fasta
Use file: 101M_1_GG.fasta.transdecoder_dir/longest_orfs.pep
~/TransDecoder/./TransDecoder.Predict -t spades_assembly/transcripts.fasta
~/busco/scripts/./run_BUSCO.py -i trinity_assembly_w_trimmomatic/Trinity.fasta -o trinity_busco --lineage_path ~/busco_datasets/diptera_odb9/ -m transcriptome -c 24
python ~/busco/scripts/generate_plot.py -wd /mnt/raid/illumina/AlexR/transcriptomes/reads/BUSCO_summaries/
ls *.fastq | parallel -j 4 'STAR --genomeDir ~/transcriptomes/ref/virilis/dvir_STAR_index/ --readFilesIn {} --runThreadN 8 --outFileNamePrefix genomeGuided_{} --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate'
