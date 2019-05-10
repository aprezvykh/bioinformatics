#!/bin/bash
RED='\033[0;31m'
NC='\033[0m'
echo "This script is used to found degenerate gRNA in CRISPR-Cas9 system"
echo "Developed by Alexander Rezvykh, aprezvykh@yandex.ru"
echo "______________________________________________________________________________________________________________________________"
echo "Arguments:"
echo "-g|--genome_fasta (filepath) - genome file in FASTA file"
echo "-s|--spacer_length (integer) - length of a spacer sequence exclude PAM site (20, for example)"
echo "-u|--unallowed_spacer_string (character) - sequence that not allowed in spacer sequence (TTT, for example)"
echo "-p|--unallowed_pam_end (character) - PAM sequence, that unallowed (AA, for example)"
echo "-t|--number_of_threads (integer) - threads number"
echo "-w|--word_size (integer) - size of the word in blast "
echo "-a|--annotation_file (filepath) - annotation file in GTF format"
echo "-d|--debug (logical, T/F) - uses only first 10 PAM sequences "
echo "-pr | -prefix (will be used in final result file)"

#banner ROPSIR

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -g|--genome_fasta)
    genome="$2"
    shift
    shift
    ;;
    -s|--spacer_length)
    spacer_length="$2"
    shift
    shift 
    ;;
    -u|--unallowed_spacer_string)
    unallowed_spacer_string="$2"
    shift
    shift 
    ;;
    -p|--unallowed_pam_end)
    unallowed_pam_end="$2"
    shift
    shift
    ;;
    -t|--numbwer_of_threads)
    threads="$2"
    shift 
    shift
    ;;
   -w|--word_size)
    word_size="$2"
    shift
    shift
    ;;
    -d|--debug)
    debug="$2"
    shift
    shift
    ;;
    -a|--annotation_file)
    annotation_file="$2"
    shift
    shift
    ;;
    -pr|--prefix)
    prefix="$2"
    shift
    shift
    ;;
esac
done

if [[ -z "$genome" ]]
	then
		echo "Genome must be specified, nowhere to look for spacer sites!"
		exit 0
	fi

if [[ -z "$prefix" ]]
        then
                echo "Prefix non set! Setting..."
                prefix=$(date +"%m-%d-%Y-%H-%M")
        fi


if [[ -z "$annotation_file" ]]
	then
		echo "Genome annotation file must be specified!"
		exit 0
	fi

if [[ -z "$spacer_length" ]]
	then
		echo "Spacer length not set! Setting default (20)"
		spacer_length=20
	fi

if [[ -z "$unallowed_spacer_string" ]]
	then
		echo "Unallowed spacer string not set! Setting default (TTT)"
		unallowed_spacer_string="TTT"
	fi

if [[ -z "unallowed_pam_end" ]]
	then
		echo "Unallowed PAM end has not set! Setting default (AA)"
		unallowed_pam_end="AA"
	fi

if [[ -z "$threads"  ]]
	then
		echo "Threads variable not set! Setting default (nproc output)"
		threads=$(nproc)
	fi

set -- "${POSITIONAL[@]}" # restore positional parameters

echo GENOME FILE IS SET TO = "${genome}"
echo SPACER LENGTH IS SET TO   = "${spacer_length}"
echo UNALLOWED STRINGS IN SPACERS IS SET TO    = "${unallowed_spacer_string}"
echo UNALLOWER PAM END IS SET TO         = "${unallowed_pam_end}"
echo THREADS IS SET TI = "${threads}"
echo SCRIPT EXECUTION DIR IS $pwd

all_ngg_sequences=potential_ngg.fasta
all_ngg_sequences_dbg=potential_dbg_ngg.fasta
all_ngg_sequences_space=potential_ngg.fasta.parsed
final_spacers=ngg.headers.fasta

cpu_n=$(nproc)
if [[ $cpu_n -lt 5 ]]
	then
		echo "${RED}WARNING! CPU CORES ON YOUR SYSTEM IS LESS THEN 4! EXECUTING OF THIS SCRIPT MAY TAKE A WHILE!${NC}"
	else
		echo "Thread number is $cpu_n"
	fi


is_blastn=$(which blastn | wc -l)
if [[ $is_blastn -eq 1 ]]
	then
		echo "blastn found!"
	elif [[ $is_blastn -ge 1 ]]
	then
		echo "Multiple blastn versions found "
	elif [[ $is_blastn -eq 0 ]]
	then
		echo -e "Cannot found blastn in system! Install blastn: ${RED}sudo apt-get install ncbi-blast+${NC}"
		exit 0
	fi

is_samtools=$(which samtools | wc -l)
if [[ $is_samtools -eq 1 ]]
        then
                echo "samtools found!"
        elif [[ $is_samtools -ge 1 ]]
        then
                echo "Multiple samtools versions found "
        elif [[ $is_samtools -eq 0 ]]
        then
                echo -e "Cannot found samtools in system! Install samtools: ${RED}sudo apt-get install samtools${NC}"
                exit 0
        fi

is_rlang=$(which Rscript | wc -l)
if [[ $is_rlang -eq 1 ]]
        then
                echo "rscript found!"
        elif [[ $is_rlang -ge 1 ]]
        then
                echo "Multiple rscript versions found "
        elif [[ $is_rlang -eq 0 ]]
        then
                echo -e "Cannot found blastn in system! Install rscript: ${RED}sudo apt-get install r-base-core${NC}"
                exit 0
        fi

is_blast2bam=$(which blastxmlparser | wc -l)
if [[ $is_blast2bam -eq 1 ]]
        then
                echo "blastxmlparser found!"
        elif [[ $is_blast2bam -eq 0 ]]
        then
                echo -e "Cannot found blastxmlparser in system! Install blastxmlparser: ${RED}sudo gem instal blastxmlparser${NC}"
                exit 0
        fi

is_rnafold=$(which RNAfold | wc -l)
if  [[ $is_rnafold -eq 1 ]]
	then
		echo "RNAfold found"
	else
		echo "RNAfold not found! Install RNAfold: ${RED}sudo apt-get install rnafold${NC}"

	fi

is_ssconvert=$(which ssconvert | wc -l)
if  [[ $is_rnafold -eq 1 ]]
        then
                echo "ssconvert found"
        else
                echo "ssconvert not found! Install RNAfold: ${RED}sudo apt-get install ssconvert${NC}"

        fi


genome_size=$(cat $genome | wc | awk '{print $3-$1}')
spacer_regexp=$(yes '.' | head -n $spacer_length | tr -d '\n')

echo "Regular expression used in search is $spacer_regexp"

cat $genome | grep -oh "$spacer_regexp.[AG][AG]" | grep -v "$unallowed_spacer_string" | grep -v "$unallowed_pam_end$" > $all_ngg_sequences

ngg_length=$(cat $all_ngg_sequences | wc -l)

echo "Genome size is $genome_size"

if [[ $ngg_length -eq 0 ]]
	then
		echo "Genome file is not valid, lengh of spacer sequences is 0! Check genome file"
		exit 0
	else
		echo "We got $ngg_length spacer sequences! Good! Let's move on!"
	fi


#add empty lines to sites file

if [ $debug = "T" ] || [ $debug = "t" ]
	then
		echo "Debug mode is true, cutting $all_ngg_sequences"
		head -n 50 $all_ngg_sequences > $all_ngg_sequences_dbg
	elif [[ $debug -eq "F" ]]
		then
			echo "Debug mod off! Do nothing!"
fi

echo "Adding empty lines to $all_ngg_sequences"

if [ $debug = "T" ] || [ $debug = "t" ]
	then
		awk ' {print;} NR % 1 == 0 { print ">"; }' $all_ngg_sequences_dbg > $all_ngg_sequences_space
	else
		awk ' {print;} NR % 1 == 0 { print ">"; }' $all_ngg_sequences > $all_ngg_sequences_space
fi

echo "Executing R script, adding fasta headers for each spacer and GC content calculation!"
echo "Some useless information may be printed below:"


./enter_fasta_headers.R $(pwd)

is_blastdb=$(ls | grep ".nin" | wc -l)

if [[ $is_blastdb -eq 0 ]]
	then
		echo "BLAST database not found! Creating..."
		makeblastdb -in $genome -dbtype nucl
	else
		echo "BLAST database found! Using existing..."
	fi

echo "Aligning spacer seqiences to reference genome! Evaluating XML blast output"
echo "XML blast"
blastn -task 'blastn-short' -db $genome -query $final_spacers -num_threads $threads -word_size $word_size -outfmt 5 -evalue 100 > blast.xml

echo "tabular blast"
blastn -task 'blastn-short' -db $genome -query $final_spacers -num_threads $threads -word_size $word_size -outfmt 6 -evalue 100 > blast.outfmt6

echo "Evaluating XML parser"
blastxmlparser --threads $threads -n 'hit.score, hsp.evalue, hsp.qseq, hsp.midline' blast.xml > blast.tsv

echo "Executing RNAfold!"
echo $final_spacers

RNAfold $final_spacers --noPS | grep ". (" | awk '{print $3}' | sed 's/)//' > energies.txt

echo "Executing final R script!"

./parse_tsv.R $(pwd) $annotation_file $prefix $threads

echo "Done!"
./purge.sh

echo "Converting to XLS! (ssconvert warning about X11 display is non-crucial, just skip it :) )"

ssconvert $prefix-result.csv $prefix-results.xls
