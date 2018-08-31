#!/bin/bash
touch data.tsv
	for f in $(seq 1 30);
		do
			df -k | grep "raid" >> data.tsv
			sleep 1
			done
			./visual.R

rm data.tsv
