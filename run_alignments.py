# run_alignments.py
# - Takes protein FASTA files and creates FASTA and NEXUS alignments using ClustalW2
# 
# Requires : ClustalW2 v2.1
# --------------------------------------------------------------------

import re
import subprocess

# Asks user to input FASTA file names (without suffix) to be aligned (separately) until nothing is input

category = input("Insert category to parse: ")
categories = []
while True:
	categories.append(category)
	category = input("Insert category to parse: ")
	if category == "":
		break

# Takes each FASTA file and aligns using ClustalW2 with FASTA output then NEXUS output

for cat in categories:
	cat_file = cat + ".fa"
	for i in range(2):
		if i == 1:
			cat_output = cat + "_alignment.fa"
			blastp_scout = subprocess.Popen(['clustalw2', '-INFILE=' + cat_file, '-ALIGN', '-TYPE=PROTEIN', '-OUTPUT=FASTA', '-OUTFILE=' + cat_output])
			blastp_scout.wait()
		if i == 2:
			cat_output = cat + "_alignment.nxs"
			blastp_scout = subprocess.Popen(['clustalw2', '-INFILE=' + cat_file, '-ALIGN', '-TYPE=PROTEIN', '-OUTPUT=NXS', '-OUTFILE=' + cat_output ])
			blastp_scout.wait()

