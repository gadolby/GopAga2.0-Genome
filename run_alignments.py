# run_alignments.py
# - Takes protein FASTA files and creates FASTA and NEXUS alignments using ClustalW2
# 
# Requires : ClustalW2 v2.1
# --------------------------------------------------------------------

import re
import sys
import subprocess

# Asks user to input FASTA file names (without suffix) to be aligned (separately) until nothing is input

with open(sys.argv[1]) as f:
	categories = f.readlines()
categories = [x.strip() for x in categories]

# Takes each FASTA file and aligns using ClustalW2 with FASTA output then NEXUS output

for cat in categories:
	cat_file = cat + ".fa"
	print(cat_file)
	for i in range(3):
		if i == 0:
			cat_output = cat + "_alignment.fa"
			clustal1 = subprocess.Popen(['clustalw2', '-INFILE=' + cat_file, '-ALIGN', '-TYPE=PROTEIN', '-MATRIX=GONNET', '-OUTPUT=FASTA', '-OUTFILE=' + cat_output])
			clustal1.wait()
		if i == 1:
			cat_output = cat + "_alignment.nxs"
			clustal2 = subprocess.Popen(['clustalw2', '-INFILE=' + cat_file, '-ALIGN', '-TYPE=PROTEIN', '-MATRIX=GONNET', '-OUTPUT=NXS', '-OUTFILE=' + cat_output ])
			clustal2.wait()
		if i == 2:
			cat_output = cat + "_alignment.phy"
			clustal3 = subprocess.Popen(['clustalw2', '-INFILE=' + cat_file, '-ALIGN', '-TYPE=PROTEIN', '-MATRIX=GONNET', '-OUTPUT=PHYLIP', '-OUTFILE=' + cat_output ])
			clustal3.wait()

