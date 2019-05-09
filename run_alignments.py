import re
import subprocess

category = input("Insert category to parse: ")
categories = []
while True:
	categories.append(category)
	category = input("Insert category to parse: ")
	if category == "":
		break

for cat in categories:
	cat_file = cat + ".fa"
	for i in range(2):
		if i == 1:
			cat_output = cat + "_alignment.fa"
			blastp_scout = subprocess.Popen(['./clustalw2', '-INFILE=' + cat_file, '-ALIGN', '-TYPE=PROTEIN', '-OUTPUT=FASTA', '-OUTFILE=' + cat_output])
			blastp_scout.wait()
		if i == 2:
			cat_output = cat + "_alignment.nxs"
			blastp_scout = subprocess.Popen(['./clustalw2', '-INFILE=' + cat_file, '-ALIGN', '-TYPE=PROTEIN', '-OUTPUT=NXS', '-OUTFILE=' + cat_output ])
			blastp_scout.wait()

