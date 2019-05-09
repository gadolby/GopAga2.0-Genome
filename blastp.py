# blastp.py
# - Pulls gene family members of a specified gene family for a specified organism using reference sequences
#
# Requires : Python v3.1, Biopython v1.7.3
# -----------------------------------------------------------------------

import re
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Identify the BLAST database for the organism to be searched

organism_db = input("What organism database do you want to pull from? ")

results_file = organism_db + "_results.txt"
blast_results = open(results_file, 'w')

# Identify the gene family text file to be used. Gene family text files must be comma-separated lists of all gene family members.

gene_family = input("What gene family would you like to pull? ")
gene_family_file = open(gene_family + ".txt", 'r')
genes = gene_family_file.read().split(",")
del(genes[-1])

# Define function to remove duplicates from id list

def remove_duplicates(ids):
    newlist = []
    for value in ids:
       if value not in newlist:
           newlist.append(value)
    return newlist

organism_file = open(gene_family + '_' + organism + ".fa", 'w+')

for gene in genes:

	# Identify which files will contain the reference sequence and which will contain specified organism sequence

	ref_file = gene + "_ref.fasta"
	id_file = gene + "_" + organism
	protein_file = id_file + ".fa"

	# Scout out the top five blastp results
	
	blastp_scout = subprocess.Popen(['blastp', '-db', organism_db, '-query', ref_file, '-out', 'blast', '-num_descriptions', '5', '-num_alignments', '0'])
	blastp_scout.wait()

	blast = open('blast', 'r')
	print (blast.read())
	blast.seek(0)
	blast_results.write(blast.read())

	# Identify which result (1-5) seems to be a best fit for this blast search. If none seem like a good fit, skip with input '0'

	best_fit = input("Which result would you like to pull? Press 0 to skip gene. ")

	blast.close()

	if int(best_fit) != 0:

		# Gather top five result ids from blastp

		blastp_gather = subprocess.Popen(['blastp', '-db', organism_db, '-query', ref_file, '-outfmt', '6 sallacc', '-out', id_file, '-max_target_seqs', '5'])
		blastp_gather.wait()

		# Open the id file, split all of the ids, pull the desired id, and rewrite the file

		id_open = open(id_file, 'r+')
		id_read = id_open.read()
		print (id_read)
		id_split = id_read.split("\n")
		id_split = remove_duplicates(id_split)
		id_split = filter(None, id_split)
		id_split = list(id_split)
		print (id_split)
		id_open.close()
		id_open = open(id_file, 'w+')
		id_open.write(id_split[int(best_fit) - 1])
		id_open.close()
		id_read = ""
		id_split = []

		# Replace rewritten id with the FASTA file it directs to

		blastdbcmd = subprocess.Popen(['blastdbcmd', '-db', organism_db, '-entry_batch', id_file, '-outfmt', '%f', '-out', protein_file])
		blastdbcmd.wait()

		# Write to organism file for compilation

		parse_fasta = SeqIO.parse(protein_file, 'fasta')
		for record in parse_fasta:
			record.description = id_file
			record.id = ""
			SeqIO.write(record, "file.fa", 'fasta')
			record_open = open('file.fa', 'r')
			organism_file.write(record_open.read())




