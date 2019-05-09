from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import cycle

# Initialize all variables and ask for input of FASTA and categories

categories = []
fasta = input("Insert name of FASTA file to categorize: ")
category = input("Insert category to parse: ")
new_fasta = category + ".fa"
category = category + "_"

# Append all categories desired by user until nothing is input

while True:
	categories.append(category)
	category = input("Insert category to parse: ")
	if category == "":
		break
	category = category + "_"

records = SeqIO.parse(fasta, 'fasta')
cat_records = []

# For all categories, search the string in each record and pull records that match

for cat in categories:
	for record in records:
		if record.id.find(cat) != -1:
			print(record.id)
			cat_records.append(record)
	records = SeqIO.parse(fasta, 'fasta')

# Write the FASTA with all records pertaining to the categories

SeqIO.write(cat_records, new_fasta, 'fasta')
