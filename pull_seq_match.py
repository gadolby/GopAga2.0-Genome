import re
from Bio import SeqIO

FASTA = raw_input("What is the FASTA file to be searched? ")
records = SeqIO.parse(FASTA, 'fasta')
patterns = ["augustus-Contig5-processed-gene-0.105-mRNA-1"
,"augustus-Contig1846-processed-gene-0.11-mRNA-1"
,"augustus-Contig1846-processed-gene-0.12-mRNA-1"
,"augustus-Contig1846-processed-gene-0.10-mRNA-1"
,"augustus-Contig356-processed-gene-0.0-mRNA-1"
,"augustus-Contig356-processed-gene-0.1-mRNA-1"]
new_records = []

for record in records:
	for pattern in patterns:
		if re.search(pattern, record.description) != None:
			new_records.append(record)

SeqIO.write(new_records, "galapagos_TLR7", 'fasta')

