# pull_id_match.py
# - Parses through GopAga2.0 genome and pulls the 26 scaffolds with the most genes
#
# Requires: Biopython v1.7.3
# ---------------------------------------------------------------------------------

import re
import sys
from Bio import SeqIO

# Asks for the FASTA that will be parsed

FASTA = str(sys.argv[1])
output = str(sys.argv[3])
records = SeqIO.parse(FASTA, 'fasta')

# Scaffold IDs for the top 26 scaffolds with most genes

with open(sys.argv[2]) as f:
	patterns = f.readlines()
patterns = [x.strip() for x in patterns]
new_records = []
print(patterns)

# Iterates through each record and searches for scaffold IDs in the patterns array

for record in records:
	for pattern in patterns:
		if re.search(pattern, record.description) != None:
			new_records.append(record)

# Writes new FASTA file with selected scaffolds

SeqIO.write(new_records, output, 'fasta')

