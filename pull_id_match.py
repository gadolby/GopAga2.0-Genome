# pull_id_match.py
# - Parses through GopAga2.0 genome and pulls the 26 scaffolds with the most genes
# ---------------------------------------------------------------------------------

import re
from Bio import SeqIO

# Asks for the FASTA that will be parsed

FASTA = raw_input("What is the FASTA file to be searched? ")
records = SeqIO.parse(FASTA, 'fasta')

# Scaffold IDs for the top 26 scaffolds with most genes

patterns = ["ScCC6lQ_1650_HRSCAF_12279"
,"ScCC6lQ_3789_HRSCAF_18069"
,"ScCC6lQ_14958_HRSCAF_36502"
,"ScCC6lQ_161315_HRSCAF_192849"
,"ScCC6lQ_11901_HRSCAF_32349"
,"ScCC6lQ_3764_HRSCAF_18001"
,"ScCC6lQ_28138_HRSCAF_52803"
,"ScCC6lQ_24292_HRSCAF_48220"
,"ScCC6lQ_16796_HRSCAF_38896"
,"ScCC6lQ_161303_HRSCAF_192837"
,"ScCC6lQ_161323_HRSCAF_192857"
,"ScCC6lQ_97_HRSCAF_1492"
,"ScCC6lQ_140_HRSCAF_2235"
,"ScCC6lQ_40749_HRSCAF_67263"
,"ScCC6lQ_161298_HRSCAF_192832"
,"ScCC6lQ_161394_HRSCAF_192928"
,"ScCC6lQ_3500_HRSCAF_17357"
,"ScCC6lQ_38465_HRSCAF_64658"
,"ScCC6lQ_1923_HRSCAF_13295"
,"ScCC6lQ_161370_HRSCAF_192904"
,"ScCC6lQ_161354_HRSCAF_192888"
,"ScCC6lQ_161421_HRSCAF_192955"
,"ScCC6lQ_9220_HRSCAF_28431"
,"ScCC6lQ_1442_HRSCAF_11293"
,"ScCC6lQ_13971_HRSCAF_35176"
,"ScCC6lQ_1417_HRSCAF_11161"]
new_records = []

# Iterates through each record and searches for scaffold IDs in the patterns array

for record in records:
	for pattern in patterns:
		if re.search(pattern, record.description) != None:
			new_records.append(record)

# Writes new FASTA file with selected scaffolds

SeqIO.write(new_records, "gaga2.0_L50scaffolds.fa", 'fasta')

