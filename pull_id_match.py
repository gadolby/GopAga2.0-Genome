import re
from Bio import SeqIO

FASTA = raw_input("What is the FASTA file to be searched? ")
records = SeqIO.parse(FASTA, 'fasta')
patterns = ["ScCC6lQ_161298_HRSCAF_192832", "ScCC6lQ_161421_HRSCAF_192955", "ScCC6lQ_16796_HRSCAF_38896", "ScCC6lQ_1923_HRSCAF_13295", "ScCC6lQ_140_HRSCAF_2235", "ScCC6lQ_24292_HRSCAF_48220", "ScCC6lQ_161394_HRSCAF_192928", "ScCC6lQ_40749_HRSCAF_67263", "ScCC6lQ_3764_HRSCAF_18001", "ScCC6lQ_161315_HRSCAF_192849"]
new_records = []

for record in records:
	for pattern in patterns:
		if re.search(pattern, record.description) != None:
			new_records.append(record)

SeqIO.write(new_records, "gaga2.0_10bestscaffolds.fa", 'fasta')

