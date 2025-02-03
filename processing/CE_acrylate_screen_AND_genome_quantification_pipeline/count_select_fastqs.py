import os
import csv

#fastq_files_dir = "/data/Unit_LMM/selberherr-group/strachan/stewart_metagenomes"
fastq_files_dir = "/data/Unit_LMM/selberherr-group/strachan/malmuthuge_metagenomes"

#csv_file = 'steward_sequence_counts.csv'
csv_file = 'malmuthuge_sequence_counts.csv'

# select_accessions = [
#     "ERR3220178_1",
#     "ERR3220179_1",
#     "ERR3220180_1",
#     "ERR3220181_1",
#     "ERR3220182_1",
#     "ERR3220183_1",
#     "ERR3220184_1",
#     "ERR3220185_1",
#     "ERR3220186_1",
#     "ERR3220187_1",
#     "ERR3220188_1",
#     "ERR3220189_1",
#     "ERR3220190_1",
#     "ERR3220191_1",
#     "ERR3220192_1",
#     "ERR3220193_1",
#     "ERR3220194_1",
#     "ERR3220195_1",
#     "ERR3220196_1",
#     "ERR3220197_1",
#     "ERR3220198_1"
# ]

select_accessions = [
    "SRR5190267_1",
    "SRR5190268_1",
    "SRR5190269_1",
    "SRR5190270_1",
    "SRR5190271_1",
    "SRR5190272_1",
    "SRR5190273_1",
    "SRR5190274_1",
    "SRR5190275_1",
    "SRR5190276_1",
    "SRR5190277_1",
    "SRR5190278_1",
    "SRR5190279_1",
    "SRR5190280_1",
    "SRR5190281_1",
    "SRR5190282_1",
    "SRR5190283_1",
    "SRR5190284_1"
]

# Prepare a list to hold the accession and sequence count
data = []

for accession in select_accessions:
    file_name = accession + ".fastq"
    file_path = os.path.join(fastq_files_dir, file_name)
    
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            # Count the total number of lines
            line_count = sum(1 for line in f)
            # Since each sequence in FASTQ format spans 4 lines
            sequence_count = line_count // 4
            data.append((accession, sequence_count))
    else:
        print(f"File {file_name} not found.")
        data.append((accession, 0))

# Write the data to a CSV file

with open(csv_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['accession', 'count'])
    for accession, count in data:
        writer.writerow([accession, count])

print(f"Sequence counts have been written to {csv_file}.")


