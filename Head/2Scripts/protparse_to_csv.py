import os

import re
from Bio import SeqIO

# Build relative file path
dirname = os.path.dirname(__file__)
dirnameFromKG = os.path.join(dirname, '../../Body/1Raw/FromKG')

# List of protein alignment files
protparseFiles = [file for file in os.listdir(dirnameFromKG) if file.endswith('final.fas')]

# Open file for writing
outfile = open(os.path.join(dirname, '../../Body/2Derived/protparse_to_csv.csv'), 'wt')

# Create dataframe
header = ["speciesName", "proteinName", "alignmentProtein"]
outfile.write(";".join(header) + "\n")

# Read protsequences into memory and then dataframe
count = 0

for filename in protparseFiles:
    fastaFile = SeqIO.parse(open("/".join([dirnameFromKG, filename])),'fasta')

    proteinName = re.search(r'sel\.(.*?).promals', filename).group(1)
    print(proteinName)
    for fasta in fastaFile:
        name, sequence = fasta.id, str(fasta.seq)
        outfile.write(";".join([name, proteinName, sequence]) + "\n")
        count += 1

outfile.close()
