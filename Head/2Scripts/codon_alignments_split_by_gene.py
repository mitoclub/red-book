import os

import pandas as pd

# Build relative file path
dirname = os.path.dirname(__file__)
dirnameFromKG = os.path.join(dirname, '../../Body/1Raw/FromKG')

protParse = pd.read_csv(os.path.join(dirname, '../../Body/2Derived/codon_alignments_from_protparse.csv'), sep=";")
protParse = protParse.sort_values(by=['proteinName']).dropna()

# sort dataframe by gene
for gene in protParse['proteinName'].unique():
    print(gene)
    result = protParse.loc[protParse['proteinName'] == gene]
    #Open outfile
    filename = gene + '.fasta'
    outfile = open(filename, 'wt')
    for species in result['speciesName']:
        alignmentProtein = result.loc[result['speciesName'] == species]['alignmentProtein'].values[0]
        outfile.write(">" + species + "\n" + alignmentProtein + "\n")
    outfile.close()
