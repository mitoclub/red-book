import os

import pandas as pd

# Build relative file path
dirname = os.path.dirname(__file__)
dirnameFromKG = os.path.join(dirname, '../../Body/1Raw/FromKG')

protParse = pd.read_csv(os.path.join(dirname, '../../Body/2Derived/codon_alignments_from_protparse.csv'), sep=";")
protParse = protParse.sort_values(by=['proteinName']).dropna()

# Files by gene
for gene in protParse['proteinName'].unique():
    print(gene)
    result = protParse.loc[protParse['proteinName'] == gene]
    #Open outfiles
    filename = gene + '.fasta'
    outfile = open(filename, 'wt')
    for species in result['speciesName']:
        alignmentProtein = result.loc[result['speciesName'] == species]['alignmentProtein'].values[0]
        outfile.write(">" + species + "\n" + alignmentProtein + "\n")
    outfile.close()
"""
outfile = open("codon_alignments_concatenated.fasta", 'wt')
count = 0

# File with all sequences concatenated
for species in protParse['speciesName'].unique():
    print(count, " ", species)
    count += 1
    # output fasta header
    header = ">" + species + "\n"
    if ["ND6"] in list(protParse.loc[(protParse['speciesName']==species)]['proteinName']):
        continue
    # By gene
    ND1 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ND1"), "alignmentProtein"].iloc[0]
    ND2 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ND2"), "alignmentProtein"].iloc[0]
    COX1 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="COX1"), "alignmentProtein"].iloc[0]
    ATP8 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ATP8"), "alignmentProtein"].iloc[0]
    ATP6 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ATP6"), "alignmentProtein"].iloc[0]
    COX3 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="COX3"), "alignmentProtein"].iloc[0]
    ND3 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ND3"), "alignmentProtein"].iloc[0]
    ND4L = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ND4L"), "alignmentProtein"].iloc[0]
    ND4 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ND4"), "alignmentProtein"].iloc[0]
    ND5 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ND5"), "alignmentProtein"].iloc[0]
    ND6 = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="ND6"), "alignmentProtein"].iloc[0]
    CYTB = protParse.loc[(protParse['speciesName']==species) & (protParse['proteinName']=="CYTB"), "alignmentProtein"].iloc[0]
    # concatenated sequences
    sequence = ND1 + ND2 + COX1 + ATP8 + ATP6 + COX3 + ND3 + ND4L + ND4 + ND5 + ND6 + CYTB + "\n"
    outfile.write(header)
    outfile.write(sequence)

outfile.close()
"""
