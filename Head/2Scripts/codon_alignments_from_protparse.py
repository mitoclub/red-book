import os
from glob import glob
import lzma

import re
from Bio import SeqIO
import pandas as pd

# Build relative file path
dirname = os.path.dirname(__file__)
dirnameFromKG = os.path.join(dirname, '../../Body/1Raw/FromKG')

protParse = pd.read_csv(os.path.join(dirname, '../../Body/2Derived/protparse_to_csv.csv'), sep=";")
genbank1 = os.path.join(dirname, '../../Body/1Raw/FromKG/mitochondrion.1.genomic.tar.xz')
genbank2 = os.path.join(dirname, '../../Body/1Raw/FromKG/mitochondrion.2.genomic.tar.xz')

# Open protparse CSV
file3 = os.path.join(dirname, '../../Body/2Derived/protparse_to_csv.csv')
protparse = pd.read_csv(file3, sep=";")

header = ["speciesName", "proteinName", "alignmentProtein"]
outfile = open(os.path.join(dirname, '../../Body/2Derived/codon_alignments_from_protparse.csv'), 'wt')
outfile.write(";".join(header) + "\n")

# Unarchive and read file
unarchived1 = lzma.open(genbank1, mode='rt')
unarchived2 = lzma.open(genbank2, mode='rt')
count = 0
# Get nucleotide and protein sequence
for rec in SeqIO.parse(unarchived1, "genbank"):
    id = rec.name # NCBI id
    speciesName = "_".join(rec.annotations['organism'].split(" "))
    # validate name
    if speciesName in list(protparse['speciesName']):
        result = protParse.loc[protParse['speciesName'] == speciesName]
        for feature in rec.features:
            if feature.type == "CDS":
                    # Protein Name
                if feature.qualifiers.get('gene') is not None:
                    featureName = feature.qualifiers.get('gene')[0].upper()
                else:
                    continue # No feature name

                featureProt = feature.qualifiers['translation'][0] # protein sequence
                featureSeq = str(feature.extract(rec.seq)) # dna sequence

                if featureName in result['proteinName'].to_list():
                    alignmentProtein = result.loc[result['proteinName'] == str(featureName)]['alignmentProtein'].values.item()
                    alignmentProteinList = list(result.loc[result['proteinName'] == str(featureName)]['alignmentProtein'].values.item())
                        # validate prot sequence
                    if alignmentProtein.replace("-", "") == featureProt:
                        codonAlign = ""
                        while len(alignmentProteinList) > 0:
                            if alignmentProteinList[0] == "-":
                                codonAlign += "---"
                                alignmentProteinList = alignmentProteinList[1::]
                            else:
                                if len(featureSeq)<3:
                                    print(speciesName)
                                codon = featureSeq[0:3]
                                codonAlign += codon
                                featureSeq = featureSeq[3::]
                                alignmentProteinList = alignmentProteinList[1::]
                        #write to df...
                        if len(alignmentProtein)*3 == len(codonAlign):
                            outfile.write(";".join([speciesName, featureName, codonAlign]) + "\n")
unarchived1.close()
# do the same for unarchived2
for rec in SeqIO.parse(unarchived2, "genbank"):
    id = rec.name # NCBI id
    speciesName = "_".join(rec.annotations['organism'].split(" "))
    # validate name
    if speciesName in list(protparse['speciesName']):
        result = protParse.loc[protParse['speciesName'] == speciesName]
        for feature in rec.features:
            if feature.type == "CDS":
                    # Protein Name
                if feature.qualifiers.get('gene') is not None:
                    featureName = feature.qualifiers.get('gene')[0].upper()
                else:
                    continue # No feature name

                featureProt = feature.qualifiers['translation'][0] # protein sequence
                featureSeq = str(feature.extract(rec.seq)) # dna sequence

                if featureName in result['proteinName'].to_list():
                    alignmentProtein = result.loc[result['proteinName'] == str(featureName)]['alignmentProtein'].values.item()
                    alignmentProteinList = list(result.loc[result['proteinName'] == str(featureName)]['alignmentProtein'].values.item())
                        # validate prot sequence
                    if alignmentProtein.replace("-", "") == featureProt:
                        codonAlign = ""
                        while len(alignmentProteinList) > 0:
                            if alignmentProteinList[0] == "-":
                                codonAlign += "---"
                                alignmentProteinList = alignmentProteinList[1::]
                            else:
                                if len(featureSeq)<3:
                                    print(speciesName)
                                codon = featureSeq[0:3]
                                codonAlign += codon
                                featureSeq = featureSeq[3::]
                                alignmentProteinList = alignmentProteinList[1::]
                        #write to df...
                        if len(alignmentProtein)*3 == len(codonAlign):
                            outfile.write(";".join([speciesName, featureName, codonAlign]) + "\n")
unarchived2.close()
outfile.close()
