import os
import pandas as pd
import lzma
from Bio import SeqIO

# Build relative file path
dirname = os.path.dirname(__file__)

# Open Genbank
file1 = os.path.join(dirname, '../../Body/1Raw/FromKG/mitochondrion.1.genomic.tar.xz')
file2 = os.path.join(dirname, '../../Body/1Raw/FromKG/mitochondrion.2.genomic.tar.xz')

unarchived1 = lzma.open(file1, mode='rt')
unarchived2 = lzma.open(file2, mode='rt')


# Open protparse CSV
file3 = os.path.join(dirname, '../../Body/2Derived/protparse_to_csv.csv')

protparse = pd.read_csv(file3, sep=";")


# Open file for writing
outfile = os.path.join(dirname, '../../Body/2Derived/taxonomy_from_KG_genbank.csv')


with open(outfile, "wt") as output:
    output.write("SpeciesName,GeneBankID,FullTaxonomy\n")
    for rec in SeqIO.parse(unarchived1, "genbank"):
        speciesName = "_".join(rec.annotations['organism'].split(" "))
        geneId = rec.name
        taxonomy = ">".join(rec.annotations['taxonomy'])
        if speciesName in list(protparse['speciesName']):
            output.write(",".join([speciesName, geneId, taxonomy]) + '\n')
    unarchived1.close()
    for rec in SeqIO.parse(unarchived2, "genbank"):
        speciesName = "_".join(rec.annotations['organism'].split(" "))
        geneId = rec.name
        taxonomy = ">".join(rec.annotations['taxonomy'])
        if speciesName in list(protparse['speciesName']):
            output.write(",".join([speciesName, geneId, taxonomy]) + '\n')
    unarchived2.close()
