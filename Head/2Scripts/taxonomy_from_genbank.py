import os
import lzma
from Bio import SeqIO

# Build relative file path
dirname = os.path.dirname(__file__)
file1 = os.path.join(dirname, '../../Body/1Raw/FromKG/mitochondrion.1.genomic.tar.xz')
file2 = os.path.join(dirname, '../../Body/1Raw/FromKG/mitochondrion.2.genomic.tar.xz')

# Open file for writing
outfile = os.path.join(dirname, '../../Body/2Derived/taxonomy_from_KG_genbank.csv')

# Unarchive and read file
unarchived1 = lzma.open(file1, mode='rt')
unarchived2 = lzma.open(file2, mode='rt')

def clean(desc):
    if "mitochondrion, complete genome" in desc:
        desc = desc.replace("mitochondrion, complete genome", "")
    if "mitochondrion,complete genome" in desc:
        desc = desc.replace("mitochondrion,complete genome", "")
    if "mitochondrion complete genome" in desc:
        desc = desc.replace("mitochondrion complete genome", "")
    if "mitochondria, complete genome" in desc:
        desc = desc.replace("mitochondria, complete genome", "")
    if "complete mitochondrial genome" in desc:
        desc = desc.replace("complete mitochondrial genome", "")
    if "mitochondrial complete genome" in desc:
        desc = desc.replace("mitochondrial complete genome", "")
    if "mitochondrial DNA, complete genome" in desc:
        desc = desc.replace("mitochondrial DNA, complete genome", "")
    if "mitochondrial DNA, complete sequence" in desc:
        desc = desc.replace("mitochondrial DNA, complete sequence", "")
    if "voucher" in desc or "strain" in desc or "culture-collection" in desc or "isolate" in desc or "haplotype" in desc:
        desc = " ".join(desc.split(" ")[0:2])
    return desc

with open(outfile, "wt") as output:
    output.write("SpeciesName,GeneBankID,FullTaxonomy\n")
    for rec in SeqIO.parse(unarchived1, "genbank"):
        # NCBI id
        id = rec.name
        taxonomy = ">".join(rec.annotations['taxonomy'])
        name = clean(rec.description)
        output.write(",".join([name, id, taxonomy]) + '\n')
    unarchived1.close()
    for rec in SeqIO.parse(unarchived2, "genbank"):
        # NCBI id
        id = rec.name
        taxonomy = ">".join(rec.annotations['taxonomy'])
        name = clean(rec.description)
        output.write(",".join([name, id, taxonomy]) + '\n')
    unarchived2.close()
