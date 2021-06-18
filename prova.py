from Bio import AlignIO

with open("sequences/P1Sequences.fasta") as handle:
    records = AlignIO.parse(handle, "fasta")

    with open("sequences/P1Sequences.phylip", "w") as output_handle:
        AlignIO.write(records, output_handle, "phylip")
