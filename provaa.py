import os

from Bio import Seq, AlignIO, SeqIO


def dotLenght(input_file):

    records = SeqIO.parse("sequences/" + input_file + ".fasta", 'fasta')
    records = list(records)  # lunghezza righe dile input per avere un range
    maxlen = max(len(record.seq) for record in records)

    # inserisce il - per avere sequenze di stessa lunghezza
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '-')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)

    # scrive su file temporaneo e fa l'allineamento
    output_file = 'sequences/{}_dot.fasta'.format(os.path.splitext(input_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    return AlignIO.read("sequences/Proj2ItaSequences_dot.fasta","fasta")


print(dotLenght("Proj2ItaSequences"))


