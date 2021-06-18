from Bio import SeqIO, Phylo
from Bio import AlignIO
from Bio import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os


def dot_lenght(input_file):
    records = SeqIO.parse("sequences/" + input_file + ".fasta", 'fasta')
    records = list(records)  # make a copy, otherwise our generator
    # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '-')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)

    # write to temporary file and do alignment
    output_file = 'sequences/{}_dot.fasta'.format(os.path.splitext(input_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')

def phylipTrascription(file_name,output_file):
    records = SeqIO.parse("sequences/" + file_name + ".fasta", "fasta")
    SeqIO.write(records, "sequences/" + output_file + ".phylip", "phylip")

dot_lenght("P2SequencesInd")
phylipTrascription("P2SequencesInd_dot","P2SequencesInd")
alignments = AlignIO.read("sequences/P2SequencesInd.phylip", "phylip")
distanceCalculator = DistanceCalculator('identity')
distanceMatrix = distanceCalculator.get_distance(alignments)

print('Distance Matrix:')
print(distanceMatrix)
print()

distanceConstructor = DistanceTreeConstructor()

NJTree = distanceConstructor.nj(distanceMatrix)
#Phylo.draw(NJTree)
print('Neighbor joining Phylogenetic Tree:')
Phylo.draw_ascii(NJTree)