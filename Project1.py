from __future__ import division

from Bio import SeqIO, Phylo
from Bio import AlignIO
from Bio import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os, glob

def writeToFile(tree, input_file):
    Phylo.write(tree, "sequences/" + input_file + '_phylogeny.txt', "phyloxml")

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

def phylipTrascription(file_name):
    records = SeqIO.parse("sequences/" + file_name + "_dot.fasta", "fasta")
    SeqIO.write(records, "sequences/" + file_name + ".phylip", "phylip")

def phylogeneticTree(alignments):
    print("|")
    alignments = AlignIO.read("sequences/" + input_file + ".phylip", "phylip")
    distanceCalculator = DistanceCalculator('identity')
    distanceMatrix = distanceCalculator.get_distance(alignments)
    distanceConstructor = DistanceTreeConstructor()
    tree = distanceConstructor.nj(distanceMatrix)
    #tree.clade[0, 1].color = "blue"
    writeToFile(tree, input_file)

def genomeAnalysis(file):
    for sequence in SeqIO.parse('sequences/' + file + ".fasta", "fasta"):
        print(sequence.seq)
        print(len(sequence), 'nucliotides')
    dnaSequence = SeqIO.read('sequences/' + file + ".fasta", "fasta")
    dna = dnaSequence.seq  # Convert DNA into mRNA Sequence
    mRNA = dna.transcribe()  # Transcribe a DNA sequence into RNA.
    print(mRNA)
    print('Size : ', len(mRNA))
    Amino_Acid = mRNA.translate(table=1, cds=False)
    print('Amino Acid', Amino_Acid)
    print("Length of Protein:", len(Amino_Acid))
    print("Length of Original mRNA:", len(mRNA))

    Proteins = Amino_Acid.split('*')  # * is translated stop codon
    df = pd.DataFrame(Proteins)
    df.describe()
    print('Total proteins:', len(df))

    def conv(item):
        return len(item)

        def to_str(item):
            return str(item)
            df['sequence_str'] = df[0].apply(to_str)

    poi_list = []
    MW_list = []
    from Bio.SeqUtils import ProtParam
    for record in Proteins[:]:
        print("\n")
        X = ProtParam.ProteinAnalysis(str(record))
        POI = X.count_amino_acids()
        poi_list.append(POI)
        MW = X.molecular_weight()
        MW_list.append(MW)
        print("Protein of Interest = ", POI)
        print("Amino acids percent =    ", str(X.get_amino_acids_percent()))
        print("Molecular weight = ", MW_list)
        print("Aromaticity = ", X.aromaticity())
        print("Flexibility = ", X.flexibility())
        print("Isoelectric point = ", X.isoelectric_point())
        print("Secondary structure fraction = ", X.secondary_structure_fraction())

    df['length'] = df[0].apply(conv)
    df.rename(columns={0: "sequence"}, inplace=True)
    df.head()  # Take only longer than 20
    functional_proteins = df.loc[df['length'] >= 20]
    print('Total functional proteins:', len(functional_proteins))
    functional_proteins.describe()

while True:

    print("|")
    input_file = str(input("| File di riferimento -> "))
    input_file = input_file.split('.')[0]
    print("|")

    try:
        print("| Controllo del file...")
        dotLenght(input_file)
        print("| Controllo terminato correttamente")
        break

    except IOError:
        print("| File non valido!")
        print("| Immettere il nome di un file formato fasta situato in 'sequences/'")
        print("|")

        directory = glob.glob("sequences/*.fasta");
        print("----> File presenti nella directory");print("|")
        for file in directory:
            file = file.split("/")[1].split('.')[0];
            print("| " + file)
        print("|");print("----")

phylogeneticTree(phylipTrascription(input_file))
genomeAnalysis(input_file)

print("| Processo terminato!")