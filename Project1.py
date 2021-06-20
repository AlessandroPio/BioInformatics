from __future__ import division

from Bio import SeqIO, Phylo
from Bio import AlignIO
from Bio import Seq
from Bio.SeqUtils import ProtParam
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os, glob
import pandas as pd

def writeToFileTree(tree, output_file):
    Phylo.write(tree, "sequences/" + output_file + "_phylotree.xml", "phyloxml")

def writeToFile(content,file):
    file.write(content + "\n")

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
    os.remove("sequences/"+ input_file + "_dot.fasta")
    alignments = AlignIO.read("sequences/" + input_file + ".phylip", "phylip")
    distanceCalculator = DistanceCalculator('identity')
    distanceMatrix = distanceCalculator.get_distance(alignments)
    distanceConstructor = DistanceTreeConstructor()
    tree = distanceConstructor.nj(distanceMatrix)
    tree.clade[0, 1].color = "blue"
    writeToFileTree(tree, input_file)
    return tree

def genomeAnalysis(file):
    def conv(item):
        return len(item)

        def to_str(item):
            return str(item)
            df['sequence_str'] = df[0].apply(to_str)

    def getProtein():
        X = ProtParam.ProteinAnalysis(str(record))
        protein_of_interest = X.count_amino_acids()
        p_o_i_list.append(protein_of_interest)
        molecular_weight = X.molecular_weight()
        molecular_weight_list.append(molecular_weight)

        writeToFile("PROTEIN OF INTEREST         N° -> " + str(protein_of_interest), file_output)
        writeToFile("AMINO ACIDS PERCENT         N° -> " + str(X.get_amino_acids_percent()), file_output)
        writeToFile("MOLECULAR WEIGHT            N° -> " + str(X.aromaticity()), file_output)
        writeToFile("FLEXIBILITY                 N° -> " + str(X.flexibility()), file_output);

        #print("Aromaticity = ", X.aromaticity())
        #print("Isoelectric point = ", X.isoelectric_point())
        #print("Secondary structure fraction = ", X.secondary_structure_fraction())

    file_output = open("sequences/" + input_file + "_phylogenetic_analysis.txt", "w+")
    cont = 0
    for sequence in SeqIO.parse('sequences/' + file + ".fasta", "fasta"):

        dnaSequence = sequence.seq
        nucleotide_of_sequence = len(sequence)
        mRNA = dnaSequence.transcribe()  # Transcribe a DNA sequence into RNA.
        mRna_lenght = len(mRNA)

        Amino_Acid = mRNA.translate(table=1, cds=False)
        #print('Amino Acid', Amino_Acid)
        #print("Length of Protein:", len(Amino_Acid))
        #print("Length of Original mRNA:", len(mRNA))

        Proteins = Amino_Acid.split('*')  # * is translated stop codon
        df = pd.DataFrame(Proteins)
        df.describe()
        total_protein = len(df)

        cont += 1
        writeToFile("SEQUENCE           N° -> " + str(sequence.id) + " (" + str(cont) + ")", file_output);writeToFile("-----------------------------------", file_output)
        writeToFile("NUCLEOTIDE         N° -> " + str(nucleotide_of_sequence), file_output)
        writeToFile("DNA                N° -> " + str(len(dnaSequence)), file_output)
        writeToFile("mRNA               N° -> " + str(mRna_lenght), file_output)
        writeToFile("PROTEIN            N° -> " + str(total_protein), file_output);writeToFile("-----------------------------------", file_output)
        writeToFile("mRNA SEQUENCE         -> " + str(mRNA), file_output)
        writeToFile("DNA  SEQUENCE         -> " + str(dnaSequence), file_output)
        writeToFile("AMINO ACID            -> " + str(Amino_Acid), file_output);writeToFile("-----------------------------------", file_output)

        p_o_i_list = []
        molecular_weight_list = []

        for record in Proteins[:]:
            try:
                getProtein()
            except:
                continue

    file_output.close()

       # df['length'] = df[0].apply(conv)
       # df.rename(columns={0: "sequence"}, inplace=True)
       # df.head()  # Take only longer than 20
       # functional_proteins = df.loc[df['length'] >= 20]
       # print('Total functional proteins:', len(functional_proteins))
       # functional_proteins.describe()

while True:

    directory = glob.glob("sequences/*.fasta");
    print("----> File presenti nella directory");
    print("|")
    for file in directory:
        file = file.split("/")[1].split('.')[0];
        print("| " + file)
    print("|");
    print("----")
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


tree = phylogeneticTree(phylipTrascription(input_file))  # creo l'albero filogenetico
genomeAnalysis(input_file)                        # parto con l'analisi genomica delle sequenze

msgbox = input("| Visualizzare L'albero filogenetico delle sequenze? (Y/N) -> ")
if msgbox.upper() == "Y":
    Phylo.draw(tree)

print("| Processo terminato!")