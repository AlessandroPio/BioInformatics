import os, platform
from itertools import combinations

from Bio import SeqIO, Seq, AlignIO, pairwise2
from Bio.SeqRecord import SeqRecord

os.chdir("sequences/Project1/")

if(platform.system() == 'Windows'):
    os.system("type Proj1VarIndia.fasta Proj1VarItalia.fasta | awk '/^>/{f=!d[$1];d[$1]=1}f' > sequences.fasta")
else:
    os.system("cat Proj1VarIndia.fasta Proj1VarItalia.fasta | awk '/^>/{f=!d[$1];d[$1]=1}f' > sequences.fasta")

def align(pair):
    al = pairwise2.align.globalms(pair[0].seq, pair[1].seq, 2, -1, -5, -2)[0]

    try:
        yield SeqRecord(Seq(al[0]), id=pair[0].id, name=pair[0].name, description=f"Score={al[2]}")
        yield SeqRecord(Seq(al[1]), id=pair[1].id, name=pair[1].name, description=f"Score={al[2]}")
    except(Exception):
        pass

def getMutation(file):
    seq_records = AlignIO.read(file + '.phylip', 'phylip')
    print("-----------------------------------")
    j = 0
    while j < len(seq_records):
        y = 0
        location = []
        cont = 0
        i = 0
        while i < len(seq_records):
            if(j != i):
                while y < len(seq_records[i]):
                    if(seq_records[i].seq[y] != seq_records[j].seq[y]):
                        cont += 1
                        location.append(y)
                    y += 1
            i += 1
        if(cont != len(seq_records[j])):
            print("Mutation Detected!   (" + str(j) + ")")
            print("Sequences Id         ->", seq_records[j].id)
            print("Lenght of sequence   ->", len(seq_records[j]))
            print("Number of mutations  ->", len(location))
            print("-----------------------------------")
        j += 1

with open("output.fasta", "w") as output:
    for pair in combinations(SeqIO.parse('sequences.fasta', "fasta"), 2):
        SeqIO.write(align(pair), output, "fasta")
#phylipTrascription(SeqIO.parse(infile, "fasta"), 2)
#phylipTrascription('sequences')
#getMutation('sequences')