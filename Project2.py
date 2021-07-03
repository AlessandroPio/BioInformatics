import os, platform
import time
from itertools import combinations

from Bio import SeqIO, Seq, AlignIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

os.chdir("sequences/Project1/")

if(platform.system() != 'Windows'): #type
    os.system("cat Proj1VarIndia.fasta Proj1VarItalia.fasta | awk '/^>/{f=!d[$1];d[$1]=1}f' > sequences.fasta")

def align(info):
    try:
        alignments = aligner.align(info[0].seq, info[1].seq)
        alignment.append(SeqRecord(Seq(str(alignments[0])), id=info[0].id, name=info[1].name, description=f"Score={str(alignments.score)}"))
        alignment.append(SeqRecord(Seq(str(alignments[1])), id=info[1].id, name=info[1].name, description=f"Score={str(alignments.score)}"))

    except(IndexError):
        pass

    #SeqIO.write(informations, "output.fasta", "fasta")

    #print("| Score                   -> " + str(alignments.score))
    #print("|")
    #print("| Id Sequence             -> " + str(info[0].id))
    #print("| Mutation Detected       -> " + str((len(info[0].seq) - int(alignments.score))))
    #print("| Lenght of sequence      -> " + str(len(info[0].seq)))
    #print("| ----------------------------")
    #print("| Id Sequence             -> " + str(info[1].id))
    #print("| Mutation Detected       -> " + str((len(info[1].seq) - int(alignments.score))))
    #print("| Lenght of sequence      -> " + str(len(info[1].seq)))

aligner = Align.PairwiseAligner(match_score=1.0)

alignment = [];cont=1
init_time = time.localtime()
print("| Start At  -> " + str(init_time.tm_hour) + ":" + str(init_time.tm_min) + ":" + str(init_time.tm_sec))
for combination in combinations(SeqIO.parse("sequences.fasta", "fasta"),2):
    print(cont)
    cont+=1
    #print("-----------------------------------")
    align(combination)
    #print("-----------------------------------")
end_time = time.localtime()
print("| End At    -> " + str(end_time.tm_hour) + ":" + str(end_time.tm_min) + ":" + str(end_time.tm_sec))

print();print()
print(len(alignment))
