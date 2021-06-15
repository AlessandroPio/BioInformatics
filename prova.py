from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIWWW

dati = SeqIO.read('sequences.fasta','fasta')
print(type(dati))
result = NCBIWWW.qblast("blastn","nt", dati.seq)
f_result = open("re.xml","w")
f_result.write(result.read())
f_result.close()
result.close()