from Bio.Seq import Seq
from Bio import SeqIO
my_seq = Seq("ACGT")                            # crea un oggetto sequenza

compl = my_seq.complement()                     # complemeto della sequenza (l'opposto)
compl_reverse = my_seq.reverse_complement()     # non capisco perche non rida il valore di myseq


# per la lettura di un file andremo ad utilizzare la classe SeqIO che ha le funzionalità per farlo
# Un file di tipo FASTA inizia la sequenza con un carattere speciale '>' seguito dalle sequenze

for seq_record in SeqIO.parse('nome_file','estensione'):    # es SeqIO.parse('ls_orchid.fasta','fasta')
    print(seq_record.campo)                                 # es seq_record.id / seq_record.seq / len(seq_record)

# il formato FASTA non specifica l'alfabeto quindi Bio.Seq è impostato in maniera
# predefinita su SingleLetterAlphabet() piuttosto generico piuttosto che su qualcosa di specifico del DNA

# Nel caso del file genbank, Bio.SeqIO ha potuto scegliere un alfabeto sensato, IUPACAmbiguousDNA(). Noterai anche
# che in questo caso è stata utilizzata una stringa più corta come seq_record.id

# Entrambi gli esempi di sopra sono un output di seq_record.id

# i tipi di formato possono essere
#
#     1° param         2° param
#   - .fasta     ->     fasta
#   - .gbk       ->     genbank

print(compl)
print(compl_reverse)

#SEQUENCE OBJECT

# Ci sono 2 importanti differenze importante tra Oggetti Seq e Stringhe python
# perche hanno differenti metodi (funzioni)

# 1)
# Sebbene l'oggetto Seq supporti molti degli stessi metodi di una stringa semplice, il metodo
# translate() differisce facendo traduzione biologica, e ci sono anche altri metodi biologicamente rilevanti
# come reverse_complement()

# 2)
# l'oggetto Seq ha un attributo importante, l'alfabeto, che è un oggetto che descrive cosa “significano”
# i singoli caratteri che compongono la stringa di sequenza e come dovrebbero essere interpretati.

# Gli alfabeti attualmente disponibili per Biopython sono definiti nel modulo Bio.Alphabet.
# Useremo qui gli alfabeti IUPAC per trattare: DNA, RNA e proteine.
# Bio.Alphabet.IUPAC fornisce definizioni di base per proteine, DNA e RNA, ma offre inoltre la possibilità
# di estendere e personalizzare le definizioni di base. Ad esempio, per le proteine, esiste una classe
# base IUPACProtein