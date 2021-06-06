# Per una sequenza nucleotidica, si può ottenere il complemento o l'inverso del complemento
# di una sequenza attraverso metodi built-in

# RIASSUNTO
# from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
# my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"   -> DEFINISCE LA SEQUENZA
# reverse_complement(my_string)                              -> FA IL COMPLEMENTO E INVERTE LA SEQUENZA
# transcribe(my_string)                                      -> RICAVA IL mRNA
# back_transcribe(my_string)                                 -> RICAVA IL DNA (?)
# translate(my_string)                                       -> RICAVA LA SEQUENZA NUCLEOTIDICA



from Bio.Seq import Seq # help

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC") #, IUPAC.unambiguous_dna
print("| La sequenza nucleotidica è         -> " + my_seq)
print("| Il  complemento della sequenza è   -> " + my_seq.complement())
print("| Mentre il complemento inverso è    -> " + my_seq.reverse_complement())
print("| Il mRNA della seq nucleotidica è   -> " + str(my_seq).replace("T","U")) #my_seq.reverse_complement().transcribe()
print("| L'ultima sequenza fa riferimento al mRNA della sequenza nucleotidica")
print("| In maniera alternativa possiamo scrivere my_seq[::-1]")
print("-------------------------------------------------------------------------------")
print("|")

protein_seq = Seq("EVRNAK")  #, IUPAC.protein
print("| La sequenza proteica è -> "+protein_seq)
print("| NON SI PUO FARE IL COMPLEMENTO DELLA SEQUENZA PROTEICA")
print("-------------------------------------------------------------------------------")
print("|")

messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGC") #, IUPAC.unambiguous_rna
print("FASE DI TRASCRIZIONE")
print("| La sequenza dell'mRNA è ->    " + messenger_rna)
print("| La sequenza nucleotidica è -> " + messenger_rna.back_transcribe())
print("-------------------------------------------------------------------------------")
print("|")

messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG") #, IUPAC.unambiguous_rna
print("FASE DI TRADUZIONE")
# Possiamo fare la traduzione in proteina partendo dalla sequenza dell'rna messaggiero
print("| La sequenza dell'mRNA è                                -> " + messenger_rna)
print("| La proteinca corrispondente alla sequenza del mRNA è   -> " + messenger_rna.translate())
print("|")
# Possiamo fare la traduziona in proteina partendo anche dalla sequenza di dna
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG") #, IUPAC.unambiguous_rna
print("| La sequenza del DNA è                                  -> " + coding_dna)
print("| La proteina corrispondente alla sequenza di dna è      -> " + coding_dna.translate())
# Da notare l'asterisco che indica il codone di stop
print("-------------------------------------------------------------------------------")
print("|")

print("TABELLA CODONI NCBI") # Tabella molto piu leggibile
from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]               # Tabella di traduzione standard
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]   # Tabella di traduzione per il DNA mitocondriale dei vertebrati
# in maniera alternativa queste 2 tabelle sono etichettate con id 1 e 2
#   ->  standard_table = CodonTable.unambiguous_dna_by_id[1]
#   ->  mito_table = CodonTable.unambiguous_dna_by_id[2]
print(standard_table)
print(mito_table)
print("|")

# Confrontare 2 sequenze non è proprio semplice perche prendendo in considerazione la sola base azotata A
# potrebbe riferirsi ad una sequenza di dna, rna o anche proteina
# Per questo BioPython utilizza gli alfabeti per specificare il 'significato' delle sequenze
# A livello di stringa è indifferente ma se specifichiamo l'alfabeto le sequenze non saranno uguali
# (non so come minchia si specifica visto che la classe non esiste piu)

# Classe UnknownSeq
# Il suo scopo è rappresentare una sequenza di cui conosciamo la lunghezza, ma non le lettere che la compongono.
# Possiamo usare anche un oggetto Seq però si occuperebbe molta memoria nel farlo
from Bio.Seq import UnknownSeq
print("UNKOWN SEQUENCE")
unk = UnknownSeq(20)
print("| L'output della sequenza sconosciuta (20) è     -> " + unk)
# unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna) si puo anche specificare il tipo di alfabeto
print("| Specificando l'alfabeto abbiamo 2 lettere per capire la sequenza anziche '?'")
print("|")
print("|---> N per sequenze nucleotidiche è lettera predefinita")
print("|---> X per le sequenze proteiche")
print("-------------------------------------------------------------------------------")
print("|")


# La classe SeqRecord è definita nel modulo Bio.SeqRecord e consente di associare alla sequenza
# caratteristiche di livello superiore come identificatori e caratteristiche
# Viene utilizzata nel modulo Bio.SeqIO
from Bio.SeqRecord import SeqRecord
# help(SeqRecord) è un esempio di manuale
print("SeqRecord e lettura files")
print("|")
print("|--> .seq                : Ritorna un oggetto sequenza (la sequenza in se per se")
print("|--> .id                 : L'ID primario utilizzato per identificata la sequenza (una stringa)") # e sarà il numero di accesso
print("|--> .name               : Un nome/id 'comune' per la sequenza (una stringa)") # può essere anche il numero di accesso ma attento ai doppioni in caso
print("|--> .description        : Una descrizione leggibile della sequenza")
print("|--> .letter_annotations : Lista di annotazioni per lettera utilizzando un dizionario (limitato) ")
print("|                          di informazioni aggiuntive sulle lettere nella sequenza")
print("|--> .annotations        : Ciò consente l'aggiunta di più informazioni 'non strutturate' alla sequenza")
print("|--> .features           : Un elenco di oggetti SeqFeature con informazioni più strutturate ")
print("|                          sulle caratteristiche di una sequenza come posizione dei geni su un genoma ")
print("|                          o domini su una sequenza proteica")
print("|--> .dbxrefs            : Un elenco di riferimenti incrociati al database come stringhe")
print("|")
print("| VEDERE SLIDE 29 PER UN ESEMPIO MA STA ANCHE NEL FILE MAIN")
print("|")
