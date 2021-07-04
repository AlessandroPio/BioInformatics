import os, glob
from os.path import exists
import time

#os.chdir("sequences/Project1/")

#if(platform.system() != 'Windows'): #type
#    os.system("cat Proj1VarIndia.fasta Proj1VarItalia.fasta | awk '/^>/{f=!d[$1];d[$1]=1}f' > sequences.fasta")
from Bio import AlignIO


def menu():
    print("|---------------------------------------------------")
    print("|")
    while True:
        directory_file = str(input("| Inserisci la directory di lavoro                -> "))
        try:
            print("| Controllo directory..")
            os.chdir(directory_file)
            print("| Directory controllata!")
            print("|");break
        except Exception:
            print("| IMMETTI UNA DIRECTORY VALIDA!")
    while True:
        print("| IMPORTANTE!!")
        clustal_directory = str(input("| Inserisci la directory di Clustal                -> "))
        try:
            print("| Controllo directory..")
            os.chdir(clustal_directory)
            print("| Directory controllata!")
            print("| Controllo eseguibile..")
            if "clustalo.exe" in glob.glob("clustalo.exe"):
                print("| Eseguibile controllato!")
                os.chdir(directory_file); break
            else:
                print("| Eseguibile Mancante!")
        except Exception:
            print("| IMMETTI UNA DIRECTORY VALIDA!")

    while True:
        print("|")
        input_file = str(input ("| Inserisci il file formato .fasta                  ->"))
        print("| Inserisci il file output .fasta")
        output_file = str(input("| (lasciare bianco se non si vuole specificarlo)    ->"))
        if not output_file:
            output_file = "SequencesAligned.fasta"
            print("| File di base 'SequenceAligned.fasta' ")
            print("|")
        if(input_file in glob.glob("*.fasta")):
            break
        else:
            print("| IMMETTI FILE INPUT .FASTA VALIDO!")
    print("|---------------------------------------------------")

    return(input_file, output_file, clustal_directory)


def getMutation(output_file):
    seq_records = AlignIO.read(output_file, 'fasta')
    print("|")
    print("| DETAILS")
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
            print("| Mutation Detected!   (" + str(j) + ")")
            print("| Sequences Id         ->", seq_records[j].id)
            print("| Lenght of sequence   ->", len(seq_records[j]))
            print("| Number of mutations  ->", len(location))
            print("|")
        j += 1
    print("|---------------------------------------------------")

def align(input_file, output_file, clustal_directory):

    print("|")
    init_time = time.localtime()
    print("| Allineamento iniziato alle  -> " + str(init_time.tm_hour) + ":" + str(init_time.tm_min) + ":" + str(
        init_time.tm_sec))
    print("|")
    print("| Sto allineando...")
    os.system(clustal_directory + "\clustalo -i " + input_file + " -o " + output_file)
    print("| Allineamento concluso!")
    print("|")
    end_time = time.localtime()
    print("| Allineamento concluso alle  -> " + str(end_time.tm_hour) + ":" + str(end_time.tm_min) + ":" + str(
        end_time.tm_sec))
    print("|")


fields = menu()
align(fields[0], fields[1], fields[2])
getMutation(fields[1])