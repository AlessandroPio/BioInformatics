import os, glob
import time
from Bio import AlignIO


def menu():
    print("|---------------------------------------------------")
    print("|")
    while True:
        directory_file = str(input("| Enter the working directory                      -> "))
        try:
            print("| Checking the directory..")
            os.chdir(directory_file)
            print("| Directory controlled!")
            print("|");break
        except Exception:
            print("| ENTER A VALID DIRECTORY!")
    while True:
        print("| IMPORTANT!!")
        clustal_directory = str(input("| Enter the Clustal directory                      -> "))
        try:
            print("| Checking the directory..")
            os.chdir(clustal_directory)
            print("| Directory controlled!");print("|")
            print("| Checking the executable..")
            if "clustalo.exe" in glob.glob("clustalo.exe"):
                print("| Executable checked!")
                os.chdir(directory_file); break
            else:
                print("| EXECUTABLE MISSING!")
        except Exception:
            print("| ENTER A VALID DIRECTORY!")

    while True:
        print("|"); cont=0
        input_file = str(input ("| Insert the .fasta format file                    -> "))
        if input_file in glob.glob("*.fasta"):
            cont += 1
        else:
            print("| FASTA FILE MISSING!")
        print("| Insert the .fasta output file")
        output_file = str(input("| (leave blank if you don't want to specify it)    -> "))
        if not output_file:
            output_file = "SequencesAligned.fasta"
            if(output_file in glob.glob("*.fasta")):
                os.remove(output_file)
            print("| Base file 'SequenceAligned.fasta' ")
            print("|")
        if(input_file in glob.glob("*.fasta")):
            cont += 1
            if (cont == 2): break
        else:
            print("| ENTER VALID .FASTA INPUT FILE!")
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
    print("| Alignment started at  -> " + str(init_time.tm_hour) + ":" + str(init_time.tm_min) + ":" + str(init_time.tm_sec))
    print("|")
    print("| Alignment in progress...")
    os.system(clustal_directory + "\clustalo -i " + input_file + " -o " + output_file)
    print("| Alignment concluded!")
    print("|")
    end_time = time.localtime()
    print("| Alignment concluded at -> " + str(end_time.tm_hour) + ":" + str(end_time.tm_min) + ":" + str(end_time.tm_sec))
    print("|")


fields = menu()
align(fields[0], fields[1], fields[2])
getMutation(fields[1])