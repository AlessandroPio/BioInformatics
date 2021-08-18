import os, glob
import time
from itertools import permutations, combinations

import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class bcolors:
    SEQUENCE = '\033[93m'
    MUTATION = '\033[91m'
    ENDC = '\033[0m'

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

def getInfo(ID, records):
    data = []; i = 0
    data.append(["ID","SEQUENCES"])
    for id in ID:
        data.append([id, records[i].seq])
        i += 1

    return data


def getMutations(output_file):
    def listToString(s):
        str1 = ""
        for ele in s:
            str1 += ele
        return str1

    fas = AlignIO.read(output_file, 'fasta')
    seq_records = np.array(fas)

    #seq_records=[SeqRecord(Seq("AgA"),id="YP_025292.1"), SeqRecord(Seq("CTA"),id="YP_025122.2"), SeqRecord(Seq("GAA"),id="YP_025122.3")]
    seq_record =np.array(seq_records)
    res = seq_record.transpose()
    i=0; currSeq = []
    while i < len(res):
        j=0
        cond = 0
        while j < len(res[i]) and not cond:
            z=0
            while z < len(res[i]):
                if(j!=z):
                    if(res[i][j] != res[i][z]):
                        currSeq.append(res[i])
                        cond = 1
                        break
                z+=1
            j+=1
        i += 1
    original= np.transpose(currSeq)

    i=0
    while i<len(original):
        print(bcolors.SEQUENCE + fas[i].id + bcolors.ENDC + " | " + bcolors.MUTATION + listToString(original[i]) + bcolors.ENDC)
        i += 1

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


#fields = menu()
#align(fields[0], fields[1], fields[2])
fields=[0,0]
fields[1]= 'Proj1SequencesAligned.fasta'
getMutations('/Users/alessandro/Desktop/Universita/3 Anno/BioInformatica/BioInformatics/sequences/Project1/' + fields[1])