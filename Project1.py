import os, glob
import time
import numpy as np
from Bio import AlignIO
from prettytable import PrettyTable

class bcolors:
    SEQUENCE = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    OK = '\033[92m'

def menu():
    print("|---------------------------------------------------")
    print("|")
    while True:
        directory_file = str(input("| Enter the working directory                      -> "))
        try:
            print("| Checking the directory..")
            os.chdir(directory_file)
            print("| " + bcolors.OK + "Working directory controlled!" + bcolors.ENDC)
            print("|");
            break
        except Exception:
            print("| " + bcolors.FAIL + "ENTER A VALID DIRECTORY!" + bcolors.ENDC);print("|")
    while True:
        print("| IMPORTANT!!")
        clustal_directory = str(input("| Enter the Clustal directory                      -> "))
        try:
            print("| Checking the directory..")
            os.chdir(clustal_directory)
            print("| " + bcolors.OK + "Directory controlled!" + bcolors.ENDC)
            print("|")
            print("| Checking the executable..")
            if "clustalo.exe" in glob.glob("clustalo.exe"):
                print("| " + bcolors.OK + "Executable checked!" + bcolors.ENDC)
                os.chdir(directory_file);
                break
            else:
                print("| " + bcolors.FAIL + "EXECUTABLE MISSING!" + bcolors.ENDC);print("|")
        except Exception:
            print("| " + bcolors.FAIL + "ENTER A VALID DIRECTORY!" + bcolors.ENDC);print("|")

    while True:
        print("|");
        cont = 0
        input_file = str(input("| Insert the .fasta format file                    -> "))
        if input_file in glob.glob("*.fasta"):
            cont += 1
        else:
            print("| " + bcolors.FAIL + "FASTA FILE MISSING!" + bcolors.ENDC);print("|")
        print("| Insert the .fasta output file")
        output_file = str(input("| (leave blank if you don't want to specify it)    -> "))
        if not output_file:
            output_file = "SequencesAligned.fasta"
            if (output_file in glob.glob("*.fasta")):
                os.remove(output_file)
            print("| " + bcolors.SEQUENCE + "Base file 'SequenceAligned.fasta' " + bcolors.ENDC)
            print("|")
        if (input_file in glob.glob("*.fasta")):
            cont += 1
            if (cont == 2): break
        else:
            print("| " + bcolors.FAIL + "ENTER VALID .FASTA INPUT FILE!" + bcolors.ENDC);print("|")
    print("|---------------------------------------------------")
    return (input_file, output_file, clustal_directory, directory_file)

def getMutations(output_file):
    def getIds(output_file):
        ids = []
        for records in output_file:
            ids.append(records.id)

        return ids

    print("|"); print("| Calculating Muations...")
    fas = AlignIO.read(output_file, 'fasta')
    seq_record = np.array(fas)
    res = seq_record.transpose()

    currSeq = []
    currMut = []; i = 0
    while i < len(res):
        j = 0
        cond = 0
        while j < len(res[i]) and not cond:
            z = 0
            while z < len(res[i]):
                if j != z:
                    if res[i][j] != res[i][z]:
                        currSeq.append(res[i])
                        currMut.append(i)
                        cond = 1
                        break
                z += 1
            j += 1
        i += 1

    return(currMut , list(currSeq), getIds(fas))

def align(input_file, output_file, clustal_directory):
    print("|")
    init_time = time.localtime()
    print("| Alignment started at  -> " + bcolors.SEQUENCE + str(init_time.tm_hour) + ":" + str(init_time.tm_min) + ":" + str(
        init_time.tm_sec) + bcolors.ENDC)
    print("|")
    print("| Alignment in progress...")
    os.system(clustal_directory + "\clustalo -i " + input_file + " -o " + output_file)
    print("| " + bcolors.OK + "Alignment concluded!" + bcolors.ENDC)
    print("|")
    end_time = time.localtime()
    print("| Alignment concluded at -> " + bcolors.SEQUENCE + str(end_time.tm_hour) + ":" + str(end_time.tm_min) + ":" + str(
        end_time.tm_sec) + bcolors.ENDC)
    print("|")

def getTable(pos, data, ids, directory_file):
    def generateCSV(table):
        file = open(directory_file + "/TableMutations.csv", "w")
        file.write(table.get_csv_string())
        file.close()

    x = PrettyTable()
    x.add_column("ID", ids)
    i=0

    while i < len(data):
        x.add_column(str(pos[i]), data[i])
        i += 1

    print("| Process ended sucsesfully!"); print("|")
    print("| Creating a csv file...")
    generateCSV(x); print("| " + bcolors.OK + "TableMutations.csv created sucsessfully inside working directory!" + bcolors.ENDC)
    print("|"); print("| " + bcolors.SEQUENCE + "Printing mutations table..." + bcolors.ENDC)
    print("|"); print(x)


fields = menu()
align(fields[0], fields[1], fields[2])
#directory_file = "/Users/alessandro/Desktop/Workspace/BioInformatics/sequences/Project1/"
#fields = "/Users/alessandro/Desktop/Workspace/BioInformatics/sequences/Project1/Proj1SequencesAligned.fasta"
data, column, ids = getMutations(fields[1])
getTable(data, column, ids, fields[3])