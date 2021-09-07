from Bio import SeqIO, Phylo
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os, glob, time
from terminaltables import AsciiTable

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
    return (input_file, output_file, clustal_directory)

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

def GetInfo(informations):
    data = []
    data.append(["ID", "SEQUENCE'S LENGTH", "DESCRIPTION"])
    for info in SeqIO.parse(informations, 'fasta'):
        data.append([info.id, len(info.seq), info.description])

    return data

def phylipTrascription(input_file):
    records = SeqIO.parse(input_file, "fasta")
    SeqIO.write(records, input_file.split(".")[0] + ".phylip", "phylip")

def NJTree(name, distanceMatrix, distanceConstructor):
    print("| NJTree Construction...")
    NJTree = distanceConstructor.nj(distanceMatrix)
    Phylo.draw(NJTree)
    Phylo.write(NJTree, name.split(".")[0] + "NJTree.xml", "phyloxml")
    print("| " + bcolors.OK + "Created " + name.split(".")[0] + "NJTree.xml on your work directory" + bcolors.ENDC)
    print("|")

def UPGMAtree(name, distanceMatrix, distanceConstructor):
    print("| UPGMAtree Construction...")
    UPGMATree = distanceConstructor.upgma(distanceMatrix)
    Phylo.draw(UPGMATree)
    Phylo.write(UPGMATree, name.split(".")[0] + "UPGMATree.xml", "phyloxml")
    print("| " + bcolors.OK + "Created" + name.split(".")[0] + "UPGMATree.xml on your work directory" + bcolors.ENDC)
    print("|")


fields = menu()
align(fields[0], fields[1], fields[2])
phylipTrascription(fields[1])

distanceConstructor = DistanceTreeConstructor()
distanceCalculator = DistanceCalculator('blosum62')
alignments = AlignIO.read(fields[1].split(".")[0] + ".phylip", "phylip")
table = AsciiTable(GetInfo(fields[0]))
print("| " + bcolors.SEQUENCE + "DETAILS" + bcolors.ENDC);print("|")
print(table.table);print("|")

distanceMatrix = distanceCalculator.get_distance(alignments)
NJTree(fields[1], distanceMatrix, distanceConstructor)
UPGMAtree(fields[1], distanceMatrix, distanceConstructor)

os.remove(fields[1].split(".")[0] + ".phylip")
print("| " + bcolors.OK + "Process finished!" + bcolors.ENDC)
print("|---------------------------------------------------")