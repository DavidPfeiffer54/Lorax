import glob, os
import re
import json 
from collections import Counter

os.chdir("loraxDataTest") #changes the directory to the one with all the data files - call the directory 'data'
TAXLABELS = re.compile("\[(\d+)\] '(.*)'") # [44] 'Unionidae_Unioninae_Cuneopsis_pisciculus'
Matrix = re.compile("(\d+)\s+((\d+( |,))+)") #  29   3 4 41 43 45 46 47 48,
cladesInput = re.compile("(.*) = (.*)") # Blue_Clams = Unionidae_Rectidentinae_Contradens_contradens, Unionidae_Rectidentinae_Solenaia_khwaenoiensis
nameToIndex = {}
indexToName = {}
matrixValues = []
cladesOfInterest = {}
cladesOfInterestIndexes = {}

def populateNameIndexDicts(lines):
    for l in lines:
        l=l.strip('\n')
        match = re.search(TAXLABELS, l)
        if match:
            index = match.group(1)
            name = match.group(2)
            
            nameToIndex[name] = index
            indexToName[index] = name


def populateMatrix(lines):
    for l in lines:
        l=l.strip('\n')
        match = re.search(Matrix, l)
        if match:
            val = match.group(1)
            lst = match.group(2)
            lst=lst.replace(',','')
            lst = [x.strip() for x in lst.split(' ')]
            matrixValues.append([val,lst])

def populateCladesOfInterest(lines):
    for l in lines:
        l=l.strip('\n')
        match = re.search(cladesInput, l)
        if match:
            cladeDescriptor = match.group(1)
            namesInClade = match.group(2)
            namesInClade = [x.strip() for x in namesInClade.split(',')]

            cladesOfInterest[cladeDescriptor] = namesInClade



def getCladeInNex():
    for name, vals in cladesOfInterest.items():
        indexes = []
        for v in vals:
            if v in nameToIndex:
                indexes.append(nameToIndex[v])
        cladesOfInterestIndexes[name] = indexes 


def compareClades(s, t):
    return Counter(s) == Counter(t)

with open("CladesOfInterest.txt", "r") as cladesFile: #open the file for reading
    lines = cladesFile.readlines() #reads in the entire file
    populateCladesOfInterest(lines)

outputFile  = open("outputFile.txt", "w")  # opens up the new file to write
outputFile.write("LOCI")
for cladeName, species in cladesOfInterest.items():
    outputFile.write(",\t"+str(cladeName))
outputFile.write("\n")

for fileName in glob.glob("*.nex"): #this will effect every file that has '.txt' - so only put the needed files in this directory
    print " ------------ "
    nameToIndex = {}
    indexToName = {}
    matrixValues = []
    cladesOfInterestIndexes = {}

    with open(fileName, "r") as nexFile: #open the file for reading
        lines = nexFile.readlines() #reads in the entire file
        populateNameIndexDicts(lines)
        populateMatrix(lines)



    getCladeInNex()




    #populateNameIndexDicts(lines)
    #populateMatrix(lines)
    #populateCladesOfInterest(clades)
    
    #print json.dumps(nameToIndex, indent=1)
    #print json.dumps(indexToName, indent=1)
    #print json.dumps(matrixValues, indent=1)
    #print json.dumps(cladesOfInterest, indent=1)
    #print json.dumps(cladesOfInterestIndexes, indent=1)

    outputFile.write(fileName)
    for clade, cladeIndexes in cladesOfInterestIndexes.items():
        confidenceLevel = "N/A"
        for mx in matrixValues:
            lst = mx[1]
            if compareClades(lst, cladeIndexes):
                confidenceLevel = str(mx[0])
                print("!!!!!!!!!!!!!!!!!!!FOUND ONE " +str(clade)+" "+ str(mx[0])+ "    " +str(lst))
        outputFile.write(",\t"+str(confidenceLevel))
    outputFile.write("\n")







