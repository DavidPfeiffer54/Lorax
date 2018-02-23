import glob, os
import re
import json 
from collections import Counter

os.chdir("loraxDataTest") #changes the directory to the one with all the data files 
TAXLABELS = re.compile("\[(\d+)\] '(.*)'") # [44] 'Unionidae_Unioninae_Cuneopsis_pisciculus'
Matrix = re.compile("(\d+)\s+((\d+( |,))+)") #  29   3 4 41 43 45 46 47 48,
cladesInput = re.compile("(.*) = (.*)") # Blue_Clams = Unionidae_Rectidentinae_Contradens_contradens, Unionidae_Rectidentinae_Solenaia_khwaenoiensis
nameToIndex = {}  # {speciesName:index}
indexToName = {}  # {index:speciesName}
matrixValues = [] # [[boostrap, indexlist],[boostrap, indexlist]]
cladesOfInterest = {} # {cladeName: speciesInCladeList}
cladesOfInterestIndexes = {} # {cladeName : splitIndexOfSpeciesInCladeList}

def populateNameIndexDicts(lines):
    '''
    populates dictionaries for converting between a speciesName and its Index within a split
    Parameters
    ----------
    lines : list
        contents of split file. each item is 1 line of the file
    Returns
    -------
    Dict
        speciesName -> index
    Dict
        index -> speciesName

    '''
    for l in lines:
        l=l.strip('\n')
        match = re.search(TAXLABELS, l)
        if match:
            index = match.group(1)
            speciesName = match.group(2)

            nameToIndex[speciesName] = index 
            indexToName[index] = speciesName
    return nameToIndex, indexToName


def populateMatrix(lines):
    '''
    finds all of the bootstrap-clade associations in the split files
    each item will is [bootstrapValue, ListOfClades]
    Parameters
    ----------
    lines : list
        contents of split file. each item is 1 line of the file
    Returns
    -------
    list of lists 
        A list of the form [['bootstrap',[index1, index2]],'bootstrap',[index1]]]  
    '''

    for l in lines:
        l=l.strip('\n')
        match = re.search(Matrix, l)
        if match:
            bootstrapSupport = match.group(1)
            speciesList = match.group(2)
            speciesList=speciesList.replace(',','')
            speciesList = [x.strip() for x in speciesList.split(' ')]
            matrixValues.append([bootstrapSupport,speciesList])
    return matrixValues

def populateCladesOfInterest(lines):
    '''
    reads clades of interest file and popluates a dictionary 
    Parameters
    ----------
    lines : list
        contents of clades file. each item is 1 line of the file
    Returns
    -------
    dict 
        A dictionary of the form {'cladeName':[speciesName, species2Name]}   
    '''
    for l in lines:
        l=l.strip('\n') #removes newlines.
        match = re.search(cladesInput, l) #look for lines that match the pattern: "(\d+)\s+((\d+( |,))+)"
        if match:
            cladeDescriptor = match.group(1)
            speciesInClade = match.group(2)
            speciesInClade = [x.strip() for x in speciesInClade.split(',')] #turns a string "1,2,3,4" into an python list [1,2,3,4]
            cladesOfInterest[cladeDescriptor] = speciesInClade
    return cladesOfInterest


def getCladeOfInterstIndexedToNex(cladesOfInterest, nameToIndex):
    ''' 
    loops through all of the clades defined in the cladesOfInterest file.
    replaces the species name, to the index for that species in this current split file
    if a species in the Clade is not found in the split, it is omitted

    Parameters
    ----------
    cladesOfInterest : dict
        A dictionary of the form {'cladeName':[speciesName, species2Name]}
    nameToIndex: dict
        a dictionary of the form {speciesName:1, species2Name:2}
    Returns
    -------
    dict 
        A dictionary of the form {'cladeName':[1, 2]}
    '''
    for cladeName, speciesList in cladesOfInterest.items(): #loop through all clade
        indexes = []
        for v in speciesList: 
            if v in nameToIndex: #check if the species named in the clade exists in this split
                indexes.append(nameToIndex[v]) 
        cladesOfInterestIndexes[cladeName] = indexes #update the new dictionary 
    return cladesOfInterestIndexes

def compareClades(l1, l2):
    '''
    compares 2 lists and returns wether or not they contain all of the same items

    Parameters
    ----------
    l1 : list
        A dictionary of the form [{'host': <hostname>,
        'time': <datetime>, 'vNode': <int>, 'oldPnode': <int>,
        'newPnode': <int>}, ...]
    l2: list
        The vNode to test crash for.

    Returns
    -------
    bool 
        wether the list

    '''
    return Counter(l1) == Counter(l2)

#set up the clades of interst file
with open("CladesOfInterest.txt", "r") as cladesFile: #open the file for reading
    lines = cladesFile.readlines() #reads in the entire file into a list - 1 line per item
    populateCladesOfInterest(lines)


#open up output file, wite out its headers
outputFile  = open("outputFile.txt", "w")  
outputFile.write("LOCI")
for cladeName, species in cladesOfInterest.items():
    outputFile.write(",\t"+str(cladeName))
outputFile.write("\n")


for fileName in glob.glob("*.nex"): #this will effect every file that has '.txt' - so only put the needed files in this directory
    #reset these values for each split file
    nameToIndex = {}
    indexToName = {}
    matrixValues = []
    cladesOfInterestIndexes = {}

    with open(fileName, "r") as nexFile: #open the file for reading
        lines = nexFile.readlines() #reads in the entire file
        nameToIndex, indexToName = populateNameIndexDicts(lines) #populate the index dictionaries
        matrixValues = populateMatrix(lines) #populate the bootstraps

    cladesOfInterestIndexes = getCladeOfInterstIndexedToNex(cladesOfInterest, nameToIndex) #convert clade filenames to indexes
    
    #print json.dumps(nameToIndex, indent=1)
    #print json.dumps(indexToName, indent=1)
    #print json.dumps(matrixValues, indent=1)
    #print json.dumps(cladesOfInterest, indent=1)
    #print json.dumps(cladesOfInterestIndexes, indent=1)

    outputFile.write(fileName)
    for clade, cladeIndexes in cladesOfInterestIndexes.items():
        confidenceLevel = "N/A" #defualt value, if we do not find the clade in the split
        for mx in matrixValues:
            bootstrapValue = mx[0]
            speciesIndexList = mx[1]
            if compareClades(speciesIndexList, cladeIndexes):
                confidenceLevel = str(bootstrapValue)
        outputFile.write(",\t"+str(confidenceLevel))
    outputFile.write("\n")







