#Coder: Mayur Pawashe
#File: Dependency required by five_primes.py, uorfs.py, kozak.py, nagnag.py
#Program: Shows usage to the other scripts

import os
import sys
import re

def printUsage():
    print "Usage: Choose one of the following scripts to run: five_primes.py, uorfs.py, kozak.py, nagnag.py"
    print "To run: python script_name.py <path_to_annotation_file> <directory_with_chromosome_sequence_files> <directory_to_output_files>"
    print "The first argument is required. If no other arguments are supplied, they are set to the current directory"

#Used for sorting utr's and cdss's regions by looking at the start and end points of the regions
def compareRecords(record1, record2):
    result = cmp(int(record1['start']), int(record2['start']))
    if result == 0:
        result = cmp(int(record1['end'], int(record2['end'])))
    
    return result

#Loads the entire chromosome file and returns a string of all the text with all the newlines ommitted
def getSequenceData(readFile):
    #skip first line
    readFile.readline()
    return str(readFile.read()).replace("\n", "")

def readSequence(sequenceData, startOffset, endOffset, direction):
    sequence = sequenceData[startOffset-1:endOffset].upper()
    
    if direction == '-':
        newSequence = ""
        for character in reversed(sequence):
            if character == 'A':
                newSequence += 'T'
            elif character == 'C':
                newSequence += 'G'
            elif character == 'G':
                newSequence += 'C'
            elif character == 'T':
                newSequence += 'A'
        
        return newSequence
    
    return sequence

def run(annotationFile="", chromesDirectory="", newChromosomeCallback=None, transcriptCallback=None, printUsageFunction=printUsage):
    currentChromosome = None
    currentChromsomeNotFound = None
    
    try:
        genomeFile = open(annotationFile, "r")
    except:
        printUsageFunction()
        sys.exit(2)
    sequenceFile = None
    
    records = None
    lastTranscriptID = None
    sequenceData = None
    newSequenceData = None
    
    while True:
        line = genomeFile.readline()
        
        #reached EOF
        if not line:
            break
        
        #ignore comments
        if line[:2] == "##":
            continue
        
        #Get the components of the line
        components = re.split(r'[\t ;]', line)
        
        if chromesDirectory != None:
            if currentChromsomeNotFound == components[0]:
                continue
            
            if not currentChromosome or components[0] != currentChromosome:
                if sequenceFile:
                    sequenceFile.close()
                    sequenceFile = None
                    
                if not os.path.exists(os.path.join(chromesDirectory, components[0] + ".fa")):
                    currentChromsomeNotFound = components[0]
                    print currentChromsomeNotFound + str(".fa does not exist, ignoring..")
                    continue
                currentChromosome = components[0]
                print "Loading " + currentChromosome + ".."
                sequenceFile = open(os.path.join(chromesDirectory, currentChromosome + ".fa"), "r")
                newSequenceData = getSequenceData(sequenceFile)
                if not sequenceData:
                    sequenceData = newSequenceData
                sequenceFile.close()
                
                if newChromosomeCallback:
                    newChromosomeCallback(currentChromosome)
                
                print "Writing Data.."
            
        transcriptID = components[12].strip('"')
        
        if lastTranscriptID != transcriptID:
            if records:
                utrs = []
                cdss = []
                stopCodon = None
                direction = records[0]['direction']
                
                for record in records:
                    if record['type'] == 'UTR':
                        utrs.append(record)
                    elif record['type'] == 'CDS':
                        cdss.append(record)
                    elif record['type'] == 'stop_codon':
                        stopCodon = record
                
                utrs.sort(compareRecords)
                cdss.sort(compareRecords)
                
                fivePrimeUTRs = []
                threePrimeUTRs = []
                
                for utr in utrs:
                    if utr['start'] <= cdss[0]['start']:
                        if direction == '+':
                            fivePrimeUTRs.append(utr)
                        else:
                            threePrimeUTRs.append(utr)
                    elif utr['start'] >= cdss[-1]['start']:
                        if direction == '+':
                            threePrimeUTRs.append(utr)
                        else:
                            fivePrimeUTRs.append(utr)
                
                if direction == '-':
                    #Read in reverse order
                    fivePrimeUTRs.reverse()
                    threePrimeUTRs.reverse()
                    cdss.reverse()
                    
                fivePrimeContent = ""
                if sequenceData:
                    for fivePrime in fivePrimeUTRs:
                        fivePrimeContent += readSequence(sequenceData, fivePrime['start'], fivePrime['end'], direction)
                
                if transcriptCallback:
                    transcriptCallback(records, sequenceData, lastTranscriptID, fivePrimeUTRs, fivePrimeContent, cdss, stopCodon, direction)
    			
                
            lastTranscriptID = transcriptID
            records = []
            if newSequenceData:
                sequenceData = newSequenceData
                newSequenceData = None
        
        records.append({'line' : line, 'chromosome' : components[0], 'type' : components[2], 'start' : int(components[3]), 'end' : int(components[4]), 'direction' : components[6], 'geneID' : components[9].strip('"'), 'geneName' : components[21].strip('"'), 'geneStatus' : components[18].strip('"'), 'geneType' : components[15].strip('"'), 'transcriptStatus' : components[27].strip('"'), 'transcriptName' : components[30].strip('"'), 'transcriptType' : components[24].strip('"'), 'levelNumber' : components[33]})
    
    genomeFile.close()

def getArguments(printUsageFunction=printUsage):
    #First argument is required
    if len(sys.argv) < 2:
        printUsageFunction()
        sys.exit(1)
        
    arguments = []
    for argument in sys.argv[1:]:
        arguments.append(argument)
    
    #3 required arguments
    while len(arguments) < 3:
        arguments.append("")
    
    return arguments

if __name__ == '__main__':
    printUsage()
