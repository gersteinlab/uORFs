#Coder: Russell Ault
#Program: Writes out uorf files
#Dependencies: gencode.py

import os
import sys
import gencode

global StartCodon
global uorfAnnotationFile
global uorfSequenceFile
global uorfUnalignedAnnotationFile
global uorfUnalignedSequenceFile

def writeUORFOutput(outputFile, transcript, fivePrimeUTRs, cdssAndStopCodon, startIndex, endIndex, direction):
    lengthAccumulator = 0 #actual length read in so far
    
    records = fivePrimeUTRs[:]
    records.extend(cdssAndStopCodon)
    
    foundFirstStart = False
    
    for record in records:
        startPosition = None
        endPosition = None
        lengthAccumulator += record['end'] - (record['start'] - 1)
        
        if not foundFirstStart and lengthAccumulator > startIndex:
            if direction == '+':
                startPosition = record['end'] - (lengthAccumulator - startIndex) + 1
            elif direction == '-':
                startPosition = record['start'] + (lengthAccumulator - startIndex) - 1
                
            foundFirstStart = True
        
        if foundFirstStart and lengthAccumulator >= endIndex:
            if direction == '+':
                endPosition = record['end'] - (lengthAccumulator - endIndex)
            elif direction == '-':
                endPosition = record['start'] + (lengthAccumulator - endIndex)
        
        
        if foundFirstStart:
            shouldBreak = False
            if not startPosition and not endPosition:
                startPosition = record['start']
                endPosition = record['end']
            elif startPosition and endPosition:
                shouldBreak = True
                if direction == '-':
                    startPosition, endPosition = endPosition, startPosition
            elif startPosition:
                if direction == '+':
                    endPosition = record['end']
                elif direction == '-':
                    endPosition = startPosition
                    startPosition = record['start']
            elif endPosition:
                shouldBreak = True
                if direction == '+':
                    startPosition = record['start']
                elif direction == '-':
                    startPosition = endPosition
                    endPosition = record['end']
                            
            outputFile.write(record['chromosome'] + "  mayur uORF  " + str(startPosition) + " " + str(endPosition) + "   .   " + direction + "   .   " + "gene_id \"" + record['geneID'] + "\"; " + "transcript_id \"" + transcript + "." + record['type'] + "\"; gene_type \"uORF\"; gene_status \"" + record['geneStatus'] + "\"; gene_name \"" + record['geneName'] + "\"; transcript_type \"uORF\"; transcript_status \"" + record['transcriptStatus'] + "\"; transcript_name \"" + record['transcriptName'] + "\"; level " + record['levelNumber'] + ";" + "\n")
            
            if shouldBreak:
                break

def findUORF(outputFile, outputSequenceFile, sequenceData, divisibleByThree, transcript, fivePrimeSequence, fivePrimeUTRs, cdss, stopCodon, direction): 
    cdsSequence = ""
    cdssAndStopCodon = cdss[:]
    
    for cds in cdssAndStopCodon:
        cdsSequence += gencode.readSequence(sequenceData, cds['start'], cds['end'], direction)
    
    stopCodonSequence = ""
    if stopCodon:
        #stop codon is always last
        cdssAndStopCodon.append(stopCodon)
        stopCodonSequence += gencode.readSequence(sequenceData, stopCodon['start'], stopCodon['end'], direction)
    
    sequence = fivePrimeSequence + cdsSequence + stopCodonSequence
    
    atgOffsetToSearch = 0
    atgSequenceLength = len(fivePrimeSequence)
    sequenceLength = len(sequence)
    
    transcriptNumber = 1
    
    while atgOffsetToSearch < atgSequenceLength - 2:
        atgIndex = fivePrimeSequence.find(StartCodon, atgOffsetToSearch)
        if atgIndex == -1:
            break
        
        atgOffsetToSearch = atgIndex + 1
        
        if atgIndex < atgSequenceLength - 2:
            endIndex = atgIndex + 3
            while endIndex < sequenceLength - 2:
                threeCharacterSequence = sequence[endIndex:endIndex+3]
                if threeCharacterSequence == 'TAG' or threeCharacterSequence == 'TAA' or threeCharacterSequence == 'TGA':
                    if not divisibleByThree or (endIndex - atgIndex) % 3 == 0:
                        uorfSequence = sequence[atgIndex:endIndex+3]
                        howFarReached = ""
                        if endIndex + 3 <= len(fivePrimeSequence):
                            howFarReached = "UTRonly"
                        elif endIndex + 3 <= len(fivePrimeSequence) + len(cdsSequence):
                            howFarReached = "CDSpartial"
                        else:
                            howFarReached = "CDSfull"
                            
                        record = fivePrimeUTRs[0]
                        outputSequenceFile.write(record["chromosome"] + "|" + record["geneID"] + "|" + transcript + ".uORF_" + StartCodon + "." + str(transcriptNumber) + "|" + direction + "|" + record["geneName"] + "|" + howFarReached + "\n")
                        outputSequenceFile.write(uorfSequence + "\n")
                        writeUORFOutput(outputFile, transcript + ".uORF_" + StartCodon + "." + str(transcriptNumber), fivePrimeUTRs, cdssAndStopCodon, atgIndex, endIndex+3, direction)
                        transcriptNumber += 1
                        break
                
                endIndex += 1

def transcriptCallback(records, sequenceData, lastTranscriptID, fivePrimeUTRs, fivePrimeContent, cdss, stopCodon, direction):
    if fivePrimeContent != "":
        findUORF(uorfAnnotationFile, uorfSequenceFile, sequenceData, True, lastTranscriptID, fivePrimeContent, fivePrimeUTRs, cdss, stopCodon, direction)
        findUORF(uorfUnalignedAnnotationFile, uorfUnalignedSequenceFile, sequenceData, False, lastTranscriptID, fivePrimeContent, fivePrimeUTRs, cdss, stopCodon, direction)

if __name__ == '__main__':
    arguments = gencode.getArguments()
    if len(sys.argv) == 5:
        if sys.argv[4][0] == '-':
            StartCodon = sys.argv[4][1:]
    else:
        StartCodon = 'ATG'
    
    uorfAnnotationFile = open(os.path.join(arguments[2], "uORFs_" + StartCodon + ".gtf"), "w") 
    uorfSequenceFile = open(os.path.join(arguments[2], "uORFs_" + StartCodon + ".fa"), "w")
    
    uorfUnalignedAnnotationFile = open(os.path.join(arguments[2], "uORFs_" + StartCodon + "_unaligned.gtf"), "w")
    uorfUnalignedSequenceFile = open(os.path.join(arguments[2], "uORFs_" + StartCodon + "_unaligned.fa"), "w")
    
    uorfSequenceFile.write("##chr|gene ID|transcript ID|strand|gene name|UTRonly/CDSpartial/CDSfull\n")
    uorfUnalignedSequenceFile.write("##chr|gene ID|transcript ID|strand|gene name|UTRonly/CDSpartial/CDSfull\n")
    
    gencode.run(arguments[0], arguments[1], None, transcriptCallback)
    
    uorfUnalignedAnnotationFile.close()
    uorfUnalignedSequenceFile.close()
    uorfSequenceFile.close()
    uorfAnnotationFile.close()
