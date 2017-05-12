#Coder: Mayur Pawashe
#Program: Writes out kozak files
#Dependencies: gencode.py

import os
import gencode

global fivePrimesKozakFile
global cdsKozakFile
global cdsKozak2File

def writeKozakOutput2(outputFile, transcript, cdssAndStopCodon, genomeRecord, startIndex, endIndex, direction):
    lengthAccumulator = 0 #actual length read in so far
    
    records = cdssAndStopCodon[:]
    records.append(genomeRecord)

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

            outputFile.write(record['chromosome'] + "  mayur kozak_cds  " + str(startPosition) + " " + str(endPosition) + "   .   " + direction + "   .   " + "gene_id \"" + record['geneID'] + "\"; " + "transcript_id \"" + transcript + "." + record['type'] + "\"; gene_type \"kozak_cds\"; gene_status \"" + record['geneStatus'] + "\"; gene_name \"" + record['geneName'] + "\"; transcript_type \"kozak_cds\"; transcript_status \"" + record['transcriptStatus'] + "\"; transcript_name \"" + record['transcriptName'] + "\"; level " + record['levelNumber'] + ";" + "\n")

            if shouldBreak:
                break

def findKozak2(outputFile, sequenceData, transcript, cdss, stopCodon, direction): 
    cdsSequence = ""
    cdssAndStopCodon = cdss[:]

    for cds in cdssAndStopCodon:
        cdsSequence += gencode.readSequence(sequenceData, cds['start'], cds['end'], direction)

    stopCodonSequence = ""
    if stopCodon:
        #stop codon is always last
        cdssAndStopCodon.append(stopCodon)
        stopCodonSequence += gencode.readSequence(sequenceData, stopCodon['start'], stopCodon['end'], direction)

    afterCDSAndStopCodonSequence = ""
    genomeRecord = {'chromosome' : cdssAndStopCodon[0]['chromosome'], 'type' : 'genome', 'transcriptName' : cdssAndStopCodon[0]['transcriptName'], 'geneID' : cdssAndStopCodon[0]['geneID'], 'geneName' : cdssAndStopCodon[0]['geneName'], 'transcriptStatus' : cdssAndStopCodon[0]['transcriptStatus'], 'geneStatus' : cdssAndStopCodon[0]['geneStatus'], 'levelNumber' : ''}
    
    sequence = cdsSequence + stopCodonSequence + afterCDSAndStopCodonSequence

    atgSequenceLength = len(cdsSequence)
    atgOffsetToSearch = 0
    sequenceLength = len(sequence)

    transcriptNumber = 1
    
    genomeLength = 300

    while atgOffsetToSearch < atgSequenceLength - 2:
        atgIndex = cdsSequence.find('ATG', atgOffsetToSearch)
        if atgIndex == -1:
            break

        atgOffsetToSearch = atgIndex + 1

        if atgIndex < atgSequenceLength - 2:
            endIndex = atgIndex + 3
            
            while endIndex < sequenceLength - 2:
                threeCharacterSequence = sequence[endIndex:endIndex+3]
                if threeCharacterSequence == 'TAG' or threeCharacterSequence == 'TAA' or threeCharacterSequence == 'TGA':
                    kozakSequence = sequence[atgIndex:endIndex+3]
                    writeKozakOutput2(outputFile, transcript + ".kozak_cds." + str(transcriptNumber), cdssAndStopCodon, genomeRecord, atgIndex, endIndex+3, direction)
                    
                    transcriptNumber += 1
                    found = True
                    break

                endIndex += 3
                
                if endIndex >= sequenceLength - 2 and sequenceLength < len(sequenceData):
                    if direction == '+':
                        afterCDSAndStopCodonSequence = gencode.readSequence(sequenceData, cdssAndStopCodon[-1]['end'] + 1, min(cdssAndStopCodon[-1]['end'] + genomeLength, len(sequenceData)), direction)
                        genomeRecord['start'] = cdssAndStopCodon[-1]['end'] + 1
                        genomeRecord['end'] = len(sequenceData)
                    elif direction == '-':
                        afterCDSAndStopCodonSequence = gencode.readSequence(sequenceData, max(cdssAndStopCodon[-1]['start'] - genomeLength, 1), cdssAndStopCodon[-1]['start'] - 1, direction)
                        genomeRecord['start'] = 1
                        genomeRecord['end'] = cdssAndStopCodon[-1]['start'] - 1
                    
                    sequence = cdsSequence + stopCodonSequence + afterCDSAndStopCodonSequence
                    sequenceLength = len(sequence)
                    genomeLength += 300

def writeKozakOutput(outputFile, record, atgIndex, atgLocalPosition, fivePrimeOrCDSLength, sequence, kozakType, transcript, transcriptVersion, direction):
    kozak = sequence[atgIndex-9:atgIndex+10]
    strength = ""
    strengthCharacter = sequence[atgIndex-3:atgIndex-2]
    if strengthCharacter == 'A':
        strength = "VS"
    elif strengthCharacter == 'G':
        strength = 'S'
    else:
        strength = 'W'
    
    forthCharacter = sequence[atgIndex+3:atgIndex+4]
    
    outputFile.write(record["chromosome"] + " " + str(fivePrimeOrCDSLength) + " " + str(atgLocalPosition) + " " + direction + " " + kozak + " " + strength + " " + forthCharacter + " gene_id " + record["geneID"] + " transcript_id " + str(transcript) + "." + kozakType + "." + str(transcriptVersion) + " gene_type " + kozakType + " gene_name " + record["geneName"] + "\n")
    
def findKozak(utrOutputFile, cdsOutputFile, transcript, sequenceData, fivePrimeSequence, fivePrimeUTRs, cdss, stopCodon, direction):
    cdsSequence = ""
    for cds in cdss:
        cdsSequence += gencode.readSequence(sequenceData, cds['start'], cds['end'], direction)
    
    beforeFivePrimeSequence = ""
    afterCDSSequence = ""
    if direction == "+":
        if fivePrimeUTRs[0]['start'] - 9 > 0:
            beforeFivePrimeSequence = gencode.readSequence(sequenceData, fivePrimeUTRs[0]['start'] - 9, fivePrimeUTRs[0]['start'] - 1, direction)
        
        if cdss[-1]['end'] + 7 <= len(sequenceData):
            afterCDSSequence = gencode.readSequence(sequenceData, cdss[-1]['end'] + 1, cdss[-1]['end'] + 7, direction)
    else:
        if fivePrimeUTRs[-1]['start'] - 9 > 0:
            beforeFivePrimeSequence = gencode.readSequence(sequenceData, fivePrimeUTRs[-1]['start'] - 9, fivePrimeUTRs[-1]['start'] - 1, direction)
        
        if cdss[0]['end'] + 7 <= len(sequenceData):
            afterCDSSequence = gencode.readSequence(sequenceData, cdss[0]['end'] + 1, cdss[0]['end'] + 7, direction)
            
    sequence = beforeFivePrimeSequence + fivePrimeSequence + cdsSequence + afterCDSSequence
    
    atgOffsetToSearch = len(beforeFivePrimeSequence)
    
    transcriptVersion = 1
    
    while atgOffsetToSearch < len(beforeFivePrimeSequence) + len(fivePrimeSequence) - 2:
        atgIndex = sequence.find('ATG', atgOffsetToSearch)
        if atgIndex == -1:
            break
        
        if atgIndex < len(beforeFivePrimeSequence) + len(fivePrimeSequence) - 2:
            writeKozakOutput(utrOutputFile, fivePrimeUTRs[0], atgIndex, atgIndex - len(beforeFivePrimeSequence) + 1, len(fivePrimeSequence), sequence, "kozak_utr", transcript, transcriptVersion, direction)
            transcriptVersion += 1
        
        atgOffsetToSearch = atgIndex + 1
    
    transcriptVersion = 1
    atgOffsetToSearch = len(beforeFivePrimeSequence) + len(fivePrimeSequence)
    
    while atgOffsetToSearch < len(beforeFivePrimeSequence) + len(fivePrimeSequence) + len(cdsSequence) - 2:
        atgIndex = sequence.find('ATG', atgOffsetToSearch)
        if atgIndex == -1:
            break
        
        if atgIndex < len(beforeFivePrimeSequence) + len(fivePrimeSequence) + len(cdsSequence) - 2:
            writeKozakOutput(cdsOutputFile, cdss[0], atgIndex, atgIndex - len(fivePrimeSequence) - len(beforeFivePrimeSequence) + 1, len(cdsSequence), sequence, "kozak_cds", transcript, transcriptVersion, direction)
            transcriptVersion += 1
        
        atgOffsetToSearch = atgIndex + 1

def transcriptCallback(records, sequenceData, lastTranscriptID, fivePrimeUTRs, fivePrimeContent, cdss, stopCodon, direction):
    if fivePrimeContent != "":
        findKozak(fivePrimesKozakFile, cdsKozakFile, lastTranscriptID, sequenceData, fivePrimeContent, fivePrimeUTRs, cdss, stopCodon, direction)
        findKozak2(cdsKozak2File, sequenceData, lastTranscriptID, cdss, stopCodon, direction)

if __name__ == '__main__':
    arguments = gencode.getArguments()
    fivePrimesKozakFile = open(os.path.join(arguments[2], "kozak_utrs.gtf"), "w")
    fivePrimesKozakFile.write("##chr length_of_utr ATG_start strand motif strength\n")
    cdsKozakFile = open(os.path.join(arguments[2], "kozak_cds.gtf"), "w")
    cdsKozak2File = open(os.path.join(arguments[2], "kozak2_cds.gtf"), "w")
    
    gencode.run(arguments[0], arguments[1], None, transcriptCallback)
    
    cdsKozakFile.close()
    fivePrimesKozakFile.close()
    cdsKozak2File.close()