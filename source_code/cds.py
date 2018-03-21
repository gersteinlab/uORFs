#Coder: Russell Ault
#Program: Writes out cds sequences
#Dependencies: gencode.py


# THERE ARE ERRORS, NOT PRINTING OUT sequences for many elements
# I PROBABLY OUGHT TO CHECK to SEE WHETHER THE CDS SEQUENCE IS EMPTY OR NOT BEFORE I WRITE IT, IN ALL THE FILES

import os
import gencode

global cdsSequenceFile
global cdsAnnotationFile

def newChromosomeCallback(newChromosome):
    cdsSequenceFile.write(">" + newChromosome + "\n")

def transcriptCallback(records, sequenceData, lastTranscriptID, fivePrimeUTRs, fivePrimeContent, cdss, stopCodon, direction):
    cdsSequence = ""
    for cds in cdss:
        cdsSequence += gencode.readSequence(sequenceData, cds['start'], cds['end'], direction)
        cdsAnnotationFile.write(cds['line'])
    
    stopCodonSequence = ""
    if stopCodon:
        stopCodonSequence += gencode.readSequence(sequenceData, stopCodon['start'], stopCodon['end'], direction)
        
        #print 'error! No transcript ID for ' + str(lastTranscriptID)
    
    if cdsSequence:    
        cdsSequenceFile.write(">" + str(records[0]['geneID']) + "|" + str(lastTranscriptID) + "|" + str(records[0]["geneName"]) + "\n")
        cdsSequenceFile.write(cdsSequence + stopCodonSequence + "\n")

if __name__ == '__main__':
    arguments = gencode.getArguments()
    
    cdsSequenceFile = open(os.path.join(arguments[2], "CDSs.fa"), "w")
    cdsAnnotationFile = open(os.path.join(arguments[2], "CDSs.gtf"), "w")
    
    gencode.run(arguments[0], arguments[1], newChromosomeCallback, transcriptCallback)
    
    cdsSequenceFile.close()
    cdsAnnotationFile.close()
