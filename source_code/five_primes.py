#Coder: Mayur Pawashe
#Program: Writes out five primes
#Dependencies: gencode.py

import os
import gencode

global fivePrimesSequenceFile
global fivePrimesAnnotationFile

def newChromosomeCallback(newChromosome):
    fivePrimesSequenceFile.write(">" + newChromosome + "\n")

def transcriptCallback(records, sequenceData, lastTranscriptID, fivePrimeUTRs, fivePrimeContent, cdss, stopCodon, direction):
    if fivePrimeContent != "":
        fivePrimesSequenceFile.write(">" + str(records[0]['geneID']) + "|" + str(lastTranscriptID) + "|" + str(records[0]["geneName"]) + "\n")
        fivePrimesSequenceFile.write(fivePrimeContent + "\n")
    
    for fivePrime in fivePrimeUTRs:
        fivePrimesAnnotationFile.write(fivePrime['line'])

if __name__ == '__main__':
    arguments = gencode.getArguments()
    
    fivePrimesSequenceFile = open(os.path.join(arguments[2], "five_primes.fa"), "w")
    fivePrimesAnnotationFile = open(os.path.join(arguments[2], "five_primes.gtf"), "w")
    
    gencode.run(arguments[0], arguments[1], newChromosomeCallback, transcriptCallback)
    
    fivePrimesSequenceFile.close()
    fivePrimesAnnotationFile.close()
