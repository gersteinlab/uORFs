import os, sys
import gencode

def transcriptCallback(records, sequenceData, lastTranscriptID, fivePrimeUTRs, fivePrimeContent, cdss, stopCodon, direction):
    #Scan intron's, ignore first cds block
    for cds in cdss[1:]:
        if direction == '+':
            startPosition = cds['start']-6
            endPosition = cds['start']+2
        elif direction == '-':
            startPosition = cds['end']-2
            endPosition = cds['end']+6
        
        #Not in bounds of sequence data, skip
        if endPosition > len(sequenceData):
            print 'once'
            print cds['start']
            print cds['end']
            print cds
            print len(sequenceData)
            continue
            
        sequence = gencode.readSequence(sequenceData, startPosition, endPosition, direction)
        
        #NAGNAGNAG
        
        lowerPosition = 'NA'
        higherPosition = 'NA'
        if sequence[1:3][0] == 'A' or sequence[1:3][1] == 'G':
            if direction == '+':
                lowerPosition = str(startPosition + 1)
            elif direction == '-':
                higherPosition = str(endPosition - 2)
        if sequence[7:9][0] == 'A' or sequence[7:9][1] == 'G':
            if direction == '+':
                higherPosition = str(startPosition + 7)
            elif direction == '-':
                lowerPosition = str(endPosition - 8)
        
        if "AG" == sequence[4:6]:
            spliceAcceptorCase = 'SpliceAcceptorIsAG'
        else:
            spliceAcceptorCase = 'SpliceAcceptorIsNotAG'
        
        agCase = None
        if "AG" == sequence[1:3] or "AG" == sequence[7:9]:
            agCase = 'FoundAGOnSides'
        elif sequence[1:3][0] == 'A' or sequence[1:3][1] == 'G' or sequence[7:9][0] == 'A' or sequence[7:9][1] == 'G':
            agCase = 'FoundEitherAOrGOnSides'
        
        if agCase:
            nagAnnotationFile.write("\t".join([records[0]['chromosome'], 'mayur', 'NAGNAG', records[0]['geneID'], records[0]['geneName'], records[0]['geneType'], lastTranscriptID, records[0]['transcriptType'], direction, str(startPosition), str(endPosition), sequence, lowerPosition, higherPosition, spliceAcceptorCase, agCase]) + "\n")

if __name__ == '__main__':
    arguments = gencode.getArguments()
    
    nagAnnotationFile = open(os.path.join(arguments[2], "nagnag_cases_test_ignore.gtf"), "w")
    
    gencode.run(arguments[0], arguments[1], None, transcriptCallback)
    
    nagAnnotationFile.close()
    