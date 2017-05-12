# Author: Russell Ault
# intersectTIS finds the matches between the uORFs in the UCSC_uorfs file
# and our uORFs.gtf file, using a transcripts file as dictionary between RefSeq
# and ENSEMBL transcripts IDs

import os
import sys

global uORFdict             # dictionary for fields of UCSC uORF file
global transdict            # transcript file dictionary for matched RefSeq and ENSEMBL Transcript IDs
global gtfdict              # dictionary for fields of uORFs.gtf file
global stringNum            # stringency number for determining how far apart two open reading frames are allowed to be

# findString_in_lst returns the index of "string" in a list of strings.
# If the string is not found, None is returned.
def findString_in_lst(lst, string):
    for n, elem in enumerate(lst):
        if elem.startswith(string):
            return n    
    return None

def intersect(Gencode_uorfs, fritsch_uorfs, RefSeq_Ensembl, matchID=False):
    uORFs_UCSC_dict = {}     # This dictionary holds keys of tuple (chromosome, chromStart, chromEnd, strand) with values of name classification from the UCSC uORFs file
    fritsch_uorfs.readline() # skips intro line
    
    # Filter uORFs.gtf file, printing matches according to start location with stringency error
    RefSeq_Ensemble_dict = {}
    if matchID:
        #RefSeq_Ensembl.readline()
        RefSeq_Ensembl.readline() # skip intro line
        for line in RefSeq_Ensembl:
            sline = line.split()
            RefSeq_Ensemble_dict[sline[5]] =  sline[1] # key is RefSeqID, value is Ensembl Transcript ID
            #print 'RefSeq: ' + sline[5] + '\tEnsembl: ' + sline[1] # DEBUGGING CODE

    # Make uORF dictionary for fast filtering of uORFs.gtf file
    for line in fritsch_uorfs:
        uORFline = line.split()
        EnsemblID = RefSeq_Ensemble_dict.get(uORFline[uORFdict['name']].split(':')[0])
        #if EnsemblID == None:
            #print 'ERROR, NO REFSEQ ID' # This line should never be run
            #sys.exit(1)
        if uORFline[uORFdict['strand']] == '+' and EnsemblID != None:
            uORFs_UCSC_dict[(EnsemblID, uORFline[uORFdict['chrom']], int(uORFline[uORFdict['chromStart']]), '+')] = uORFline[uORFdict['name']]
        elif uORFline[uORFdict['strand']] == '-' and EnsemblID != None:
            uORFs_UCSC_dict[(EnsemblID, uORFline[uORFdict['chrom']], int(uORFline[uORFdict['chromEnd']]), '-')] = uORFline[uORFdict['name']]
            
    # chrompos = 'chromStart' #or 'chromEnd'
    inUORF = False
    prevID = None
    
    for line in Gencode_uorfs:
        genUORFline = line.split()
        transID = genUORFline[gtfdict['Ensembl_Transcript_ID']].rstrip('UTRCDSstop_codon";') # Ensembl_Transcript_ID with leading "
        #transID = rawtransID[:rawtransID.find('uORF')] + rawtransID.rstrip('UTRCDSstop_codon')[-1]
        if transID != prevID:
            
            inUORF = False
    
            if genUORFline[gtfdict['strand']] == '+':
                chrompos = 'chromStart'
                strand = '+'
            elif genUORFline[gtfdict['strand']] == '-':
                chrompos = 'chromEnd'
                strand = '-'
            modtransID = transID.strip('"')[:transID.find('.')-1]
    
            for n in xrange(0, stringNum + 1):
                poskey = (modtransID, genUORFline[gtfdict['chrom']], int(genUORFline[gtfdict[chrompos]]) + int(n)-1, strand) # Minus one is to account for RefSeq 0-based numbering
                negkey = (modtransID, genUORFline[gtfdict['chrom']], int(genUORFline[gtfdict[chrompos]]) - int(n), strand)
                
                posRefSeq = uORFs_UCSC_dict.get(poskey)
                negRefSeq = uORFs_UCSC_dict.get(negkey)
                if posRefSeq:
                    if not matchID or RefSeq_Ensemble_dict.get(posRefSeq[:posRefSeq.find(':')]) == modtransID :# matchID option off or ID in dictionary
                        print 'Intersection with Fritsch uORFs at:\t' + poskey[1] + '\t' + chrompos + ': ' + str(poskey[2]) \
                              + '\tstrand: ' + poskey[3] + '\tRefSeqID_and_classification: ' + posRefSeq
                        print line,
                        inUORF = True
                        break
                elif negRefSeq:
                    if not matchID or RefSeq_Ensemble_dict.get(negRefSeq[:negRefSeq.find(':')]) == modtransID:
                        print 'Intersection with Fritsch uORFs at:\t' + negkey[1] + '\t' + chrompos + ': ' + str(negkey[2]) \
                              + '\tstrand: ' + negkey[3] + '\tRefSeqID_and_classification: ' + negRefSeq
                        print line,
                        inUORF = True
                        break
            
        elif inUORF:
            print line,
        
        prevID = transID            
        
        #if inUORF:
            #if transID != prevID:
             #   inUORF = False
            #else
                
        #if genUORFline[gtfdict['strand']] == '+':
         #   chrompos = 'chromStart'
        #else:
        #    chrompos = 'chromEnd'
        #if not inUORF:

if __name__ == '__main__':
    uORFdict = {'chrom': 0, 'chromStart': 1, 'chromEnd': 2, 'name': 3, 'score': 4,
                'strand': 5, 'thickStart': 6, 'thickEnd': 7, 'reserved': 8, 'blockCount': 9,
                'blockSizes': 10, 'chromStarts': 11}
    transdict = {'Ensembl_Gene_ID': 0, 'Ensembl_Transcript_ID': 1, 'chrom': 2, 'Transcript_Start': 3,
                 'Transcript_End': 4, 'RefSeq_mRNA_ID': 5}
    gtfdict = {'chrom': 0, 'chromStart': 3, 'chromEnd': 4, 'strand': 6, 'Ensembl_Transcript_ID': 11}
    if len(sys.argv) !=  6:
        print "Incorrect number of arguments"
        print "Make sure -m flag is included"
        sys.exit(1)
    stringNum = int(sys.argv[1])
    uORFs_gtf_file = sys.argv[2]
    UCSC_uorfs_file = sys.argv[3]
    RefSeq_ENSEMBL_file = sys.argv[4]
    if len(sys.argv) == 6:
        matchID = True
    else:
        matchID = False
    
    with open(uORFs_gtf_file) as our_uORFs_file: 
        with open(UCSC_uorfs_file) as fritsch_uORFs_file:
            with open(RefSeq_ENSEMBL_file) as dict_file:
                intersect(our_uORFs_file, fritsch_uORFs_file, dict_file, matchID)
