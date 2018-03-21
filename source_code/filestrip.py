# Author: Russell Ault
# filestrip extracts RefSeq trancsript IDs from the puroTranscript file of the 
# Fritsch et al 2012 study and strips non uORF annotations from the puroClassified file
# of the same study.

# python filestrip.py -[tr,u, s, XY, r, rev, su, rX, fri] <absolute-path-to-file> [<optional-second-fileof-uORFs>]
# -tr is for extracting transcripts; -u is for stripping non uORFs
# -s swaps columns 1 and 2 with columns 3 and 4 of the input fields
# -rev currently reverts this swapping and adds 'chr' to each chromosome number
# -XY turns chromosomes X and Y into chromosomes 23 and 24 in the input file
# -rX reverts chromosomes 23 and 24 back to X and Y
# -fri is for extracting uORF IDs that are translated according to Fritsch et al.(2012) data

import os
import sys

global studySource          # lee or fritsch
global uORFdict
global transdict            # transcript file dictionary for matched RefSeq and ENSEMBL Transcript IDs

def extractTranscripts(file):
    file.readline()         # skips intro line
    for line in file:
        fields = line.split()
        RefSeqID = fields[uORFdict['name']]
        stripRefSeqID = RefSeqID[:RefSeqID.find(':')]
        print stripRefSeqID

def stripNonUORF(file):
    for line in file:
        if 'ORF' in line:
            print line,

def sortChange(file):
    for line in file:
        fields = line.split()
        print fields[2] + '\t' + fields[3] + '\t' + fields[0] + '\t' + fields[1] + '\t' + fields[4] + '\t' + fields[5]
def revert(file):
    for line in file:
        fields = line.split()
        print fields[2] + '\t' + fields[3] + '\t' + "chr" + fields[0] + '\t' + fields[1] + '\t' + fields[4] + '\t' + fields[5]
            
def alterXY(file):
    for line in file:
        fields = line.split()
        if fields[0] == 'X':
            print '23' + '\t' + fields[1]+ '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5]
        elif fields[0] == 'Y':
            print '24' + '\t' + fields[1]+ '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5]
        else :
            print fields[0] + '\t' + fields[1]+ '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5]

def reverseXY(file):
    for line in file:
        fields = line.split()
        if fields[0] == '23':
            print 'X' + '\t' + fields[1]+ '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5]
        elif fields[0] == '24':
            print 'Y' + '\t' + fields[1]+ '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5]
        else :
            print fields[0] + '\t' + fields[1]+ '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5]
            
def eliminateNonUORFtranscripts(transcriptfile, uORFfile):
    uORFfile.readline()
    uORFfile.readline()                 # skips first two intro lines
    transcriptfile.readline()           # skips first intro line
    
    ID_dict = {}
    for line in uORFfile:
        uORFline = line.split()
        # Puts the transcript ID for each uORF in the ID dictionary
        ID_dict.setdefault(uORFline[uORFdict['name']].rstrip(':ovuORF_ATCG_conserved'), True)
    for line in transcriptfile:
        transline = line.split()
        if ID_dict.get(transline[transdict['RefSeq_mRNA_ID']]): # None is returned if the ID is not in the dictionary
            print line,
    fileUORF = open("help_compare_uORF", "w") 
    for key in ID_dict.keys():
        fileUORF.write(key + '\n')
    fileUORF.close()

def fritschUORF(translated_uORF_file, codon):
    line = translated_uORF_file.readline()
    if studySource == 'fritsch':
        while line:
            if 'Intersection' in line and "_" + codon in line:
                line = translated_uORF_file.readline()
                fields = line.split()
                print fields[uorf_gtf_dict['Ensembl_Transcript_ID']].strip('";UTRCDS.')
            line = translated_uORF_file.readline()
            
    elif studySource == 'lee':
        while line:
            if 'Intersection' in line and "5'UTR_" + codon in line:
                line = translated_uORF_file.readline()
                fields = line.split()
                print fields[uorf_gtf_dict['Ensembl_Transcript_ID']].strip('";UTRCDS.')
            line = translated_uORF_file.readline()

    elif studySource == 'gao':
        while line:
            if 'Intersection' in line and ":" + codon in line:
                line = translated_uORF_file.readline()
                fields = line.split()
                print fields[uorf_gtf_dict['Ensembl_Transcript_ID']].strip('";UTRCDS.')
            line = translated_uORF_file.readline()
    

if __name__ == '__main__':
    opt = sys.argv[1]
    filepath = sys.argv[2]
    if 'lee' in filepath:
        studySource = 'lee'
    elif 'fritsch' in filepath:
        studySource = 'fritsch'
    elif 'gao' in filepath:
        studySource = 'gao'
    uORFdict = {'chrom': 0, 'chromStart': 1, 'chromEnd': 2, 'name': 3, 'score': 4,
                   'strand': 5, 'thickStart': 6, 'thickEnd': 7, 'reserved': 8, 'blockCount': 9,
                   'blockSizes': 10, 'chromStarts': 11}
    transdict = {'Ensembl_Gene_ID': 0, 'Ensembl_Transcript_ID': 1, 'chrom': 2, 'Transcript_Start': 3,
                 'Transcript_End': 4, 'RefSeq_mRNA_ID': 5}
    uorf_gtf_dict = {'chrom': 0, 'chromStart': 3, 'chromEnd': 4, 'strand': 6, 'Ensembl_Transcript_ID': 11}
    with open(filepath) as file:
        if opt == '-tr':
            extractTranscripts(file)
        elif opt == '-u':
            stripNonUORF(file)
        elif opt == '-s':
            sortChange(file)
        elif opt == '-XY':
            alterXY(file)
        elif opt == '-rX':
            reverseXY(file)
        elif opt == '-rev':             # 'rev' for revert
            revert(file)
        elif opt == '-su':
            if len(sys.argv) == 4:
                with open(uORFfilepath) as uORFfile:
                    eliminateNonUORFtranscripts(file, uORFfile)
        elif opt == '-fri':
            codon = sys.argv[3][1:]
            fritschUORF(file, codon)
        else:
            print >> sys.stderr, "No option selected for filestrip.py"
        
            
        
