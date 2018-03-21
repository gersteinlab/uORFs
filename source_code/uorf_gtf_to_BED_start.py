#! /usr/bin/env python

# Author: Russell Ault

# uorf_gtf_to_BED_start.py converts a uORF GTF file to a BED file containing start codon
# coordinates for each uORF. This file will be used in intersecting properly with experimental uORFs

# python uorf_gtf_to_BED_start.py CODON <uORF-gtf-file> <out-put-directory>

from uorf_gtf_to_GERP import codon_base_positions

import os
import sys

global Codon    # start codon for the uORFs


# COMMENT HERE!e
def output_to_BED(fields_lst, uORF_ID, uORF_BED_file):
    chromosome = fields_lst[0][0]
    strand = fields_lst[0][6]
    start_positions = [(position-1) for position in codon_base_positions(fields_lst, strand, True)]    # converts to 0-based BED format
            
    # Output to BED file
    for base_position in start_positions:
        uORF_BED_file.write('\t'.join([chromosome, str(base_position), str(base_position+1),
                                        uORF_ID, '0', strand]) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print "Usage: python uorfgtfToGERP.py CODON <uORF-gtf-file> <out-put-directory>"
        print "Incorrect number of arguments supplied"
        print "Exiting..."
        sys.exit(1)
    
    Codon = sys.argv[1]
    uORF_GTF_filename = sys.argv[2]
    output_path = sys.argv[3]
    
    uORF_GTF_file = open(uORF_GTF_filename, 'r')
    uORF_BED_file = open(os.path.join(output_path, "uORFs_%s_start.bed" % (Codon)), "w")
    
    # Header line
    uORF_BED_file.write('#chr\tchromStart\tchromEnd\tuORF-ID\tfake-score\tstrand\n')
    
    prevUORF_ID = None
    fields_lst = []
    line = uORF_GTF_file.readline()  
    while line:
        fields = line.split()
        uORF_ID = fields[11].rstrip(';').strip('"').rstrip('.UTRCDSstop_codon')
        if uORF_ID != prevUORF_ID and prevUORF_ID:                # first line of a new uORF descriptor
            output_to_BED(fields_lst, prevUORF_ID, uORF_BED_file)
            fields_lst = []
            
        fields_lst.append(fields)
        prevUORF_ID = uORF_ID
        line = uORF_GTF_file.readline()
        
        if not line:                        # Process last uORF
            output_to_BED(fields_lst, prevUORF_ID, uORF_BED_file)
            break
    
    # Close files
    uORF_GTF_file.close()
    uORF_BED_file.close()
