# Author: Russell Ault
# lee_to_BED converts the translation initiation data file from the Lee et al 2012
# ribosome profiling paper to a BED file format similar to that of the Fritsch et al 2012 paper,
# for use in intersection with mayur's uORFs.

# Dependencies: gencode.py

# python leetofritsch.py <path_to_annotation_file> <directory_with_chromosome_sequence_files> <lee_TIS_file> <RefSeq_Ensembl_file>

import os
import sys
import gencode

global leeToFritschFile
global tooFarleeuorfsFile


global RefSeq_dict
global lee_TIS_dict

# Finds the integer in lst closest to num
def closest_match(lst, num):
	close_num = sys.maxint
	small_diff = sys.maxint
	#print "closest position so far:\t" + str(close_num)
	#print "smallest difference so far:\t" + str(small_diff)
	for n, value in enumerate(lst):
		if abs(value - num) < small_diff:
			small_diff = abs(value - num)
			close_num = value
			#print "closest position so far:\t" + str(close_num)
			#print "smallest difference so far:\t" + str(small_diff)
	#print "smallest_difference at end:\t" + str(small_diff)
	return close_num
	
def writeUORFOutput(uORF_ID, fivePrimeUTRs, cdssAndStopCodon, startIndex, endIndex, direction, score):
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
							
			leeToFritschFile.write('\t'.join([record['chromosome'], str(startPosition - 1), str(endPosition),	# minus one converts to BED file 0-based format
											  uORF_ID, score, direction]) + '\n')
			#record['chromosome'] + '\t' + str(startPosition) + '\t' + str(endPosition) + '\t' + \
			 #					  RefSeqID + ':' + annotation + '_' + codon + '\t' + RLTM_RCHX	+ '\t' + \
			  #					 direction
			if shouldBreak:
				break

"""def conversion_parts(fivePrimeUTRs, cdssAndStopCodon, startIndex, endIndex, direction):
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
							
			return startPosition, endPosition
	#print "WHAT HAS HAPPENED!?"
	
	return 'bad', 'bad_uORF'			# Reached if the uORF coordinates don't match the transcript 
	
	#print "StartIndex: " + str(startIndex) + '\tendIndex: ' + str(endIndex)
	#print records[0]['start']
	#sys.exit(1
	#return 1, 2 """
	
# COMMENT HERE?
def convert(sequenceData, lastTranscriptID, fivePrimeSequence, fivePrimeUTRs, cdss, stopCodon, direction):
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
	
		
	EnsemblID = lastTranscriptID.split('.')[0]  # transcript stripped of version number to match RefSeq dictionary transcript
	uORF_number = 1
	while lee_TIS_dict.get((EnsemblID, uORF_number)):
		RefSeqID, position_to_aTIS, annotation, RLTM_RCHX, codon, uORF_number = lee_TIS_dict.get((EnsemblID, uORF_number))
		position = position_to_aTIS + sequence.find(cdsSequence)
		
		# Now find the closest codon in the transcript sequence
		offSet = 0
		lst_pos = []
		codon_position = sequence.find(codon, offSet)
		offSet = codon_position + 1
		while codon_position != -1:
			lst_pos.append(codon_position)
			codon_position = sequence.find(codon, offSet)
			offSet = codon_position + 1
		best_position = closest_match(lst_pos, position)
		
		
		#startPosition, endPosition = conversion_parts(fivePrimeUTRs, cdssAndStopCodon, best_position - 1, best_position + 3, direction)
		
		record = fivePrimeUTRs[0]
		uORF_ID = RefSeqID + ':' + annotation + '_' + codon
		
		if abs(best_position - position) <= 5:	   # 5 is an arbitrary threshold	
			writeUORFOutput(uORF_ID, fivePrimeUTRs, cdssAndStopCodon, best_position, best_position + 3, direction, RLTM_RCHX)
			# 5 is an arbitrary threshold for matching based on me seeing that most of the elements were off by 1 or 0, which is what I expect
			#leeToFritschFile.write(record['chromosome'] + '\t' + str(startPosition) + '\t' + str(endPosition) + '\t' + \
			 #					  RefSeqID + ':' + annotation + '_' + codon + '\t' + RLTM_RCHX	+ '\t' + \
			  #					 direction + '\n')
		else:
			tooFarleeuorfsFile.write(record['chromosome'] + '\t' + str(best_position-1) + '\t' + str(best_position + 3) + '\t' + \
									 RefSeqID + ':' + annotation + '_' + codon + '\t' + RLTM_RCHX	+ '\t' + \
									 direction + '\t' + str(best_position - position) + '\n')   
		uORF_number += 1
		

def transcriptCallback(records, sequenceData, lastTranscriptID, fivePrimeUTRs, fivePrimeContent, cdss, stopCodon, direction):
	if fivePrimeContent != "":
		convert(sequenceData, lastTranscriptID, fivePrimeContent, fivePrimeUTRs, cdss, stopCodon, direction)
		

if __name__ == '__main__':
	transdict = {'Ensembl_Gene_ID': 0, 'Ensembl_Transcript_ID': 1, 'chrom': 2, 'Transcript_Start': 3,
				 'Transcript_End': 4, 'RefSeq_mRNA_ID': 5}
				 
	if len(sys.argv) != 5:
		print "Incorrect number of arguments"
		sys.exit(1)
	
	RefSeq_Ensembl_dict = {}
	with open(sys.argv[4]) as RefSeq_Ensembl_file:
		RefSeq_Ensembl_file.readline()	  # skip intro line
		for line in RefSeq_Ensembl_file:
			fields = line.split()
			RefSeq_Ensembl_dict[fields[transdict['RefSeq_mRNA_ID']]] = fields[transdict['Ensembl_Transcript_ID']]
										   
	# Dictionary for creating unique uORF identifiers for multiple uORF transcripts
	Ensembl_count_dict = {}
	for EnsemblID in RefSeq_Ensembl_dict.values():
		Ensembl_count_dict[EnsemblID] = 1
	lee_TIS_dict = {}
	with open(sys.argv[3]) as lee_TIS_file:
		firstline = lee_TIS_file.readline()				 # replicate intro line
		for line in lee_TIS_file:
			if "3'UTR" not in line:
				fields = line.split()
				RefSeqID = fields[0]
				position_to_aTIS = int(fields[4])
				annotation = fields[5]
				RLTM_RCHX = fields[6]
				codon = fields[12]
				EnsemblID = RefSeq_Ensembl_dict.get(RefSeqID)

				if EnsemblID:
					uORF_number = Ensembl_count_dict.get(EnsemblID)
					lee_TIS_dict.setdefault((EnsemblID, uORF_number),
											(RefSeqID,
											position_to_aTIS,
											annotation,
											RLTM_RCHX,
											codon,
											uORF_number))
					Ensembl_count_dict[EnsemblID] = uORF_number + 1

	leeToFritschFile = open("lee_uORF2.bed", "w") 
	leeToFritschFile.write('#chrom\tchromStart\tchromEnd\tname\tRLTM-RCHX\tstrand\n')			   # header line
   
	tooFarleeuorfsFile = open("lee_toofar_uORFs.bed", "w")
	tooFarleeuorfsFile.write("#These are the bad uorfs that don't match too far from corresponding codons in the GENCODE transcripts\n" +
							 "#chrom\tchromStart\tchromEnd\tname\tRLTM-RCHX\tstrand\tdistance-from-nearest-GENCODE-codon\n")
	
	gencode.run(sys.argv[1], sys.argv[2], None, transcriptCallback)
	
	leeToFritschFile.close()
	tooFarleeuorfsFile.close()
