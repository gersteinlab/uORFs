import os
import sys
import operator

global uORFs_codons
global Codon_Frequency
global Rarest_Codons
global secstruct_start_positions
global secstruct_uORF_end_positions
global RNAseq_window_length	

def makedict(seq_file, ID_field=1, chrline=False, nintro=0, uORF_option=False):
	seq_dict = {}			   # This dictionary holds transcript ID keys and the DNA sequence as the value
	for n in xrange(0, nintro):
		seq_file.readline()	 # skip intro line
	seq_line = seq_file.readline()
	while seq_line:
		if chrline:
			if seq_line.startswith('>chr') or seq_line.startswith('#'):						# skip new chromosome lines
				seq_line = seq_file.readline()
			if not seq_line:
				break
		fields = seq_line.split('|')
		#print "ID_field:\t" + str(ID_field)
		#if len(fields) < 3:
			#print seq_line
			#sys.exit(1)
		trans_ID = fields[ID_field]
		if uORF_option:
			seq_dict.setdefault(trans_ID,
								True)
			seq_line = seq_file.readline()
			seq_line = seq_file.readline()
			continue
		seq_line = seq_file.readline()
		sequence = seq_line.split()
		seq_dict[trans_ID] = sequence[0]
		seq_line = seq_file.readline()
	return seq_dict


def make_profile_dict(IDs_file, key_pos, value):
	profile_dict = {}
	for line in IDs_file:
		profile_dict[line.split()[key_pos]] = value
	return profile_dict
		

def make_cons_dict(file, key_pos, nvalues, nintro_lines):
	cons_dict = {}
	while nintro_lines > 0:
		file.readline()			# skip intro lines
		nintro_lines -= 1
	for line in file:
		fields = line.split()
		if nvalues == 1:
			cons_dict[fields[key_pos]] = fields[1]
		else:
			cons_dict[fields[key_pos]] = tuple(fields[1:])
	return cons_dict

def make_GTEX_dict(file, key_pos):
	GTEX_dict = {}
	file.readline()  # skip intro line
	currentline = file.readline()
	while currentline:
		if currentline == '\n':
			break
		if not currentline:
			break
		fields = currentline.split()
		GTEX_dict[fields[key_pos]] = tuple(fields[1:])
		currentline = file.readline()
	return GTEX_dict

def make_SNP_dict(file, key_pos):
	SNP_dict = {}
	file.readline()  # skip intro line
	currentline = file.readline()
	while currentline:
		if currentline == '\n':
			break
		if not currentline:
			break
		fields = currentline.split()
		SNP_dict[fields[key_pos]] = tuple(fields[1:])
		currentline = file.readline()
	return SNP_dict

def make_Noderer_Score_dict(file):
	Noderer_dict = {}
	file.readline() #skip intro line
	file.readline() #skip second intro line
	currentline = file.readline()
	while currentline:
		if currentline == '\n':
			break
		if not currentline:
			break
		fields = currentline.split()
		Noderer_dict[fields[0]] = fields[1]
		currentline = file.readline()
	return Noderer_dict
	
"""
	{}
	GERP_start_stop.readline()					# skip intro line
	for line in GERP_start_stop:
		fields = line.split()
		uORF_ID, start_score, stop_score = fields[0], fields[1], fields[2]
		GERP_score_dict[uORF_ID] = (start_score, stop_score)
	
	GERP_elements_dict = {}
	GERP_elements.readline()					# skip intro line
	for line in GERP_elements:
		fields = line.split()
		uORF_ID, base_percent = fields[0], fields[1]
		GERP_elements_dict[uORF_ID] = base_percent
	
	PhyloCSF_dict = {}
	for line in PhyloCSF:
		fields = line.split()
		uORF_ID, cons_score = fields[0], fields[1] """


def percent_base(DNA_seq, base):
	return float(DNA_seq.count(base)) / float(len(DNA_seq))


def find_float(line):
	float_lst = []
	in_float = False
	for char in line:
		if not in_float and char == '-' or char == '0':
			in_float = True
			float_lst.append(char)
		elif in_float:
			if char == ')':
				break
			else:
				float_lst.append(char)
	return ''.join(float_lst)

 
# COMMENT HERE!	
def get_thermo_mean(energy_file):
	dG_values = []
	
	#energy_file = open(filename, 'r')
	line = energy_file.readline()
	is_odd_line = True
	while line:
		if is_odd_line:
			line = energy_file.readline()
			is_odd_line = False
			continue
		
		dG_values.append(find_float(line))
		is_odd_line = True
		line = energy_file.readline()
	energy_file.close()
	
	dG_mean = 0.0
	for value in dG_values:
		dG_mean += float(value)
	if dG_mean != 0.0:
		dG_mean /= float(len(dG_values))
	return dG_mean


def get_RNA_dGs(mRNA_transcript, mRNA_len, wind_len, DNA_index, Start=True):
	RNAseqs = []
	RNApath = '/home2/rca26/Programs/bin/'
	positions = secstruct_start_positions if Start else secstruct_uORF_end_positions
	for start in positions:
		index = start + DNA_index
		if index > 0 and index + wind_len < mRNA_len: 
			RNAseqs.append(mRNA_transcript[index:index + wind_len])
		else:
			RNAseqs.append('N_A')
	dGs = []
	for seq in RNAseqs:
		if seq == 'N_A':
			dGs.append('N_A')
		else:
			dGs.append(get_thermo_mean(os.popen('echo %s | %sRNAfold --noPS' % (seq, RNApath), 'r')))
	return [str(i) for i in dGs]
	

def get_uORF_count(uORFs_ID_dict, codon, uORF_sampler_stripped):
	ID_number = 1
	uORFs_in_transcript = 0
	uORF_sampler = uORF_sampler_stripped + codon + '.'
	ID_sampler = uORF_sampler + str(ID_number)
	while uORFs_ID_dict.get(ID_sampler):
		uORFs_in_transcript += 1
		ID_number += 1
		ID_sampler = uORF_sampler + str(ID_number)
	return uORFs_in_transcript
	

def statoutput(uORF_ID, uORF_DNA_seq, uORFs_ID_dict, uORFs_ATG_ID_dict, uORFs_CTG_ID_dict, uORFs_TTG_ID_dict, uORFs_GTG_ID_dict, uORFs_ATC_ID_dict, uORFs_ATT_ID_dict, uORFs_ATA_ID_dict, uORFs_ACG_ID_dict, uORFs_AAG_ID_dict, uORFs_AGG_ID_dict, fiveprime_seq_dict, 
			   CDS_seq_dict, GERP_score_dict, GERP_elements_dict, PhyloCSF_dict, GTEX_dict, SNP_dict, Noderer_dict, fritsch_uORF_dict=False, 
			   lee_uORF_dict=False, gao_uORF_dict=False, lee_aTIS_dict=False, fritsch_aTIS_dict=False, gao_aTIS_dict=False):
	
	transcript_ID = uORF_ID.rstrip('0123456789.').rstrip('uORF._ATGC')
	fiveprime_DNA = fiveprime_seq_dict.get(transcript_ID)
	CDS_DNA = CDS_seq_dict.get(transcript_ID)
	
	if fiveprime_DNA and CDS_DNA:
		mRNA_transcript = fiveprime_DNA + CDS_DNA 
		CDS_DNA_start = mRNA_transcript.index(CDS_DNA)
		uORF_DNA_start = mRNA_transcript.index(uORF_DNA_seq)
		cap_to_uORF_distance = uORF_DNA_start
		uORF_length = len(uORF_DNA_seq)
		CDSstart_to_uORFstart = CDS_DNA_start - uORF_DNA_start
		CDSstart_to_uORFend = CDS_DNA_start - (uORF_DNA_start + uORF_length)
		percent_A = percent_base(uORF_DNA_seq, 'A')
		percent_T = percent_base(uORF_DNA_seq, 'T')
		percent_G = percent_base(uORF_DNA_seq, 'G')
		percent_C = percent_base(uORF_DNA_seq, 'C')

		#Note: the Noderer context score below.

		uORF_Noderer_context = mRNA_transcript[(uORF_DNA_start-6):(uORF_DNA_start)].upper() + 'AUG' + uORF_DNA_seq[3:5].upper()
		uORF_Noderer_context = uORF_Noderer_context.replace('T', 'U')
		if not Noderer_dict.get(uORF_Noderer_context):
			uORF_Noderer_score = 0
		else:
			uORF_Noderer_score = Noderer_dict.get(uORF_Noderer_context)
		
		uORF_context = ('strong' if (uORF_DNA_seq[3].upper() == 'G'
									 and (mRNA_transcript[uORF_DNA_start-3].upper() == 'A'
										  or mRNA_transcript[uORF_DNA_start-3].upper() == 'G'))
								 else ('weak' if (uORF_DNA_seq[3].upper() == 'G'
												  or (mRNA_transcript[uORF_DNA_start-3].upper() == 'A'
										  			  or mRNA_transcript[uORF_DNA_start-3].upper() == 'G'))
										  	  else 'none'))
		UTR_size = len(fiveprime_DNA)
		
		uORF_codon = uORF_DNA_seq[:3]
		uORF_sampler = uORF_ID.rstrip('.0123456789').rstrip('ATCG')

		ATG_uORFs_in_transcript = get_uORF_count(uORFs_ATG_ID_dict, 'ATG', uORF_sampler)
		CTG_uORFs_in_transcript = get_uORF_count(uORFs_CTG_ID_dict, 'CTG', uORF_sampler)
		TTG_uORFs_in_transcript = get_uORF_count(uORFs_TTG_ID_dict, 'TTG', uORF_sampler)
		GTG_uORFs_in_transcript = get_uORF_count(uORFs_GTG_ID_dict, 'GTG', uORF_sampler)
		ATC_uORFs_in_transcript = get_uORF_count(uORFs_ATC_ID_dict, 'ATC', uORF_sampler)
		ATT_uORFs_in_transcript = get_uORF_count(uORFs_ATT_ID_dict, 'ATT', uORF_sampler)
		ATA_uORFs_in_transcript = get_uORF_count(uORFs_ATA_ID_dict, 'ATA', uORF_sampler)
		ACG_uORFs_in_transcript = get_uORF_count(uORFs_ACG_ID_dict, 'ACG', uORF_sampler)
		AAG_uORFs_in_transcript = get_uORF_count(uORFs_AAG_ID_dict, 'AAG', uORF_sampler)
		AGG_uORFs_in_transcript = get_uORF_count(uORFs_AGG_ID_dict, 'AGG', uORF_sampler)

		total_start_codons = (ATG_uORFs_in_transcript + CTG_uORFs_in_transcript)

		same_codon_uORFs_in_transcript = get_uORF_count(uORFs_ID_dict, uORF_codon, uORF_sampler)
		
		codon_frequency_statfile = open("~/ribosome_profiling/uORF_frequency_stats_" + uORF_codon + ".txt", 'r')
		codon_frequency_statfile.readline() #skip the header line
		codon_frequency = codon_frequency_statfile.readline().split()[1]

		"""for codon in uORFs_codons:
			ID_number = 1
			uORF_sampler_curr = uORF_sampler + codon + '.'
			ID_sampler = uORF_sampler_curr + str(ID_number)
			while uORFs_ID_dict.get(ID_sampler):
				uORFs_in_transcript += 1
				if codon == uORF_codon:
					same_codon_uORFs_in_transcript += 1
				ID_number += 1
				ID_sampler = uORF_sampler_curr + str(ID_number)"""
		
		# Translation Initiation by Ribosome Profiling (Lee and Fritsch studies so far [2012])		
		is_fritsch_translated = 'no'
		is_lee_translated = 'no'
		is_gao_translated = 'no'
		is_lee_aTIS_translated = 'N_A'
		is_fritsch_aTIS_translated = 'N_A'
		is_gao_aTIS_translated = 'N_A'
		
		if fritsch_uORF_dict and lee_uORF_dict and gao_uORF_dict:
			is_fritsch_translated = fritsch_uORF_dict.get(uORF_ID, 'no')
			is_lee_translated = lee_uORF_dict.get(uORF_ID, 'no')
			is_gao_translated = gao_uORF_dict.get(uORF_ID, 'no')
		elif fritsch_uORF_dict:
			is_fritsch_translated = 'yes'
		elif lee_uORF_dict:
			is_lee_translated = 'yes'
		elif gao_uORF_dict:
			is_gao_translated = 'yes'
			
		if lee_aTIS_dict.get(transcript_ID.rstrip('0123456789').rstrip('.')):
			is_lee_aTIS_translated = 'yes'
		elif is_lee_translated == 'yes':
			is_lee_aTIS_translated = 'no'
		if fritsch_aTIS_dict.get(transcript_ID.rstrip('0123456789').rstrip('.')):
			is_fritsch_aTIS_translated = 'yes'
		elif is_fritsch_translated == 'yes':
			is_fritsch_aTIS_translated = 'no'
		if gao_aTIS_dict.get(transcript_ID.rstrip('0123456789').rstrip('.')):
			is_gao_aTIS_translated = 'yes'
		elif is_gao_translated == 'yes':
			is_gao_aTIS_translated = 'no'
		
		# RNA Secondary Structure
		wind_len = RNAseq_window_length
		mRNA_len = len(mRNA_transcript)

		uORF_DNA_end = uORF_DNA_start + uORF_length
		uORF_start_dGs = get_RNA_dGs(mRNA_transcript, mRNA_len, wind_len, uORF_DNA_start, True)
		uORF_end_dGs = get_RNA_dGs(mRNA_transcript, mRNA_len, wind_len, uORF_DNA_end, False)
		CDS_start_dGs = get_RNA_dGs(mRNA_transcript, mRNA_len, wind_len, CDS_DNA_start, True)
		
		# Evolutionary Conservation
		if not GERP_score_dict.get(uORF_ID):
			print "ERROR! uORF_ID not found in GERP score dictionary"
			print "Exiting..."
			sys.exit(1)
		GERP_start_score, GERP_stop_score = GERP_score_dict.get(uORF_ID)
		if not GERP_elements_dict.get(uORF_ID):
			print "ERROR! uORF_ID not found in GERP elements dictionary"
			print "Exiting..."
			sys.exit(1)
		GERP_elements_perc = GERP_elements_dict.get(uORF_ID)
		PhyloCSF_score = PhyloCSF_dict.get(uORF_ID)
				
		# Rare Codon analysis
		rare_codon_count = {}
		for codon in Rarest_Codons:
			rare_codon_count[codon] = 0
		j = 0
		while j < uORF_length:
			current_codon = uORF_DNA_seq[j:j+3]
			if current_codon in Rarest_Codons:
				rare_codon_count[current_codon] += 1
			j += 3
		rare_codon_result = []
		for codon in Rarest_Codons:
			rare_codon_result.append(str(rare_codon_count[codon]))
		
		# + PhyloCSF

		# GTEX_information:
		if not GTEX_dict.get(transcript_ID):
			GTEX_tissues_data = [0] * 32
		else:
			GTEX_tissues_data = GTEX_dict.get(transcript_ID)
			GTEX_length = len(GTEX_tissues_data)
		Thyroid_ex = GTEX_tissues_data[0]
		Testis_ex = GTEX_tissues_data[1]
		Cervix_Uteri_ex = GTEX_tissues_data[2]
		Adipose_Tissue_ex = GTEX_tissues_data[3]
		Breast_ex = GTEX_tissues_data[4]
		Vagina_ex = GTEX_tissues_data[5]
		Ovary_ex = GTEX_tissues_data[6]
		Stomach_ex = GTEX_tissues_data[7]
		Fallopian_Tube_ex = GTEX_tissues_data[8]
		Bone_Marrow_ex = GTEX_tissues_data[9]
		Spleen_ex = GTEX_tissues_data[10]
		Bladder_ex = GTEX_tissues_data[11]
		Blood_ex = GTEX_tissues_data[12]
		Colon_ex = GTEX_tissues_data[13]
		Prostate_ex = GTEX_tissues_data[14]
		Pancreas_ex = GTEX_tissues_data[15]
		Blood_Vessel_ex = GTEX_tissues_data[16]
		Liver_ex = GTEX_tissues_data[17]
		Heart_ex = GTEX_tissues_data[18]
		Small_Intestine_ex = GTEX_tissues_data[19]
		Uterus_ex = GTEX_tissues_data[20]
		Pituitary_ex = GTEX_tissues_data[21]
		Muscle_ex = GTEX_tissues_data[22]
		Nerve_ex = GTEX_tissues_data[23]
		Adrenal_Gland_ex = GTEX_tissues_data[24]
		Brain_ex = GTEX_tissues_data[25]
		Salivary_Gland_ex = GTEX_tissues_data[26]
		Lung_ex = GTEX_tissues_data[27]
		Skin_ex = GTEX_tissues_data[28]
		Esophagus_ex = GTEX_tissues_data[29]
		Kidney_ex = GTEX_tissues_data[30]
		Tissue_Entropy_ex = GTEX_tissues_data[31]

		# SNP information:
		SNP_data = SNP_dict.get(uORF_ID)
		SNP_number = float(SNP_data[0])
		SNPs_per_length = '{0:.14f}'.format(SNP_number / uORF_length)
		heterozygosity = float(SNP_data[1])
		heterozygosity_per_length = '{0:.14f}'.format(heterozygosity / uORF_length)

		# print the results
		print '\t'.join([uORF_ID, str(cap_to_uORF_distance), str(uORF_length), str(CDSstart_to_uORFstart), 
			  str(CDSstart_to_uORFend), str(percent_A), str(percent_T), str(percent_G), str(percent_C),
			  str(uORF_context), str(UTR_size), str(ATG_uORFs_in_transcript), str(CTG_uORFs_in_transcript), str(TTG_uORFs_in_transcript), str(GTG_uORFs_in_transcript), str(ATC_uORFs_in_transcript), str(ATT_uORFs_in_transcript), str(ATA_uORFs_in_transcript), str(ACG_uORFs_in_transcript), str(AAG_uORFs_in_transcript), str(AGG_uORFs_in_transcript), str(same_codon_uORFs_in_transcript), str(total_start_codons), str(codon_frequency)] + \
			  rare_codon_result + uORF_start_dGs + uORF_end_dGs + CDS_start_dGs + \
			  [str(GERP_start_score), str(GERP_stop_score), str(GERP_elements_perc),
			  str(Thyroid_ex), str(Testis_ex), str(Cervix_Uteri_ex), str(Adipose_Tissue_ex), 
			  str(Breast_ex), str(Vagina_ex), str(Ovary_ex), str(Stomach_ex), 
			  str(Fallopian_Tube_ex), str(Bone_Marrow_ex), str(Spleen_ex), str(Bladder_ex), 
			  str(Blood_ex), str(Colon_ex), str(Prostate_ex), str(Pancreas_ex), 
			  str(Blood_Vessel_ex), str(Liver_ex), str(Heart_ex), str(Small_Intestine_ex), 
			  str(Uterus_ex), str(Pituitary_ex), str(Muscle_ex), str(Nerve_ex), str(Adrenal_Gland_ex), 
			  str(Brain_ex), str(Salivary_Gland_ex), str(Lung_ex), str(Skin_ex), 
			  str(Esophagus_ex), str(Kidney_ex), str(Tissue_Entropy_ex),
			  str(SNPs_per_length), str(heterozygosity_per_length), str(SNP_number), 
			  str(heterozygosity), str(uORF_Noderer_score),
			  str(is_fritsch_translated), str(is_lee_translated), str(is_gao_translated),
			  str(is_fritsch_aTIS_translated), str(is_lee_aTIS_translated), str(is_gao_aTIS_translated)])

	else:
		print 'Error, 5\'UTR or CDS sequence not found for transcript ID ' + transcript_ID + '!'	
		print 'Exiting...'
		sys.exit(1)

def stats(uORFs_seq, uORFs_ATG_seq, uORFs_CTG_seq, uORFs_TTG_seq, uORFs_GTG_seq, uORFs_ATC_seq, uORFs_ATT_seq, uORFs_ATA_seq, uORFs_ACG_seq, uORFs_AAG_seq, uORFs_AGG_seq, fiveprime_seq, cds_seq, fritsch_uORF_IDs, lee_uORF_IDs, gao_uORF_IDs, 
		  lee_aTISs_IDs, fritsch_aTISs_IDs, gao_aTISs_IDs, GERP_start_stop, GERP_elements, PhyloCSF, SNPs_uORFs_intersect, GTEX_tissue_data, Noderer_data, partition_number):

	# make sequence dictionaries
#	print("making uORFs_ID_dict...")
	uORFs_ID_dict = makedict(uORFs_seq, 2, True, 1, True)
#	print("making uORFs_*CODON*_ID_dict...")
	uORFs_ATG_ID_dict = makedict(uORFs_ATG_seq, 2, True, 1, True)
	uORFs_CTG_ID_dict = makedict(uORFs_CTG_seq, 2, True, 1, True)
	uORFs_TTG_ID_dict = makedict(uORFs_TTG_seq, 2, True, 1, True)
	uORFs_GTG_ID_dict = makedict(uORFs_GTG_seq, 2, True, 1, True)
	uORFs_ATA_ID_dict = makedict(uORFs_ATA_seq, 2, True, 1, True)
	uORFs_ATC_ID_dict = makedict(uORFs_ATC_seq, 2, True, 1, True)
	uORFs_ATT_ID_dict = makedict(uORFs_ATT_seq, 2, True, 1, True)
	uORFs_ACG_ID_dict = makedict(uORFs_ACG_seq, 2, True, 1, True)
	uORFs_AAG_ID_dict = makedict(uORFs_AAG_seq, 2, True, 1, True)
	uORFs_AGG_ID_dict = makedict(uORFs_AGG_seq, 2, True, 1, True)
#	print("making fiveprime_seq_dict...")
	fiveprime_seq_dict = makedict(fiveprime_seq, 1, True)
#	print("making CDS_seq_dict...")
	CDS_seq_dict = makedict(cds_seq, 1, True)
	
	# Populate dictionaries for translated uORFs and aTISs
#	print("making fritsch_uORF_dict...")
	fritsch_uORF_dict = make_profile_dict(fritsch_uORF_IDs, 0, 'yes')
#	print("making lee_uORF_dict...")
	lee_uORF_dict = make_profile_dict(lee_uORF_IDs, 0, 'yes')
#	print("making gao_uORF_dict...")
	gao_uORF_dict = make_profile_dict(gao_uORF_IDs, 0, 'yes')
#	print("making lee_aTIS_dict...")
	lee_aTIS_dict = make_profile_dict(lee_aTISs_IDs, 0, True)
#	print("making fritsch_aTIS_dict...")
	fritsch_aTIS_dict = make_profile_dict(fritsch_aTISs_IDs, 0, True)
#	print("making gao_aTIS_dict...")
	gao_aTIS_dict = make_profile_dict(gao_aTISs_IDs, 0, True)
	
#	print("making GTEX_dict...")
	GTEX_dict = make_GTEX_dict(GTEX_tissue_data, 0)

#	print("making SNP_dict...")
	SNP_dict = make_SNP_dict(SNPs_uORFs_intersect, 0)

#	print("making Noderer_Score_dict...")
	Noderer_dict = make_Noderer_Score_dict(Noderer_data)
		
	"""fritsch_uORF_dict = {}
	for line in fritsch_uORF_IDs:
		fritsch_uORF_dict[line.split()[0]] = 'yes'S
	lee_uORF_dict = {}
	for line in lee_uORF_IDs:
		lee_uORF_dict[line.split()[0]] = 'yes'
	lee_aTIS_dict = {}
	for line in lee_aTISs_IDs:
		lee_aTIS_dict[line.split()[0]] = True
	fritsch_aTIS_dict = {}
	for line in fritsch_aTISs_IDs:
		fritsch_aTIS_dict[line.split()[0]] = True"""
	
	# Make Conservation Information Dictionaries
#	print("making GERP_score_dict...")
	GERP_score_dict = make_cons_dict(GERP_start_stop, 0, 2, 1)
#	print("making GERP_elements_dict...")
	GERP_elements_dict = make_cons_dict(GERP_elements, 0, 1, 1)
#	print("making PhyloCSF_dict...")
	PhyloCSF_dict = make_cons_dict(PhyloCSF, 0, 1, 0)

	"""
	{}
	GERP_start_stop.readline()					# skip intro line
	for line in GERP_start_stop:
		fields = line.split()
		uORF_ID, start_score, stop_score = fields[0], fields[1], fields[2]
		GERP_score_dict[uORF_ID] = (start_score, stop_score)
	
	GERP_elements_dict = {}
	GERP_elements.readline()					# skip intro line
	for line in GERP_elements:
		fields = line.split()
		uORF_ID, base_percent = fields[0], fields[1]
		GERP_elements_dict[uORF_ID] = base_percent
	
	PhyloCSF_dict = {}
	for line in PhyloCSF:
		fields = line.split()
		uORF_ID, cons_score = fields[0], fields[1] """
		
	
	rarest_codons_fields = ['Number-codon-' + codon for codon in Rarest_Codons]
	
	# Automatically generate field names for RNA secondary structure deltaG values
	uORF_start_fields = ['mRNA_dG_uORF_start_' + str(i) + '-' + str(i+RNAseq_window_length-1) 
						 for i in secstruct_start_positions]
	uORF_end_fields = ['mRNA_dG_uORF_end_' + str(i) + '-' + str(i+RNAseq_window_length-1) 
						 for i in secstruct_uORF_end_positions]
	CDS_start_fields = ['mRNA_dG_CDS_start_' + str(i) + '-' + str(i+RNAseq_window_length-1) 
						 for i in secstruct_start_positions]	
	
	# Table header
	print '#uORF-transcriptID\tcap-to-uORF-distance\tuORF-length\t'		   + \
		  'CDS-start-to-uORF-start\tCDS-start-to-uORF-end\t%A\t%T\t%G\t%C\t'  + \
		  'uORF-context\t5\'UTR-size\tuORF-number-ATG\tuORF-number-CTG\tuORF-number-TTG\tuORF-number-GTG\tuORF-number-ATC\tuORF-number-ATT\tuORF-number-ATA\tuORF-number-ACG\tuORF-number-AAG\tuORF-number-AGG\tuORF-number-same-codon\tuORF-total-start-codon\tcodon-frequency\t'  + '\t'.join(rarest_codons_fields) + '\t' + \
		  '\t'.join(uORF_start_fields + uORF_end_fields + CDS_start_fields) + '\t' + \
		  'GERP-startcodon-score\tGERP-stopcodon-score\tGERP-elements-percent\t' + \
		'Thyroid\tTestis\tCervix Uteri\tAdipose Tissue\tBreast\tVagina\tOvary\tStomach\t' + \
		'Fallopian Tube\tBone Marrow\tSpleen\tBladder\tBlood\tColon\tProstate\tPancreas\t' + \
		'Blood Vessel\tLiver\tHeart\tSmall Intestine\tUterus\tPituitary\tMuscle\t' + \
		'Nerve\tAdrenal Gland\tBrain\tSalivary Gland\tLung\tSkin\tEsophagus\tKidney\tTissue Entropy\t' + \
		'SNPs/length\tHeterozygosity/length\tSNPs\tHeterozygosity\tNoderer Context Score\t' + \
		'Fritsch-translated\tLee-translated\tGao-translated\tFritsch-aTIS-translated\tLee-aTIS-translated\tGao-aTIS-translated'
	
	# Output information for each genomic uORF
	iteration = int(partition_number)
	uORFs_seq.seek(0)
	all_data = uORFs_seq.readlines()		# skip intro line
	amount = (num_lines-1)/16
	amount_new = amount*2
	remaining_amount = (num_lines-1)%16
	i=1+((iteration-1)*amount_new)
	upper_limit = amount_new*iteration
	if iteration == 8:
		upper_limit = amount_new*iteration + remaining_amount
	data = all_data[i]
	while data:
		while (i <= upper_limit):
			fields = data.split('|')
			uORF_ID = fields[2]
			seq_line = all_data[i+1]
			#uORF_DNA_seq is the sequence of the uORF, as pulled from the fasta file.
			uORF_DNA_seq = seq_line.split()[0]
			statoutput(uORF_ID, uORF_DNA_seq, uORFs_ID_dict, uORFs_ATG_ID_dict, uORFs_CTG_ID_dict, uORFs_TTG_ID_dict, uORFs_GTG_ID_dict, uORFs_ATC_ID_dict, uORFs_ATT_ID_dict, uORFs_ATA_ID_dict, uORFs_ACG_ID_dict, uORFs_AAG_ID_dict, uORFs_AGG_ID_dict, fiveprime_seq_dict, 
					   CDS_seq_dict, GERP_score_dict, GERP_elements_dict, PhyloCSF_dict, GTEX_dict, SNP_dict, Noderer_dict, fritsch_uORF_dict, 
					   lee_uORF_dict, gao_uORF_dict, lee_aTIS_dict, fritsch_aTIS_dict, gao_aTIS_dict)
			data = all_data[i+2]
			if data and data.startswith('#'):
				data = all_data[i+3]
				i += 4
			i += 2
		break

# COMMENT HERE!
def get_rare_codons(ncodons):
	sorted_frequency_table = sorted(Codon_Frequency, key=lambda codon: int(codon[2]))[3:]		# slice out the three stop codons
	return [item[0] for item in sorted_frequency_table][:ncodons]
	
if __name__ == '__main__':
	if len(sys.argv) != 27:
		print "Incorrect number of arguments"
		print(len(sys.argv))
		sys.exit(1)
	uORFs_codons = ['ACC', 'ATG', 'AAG', 'ACG', 'ATC', 'AAC', 'ATA', 'AGG', 
					'ACT', 'AGC', 'ACA', 'AGA', 'ATT', 'CTG', 'CTA', 'CTC', 
					'AAA', 'CCG', 'CAG', 'CTT', 'GGT', 'CGA', 'CCA', 'TCT', 
					'CGG', 'TTT', 'GGG', 'TAG', 'GGA', 'GGC', 'TAC', 'TTC', 
					'TCG', 'TTG', 'TCC', 'GAA', 'TCA', 'GCA', 'GTA', 'GCC', 
					'GTC', 'GCG', 'GTG', 'GAG', 'GTT', 'GCT', 'GAC', 'CGT', 
					'CGC']
	
	# From GenScript: [Codon, frequency per thousand, number]			
	Codon_Frequency = [['TTT', '16.9', '336562'], ['TCT', '14.6', '291040'], ['TAT', '12.0', '239268'], 
					   ['TGT', '9.9', '197293'], ['TTC', '20.4', '406571'], ['TCC', '17.4', '346943'], 
					   ['TAC', '15.6', '310695'], ['TGC', '12.2', '243685'], ['TTA', '7.2', '143715'], 
					   ['TCA', '11.7', '233110'], ['TAA', '0.7', '14322'], ['TGA', '1.3', '25383'], 
					   ['TTG', '12.6', '249879'], ['TCG', '4.5', '89429'], ['TAG', '0.5', '10915'], 
					   ['TGG', '12.8', '255512'], ['CTT', '12.8', '253795'], ['CCT', '17.3', '343793'], 
					   ['CAT', '10.4', '207826'], ['CGT', '4.7', '93458'], ['CTC', '19.4', '386182'], 
					   ['CCC', '20.0', '397790'], ['CAC', '14.9', '297048'], ['CGC', '10.9', '217130'], 
					   ['CTA', '6.9', '138154'], ['CCA', '16.7', '331944'], ['CAA', '11.8', '234785'], 
					   ['CGA', '6.3', '126113'], ['CTG', '40.3', '800774'], ['CCG', '7.0', '139414'], 
					   ['CAG', '34.6', '688316'], ['CGG', '11.9', '235938'], ['ATT', '15.7', '313225'], 
					   ['ACT', '12.8', '255582'], ['AAT', '16.7', '331714'], ['AGT', '11.9', '237404'], 
					   ['ATC', '21.4', '426570'], ['ACC', '19.2', '382050'], ['AAC', '19.5', '387148'], 
					   ['AGC', '19.4', '385113'], ['ATA', '7.1', '140652'], ['ACA', '14.8', '294223'], 
					   ['AAA', '24.0', '476554'], ['AGA', '11.5', '228151'], ['ATG', '22.3', '443795'], 
					   ['ACG', '6.2', '123533'], ['AAG', '32.9', '654280'], ['AGG', '11.4', '227281'], 
					   ['GTT', '10.9', '216818'], ['GCT', '18.6', '370873'], ['GAT', '22.3', '443369'], 
					   ['GGT', '10.8', '215544'], ['GTC', '14.6', '290874'], ['GCC', '28.5', '567930'], 
					   ['GAC', '26.0', '517579'], ['GGC', '22.8', '453917'], ['GTA', '7.0', '139156'], 
					   ['GCA', '16.0', '317338'], ['GAA', '29.0', '577846'], ['GGA', '16.3', '325243'], 
					   ['GTG', '28.9', '575438'], ['GCG', '7.6', '150708'], ['GAG', '40.8', '810842'], 
					   ['GGG', '16.4', '326879']]
	
	rare_codon_number = 6
	Rarest_Codons = get_rare_codons(rare_codon_number)
	
	# Define global variables
	RNAseq_window_length = 40
	secstruct_start_positions = [-20, 0, 20, 40, 60, 80, 100]
	secstruct_uORF_end_positions = [-40, -20, 0, 20, 40]
   
	i = 1
	uORFs_seq_file = sys.argv[i]
	fiveprime_seq_file = sys.argv[i+1]
	cds_seq_file = sys.argv[i+2]
	fritsch_uORF_IDs_file = sys.argv[i+3]
	lee_uORF_IDs_file = sys.argv[i+4]
	gao_uORF_IDs_file = sys.argv[i+5]
	lee_aTISs_IDs_file = sys.argv[i+6]
	fritsch_aTISs_IDs_file = sys.argv[i+7]
	gao_aTISs_IDs_file = sys.argv[i+8]
	GERP_start_stop_file = sys.argv[i+9]
	GERP_elements_file = sys.argv[i+10]
	PhyloCSF_file = sys.argv[i+11]
	uORFs_ATG_seq_file = sys.argv[i+12]
	uORFs_CTG_seq_file = sys.argv[i+13]
	uORFs_TTG_seq_file = sys.argv[i+14]
	uORFs_GTG_seq_file = sys.argv[i+15]
	uORFs_ATC_seq_file = sys.argv[i+16]
	uORFs_ATT_seq_file = sys.argv[i+17]
	uORFs_ATA_seq_file = sys.argv[i+18]
	uORFs_ACG_seq_file = sys.argv[i+19]
	uORFs_AAG_seq_file = sys.argv[i+20]
	uORFs_AGG_seq_file = sys.argv[i+21]
	SNPs_uORFs_intersect_file = sys.argv[i+22]
	GTEX_tissue_data_file = sys.argv[i+23]
	Noderer_data_file = sys.argv[i+24]
	partition_number = sys.argv[i+25]

	uORFs_seq_linecount = open(uORFs_seq_file)
	num_lines = sum(1 for line in uORFs_seq_linecount)
	uORFs_seq_linecount.close()

	uORFs_seq = open(uORFs_seq_file, 'r')
	fiveprime_seq = open(fiveprime_seq_file, 'r')
	cds_seq = open(cds_seq_file, 'r')
	fritsch_uORF_IDs = open(fritsch_uORF_IDs_file, 'r')
	lee_uORF_IDs = open(lee_uORF_IDs_file, 'r')
	gao_uORF_IDs = open(gao_uORF_IDs_file, 'r')
	lee_aTISs_IDs = open(lee_aTISs_IDs_file, 'r')
	fritsch_aTISs_IDs = open(fritsch_aTISs_IDs_file, 'r')
	gao_aTISs_IDs = open(gao_aTISs_IDs_file, 'r')
	GERP_start_stop = open(GERP_start_stop_file, 'r')
	GERP_elements = open(GERP_elements_file, 'r')
	PhyloCSF = open(PhyloCSF_file, 'r')
	uORFs_ATG_seq = open(uORFs_ATG_seq_file, 'r') 
	uORFs_CTG_seq = open(uORFs_CTG_seq_file, 'r')
	uORFs_TTG_seq = open(uORFs_TTG_seq_file, 'r')
	uORFs_GTG_seq = open(uORFs_GTG_seq_file, 'r')
	uORFs_ATC_seq = open(uORFs_ATC_seq_file, 'r')
	uORFs_ATT_seq = open(uORFs_ATT_seq_file, 'r')
	uORFs_ATA_seq = open(uORFs_ATA_seq_file, 'r')
	uORFs_ACG_seq = open(uORFs_ACG_seq_file, 'r')
	uORFs_AAG_seq = open(uORFs_AAG_seq_file, 'r')
	uORFs_AGG_seq = open(uORFs_AGG_seq_file, 'r')
	SNPs_uORFs_intersect = open(SNPs_uORFs_intersect_file, 'r')
	GTEX_tissue_data = open(GTEX_tissue_data_file, 'r')
	Noderer_data = open(Noderer_data_file, 'r')

	stats(uORFs_seq, uORFs_ATG_seq, uORFs_CTG_seq, uORFs_TTG_seq, uORFs_GTG_seq, uORFs_ATC_seq, uORFs_ATT_seq, uORFs_ATA_seq, uORFs_ACG_seq, uORFs_AAG_seq, uORFs_AGG_seq, fiveprime_seq, cds_seq, fritsch_uORF_IDs, 
	lee_uORF_IDs, gao_uORF_IDs, lee_aTISs_IDs, fritsch_aTISs_IDs, gao_aTISs_IDs, 
	GERP_start_stop, GERP_elements, PhyloCSF, SNPs_uORFs_intersect, GTEX_tissue_data, Noderer_data, partition_number)
