## Author: Patrick McGillivray

## The following README file documents the procedure for generating a genome-wide set
## of predicted positive (translated) uORFs, as described in the manuscript “A
## comprehensive catalog of predicted functional upstream open reading frames” by Patrick ## McGillivray, Russell Ault, Mayur Pawashe, Robert Kitchen, Suganthi Balasubramanian,
## and Mark Gerstein

## Acquire translated uORF locations from source data:

Lee, S., et al. "Global mapping of translation initiation sites in mammalian cells at single-nucleotide resolution." Proceedings of the National Academy of Sciences 109.37 (2012): E2424-E2432.

1. Download Supplemental Data 1: TIS positions identified in HEK293 cells
2. grep command for 5’UTR positions > Lee2012_transcriptIDs

Gao, X., et al. "Quantitative profiling of initiating ribosomes in vivo." Nature methods 12.2 (2015): 147-153.

1. Download Supplementary Table 1: HEK293 QTI-seq data sets
2. grep command for upstream to start site positions > Gao2014_transcriptIDs

Fritsch, C., et al. "Genome-wide search for novel human uORFs and N-terminal protein extensions using ribosomal footprinting." Genome research 22.11 (2012): 2208-2218.

1. Download hub_7251_puroTranscript_fromUCSCbrowser from http://gengastro.1med.uni-kiel.de/suppl/footprint/ 
2. grep command for upstream to start site positions > puroTranscriptIDs.txt 

## Converting from RefSeq to Ensembl format using Ensembl BioMart (GRCH37):

The files Lee2012_transcriptIDs (for Lee et al. data), puroTranscriptIDs.txt (for the Fritsch et al data), and Gao2014_transcriptIDs (for the Gao data), contain the transcript data for these three experiments in RefSeq format. These must be converted to Ensembl format .tsv files, using BioMart.

IDs can be uploaded under “FILTER” as RefSeq mRNA IDs.
The output must include the following attributes (selected from the biomart browser under “ATTRIBUTES” – under “GENE”, with the exception of RefSeq mRNA, that is under “EXTERNAL”): Ensembl Gene ID; Ensembl Transcript ID; Chromosome Name; Transcript Start (bp); Transcript End (bp); RefSeq mRNA [e.g. NM_001195597]

Unique output should be selected, as TSV.

Output should be retitled as RefSeq_ENSEMBL_Lee_match.txt, RefSeq_ENSEMBL_match.txt, and RefSeq_ENSEMBL_Gao_match.txt, for Lee, Fritsch, and Gao experiments respectively.

## Update the aTIS dictionaries:

Extract all of the lines from classified purotranscripts (Fritsch) that contain the string ‘aTIS’, using a grep command, or a similar tool. Convert the RefSeq_IDs for these transcripts to ENSEMBL ID using biomart (GRCH37). Put these IDs in the following directories:

Lee:

~/ribosome_profiling/globalmapping_uORFs_paper_Leeetal2012/TIS_lee_aTIS_EnsemblIDs

Fritsch: 

~/ribosome_profiling/novelhuman_uORFs_paper_Fritschetal2012/UCSCBrowser_identifiedTISs/TIS_fritsch_aTIS_EnsemblIDs

Gao:

~/ribosome_profiling/uORFs_paper_Gaoetal2014/TIS_gao_aTIS_EnsemblIDs

*NOTE: none of the above commands are specific for a certain start codon. They are general commands. This code does not need to be modified, according to start codon.

## Create uORF_ATG.fa + uORFs_ATG.gtf files:

Using hg19:

python ~/python_scripts/uorfs.py /gencode.v19.annotation.gtf /$directory_containing_chromosomal_fasta_files ~/uORFs/ -ATG

The following code, runs a script, that loops through all possible start codons:

python ~/python_scripts/all_start_codons_uORFs.py

Using GRCh38:

python ~/python_scripts/uorfs.py ~/gencode.v21.annotation.gtf ~/$directory_containing_chromosomal_fasta_files ~/uORFs/ -ATG

*Note that the fasta files must correspond to the GENCODE version used.

Get the fasta files from the following:

ftp hgdownload.cse.ucsc.edu
username: anonymous
password: $user_email

‘prompt’ toggles off interactive mode.

cd /goldenPath/hg38/chromosomes/ --- for Gencode v21
cd /goldenPath/hg19/chromosomes/ --- for Gencode v19

mget -a

gunzip *.gz

## Create the uORF_ATG.bed file:

python ~/python_scripts/uorf_gtf_to_bed.py ATG

Or for all possible start codons by the following script:

python ~/python_scripts/all_start_codons_uORFs.py

## Create the CDS.fa file:

python ~/python_scripts/cds.py ~/gencode.v19.annotation.gtf ~/$directory_containing_chromosomal_fasta_files ~/CDS/

Note: This file has no dependence, on the start codon, it simply retrieves CDS sequences.

## Create the five_primes.fa file:

python ~/python_scripts/five_primes.py ~/gencode.v19.annotation.gtf /$directory_containing_chromosomal_fasta_files ~/five_primes/

This file has no dependence on the start codon, it simply retrieves 5’ sequences.

### Extract SNP data, from the 1000 Genomes project (GENCODE release version is important here):

For each of the possible start codons, run the following script:

python /net/gerstein/pdm32/python_scripts/SNPdataextract.py > heteroSNPfile.txt

Note that the alignments are based on GRCh37. The fasta file used can be found here:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

## Creating the tissue expression/entropy file:

python ~/GTEXexpressiondatalooping.py > GTEX_tissue_meanexpression_entropy.tsv

## Create the uORF stats table:

python all_start_codons_uORFs_partition.py
bash partitioned_uORF_stats_most_frequent

python ~/python_scripts/join_uORFs_partition.py ATG
python ~/python_scripts/join_uORFs_partition.py CTG
python ~/python_scripts/join_uORFs_partition.py GTG
python ~/python_scripts/join_uORFs_partition.py TTG
python ~/python_scripts/join_uORFs_partition.py ACG
python ~/python_scripts/join_uORFs_partition.py AAG
python ~/python_scripts/join_uORFs_partition.py AGG
python ~/python_scripts/join_uORFs_partition.py ATC
python ~/python_scripts/join_uORFs_partition.py ATA
python ~/python_scripts/join_uORFs_partition.py ATT

## Run analysis until discretization:

bash uORF_table_to_end_analysis_allcodons
