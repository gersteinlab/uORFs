## Noderer context file to be retrieved from supplement of DOI: 10.15252/msb.20145136

#! /bin/bash

# Check if there are 3 arguments
if [[ $4 || ! $3 ]]
then
	exit 0
fi

CODON=$1
FRITSCHUORFIDS=$2
LEEUORFIDS=$3
HOME="~/"
PHYLOCSF="${HOME}/PhyloCSF/uORFs_ATG.csf"
ATGSEQ="${HOME}/uORFs/uORFs_ATG.fa"
CTGSEQ="${HOME}/uORFs/uORFs_CTG.fa"
TTGSEQ="${HOME}/uORFs/uORFs_TTG.fa"
GTGSEQ="${HOME}/uORFs/uORFs_GTG.fa"
SNPS="${HOME}/SNPs_intersections/${CODON}_heteroSNPfile.txt"
TISSUE="${HOME}/python_scripts/GTEX_tissue_meanexpression_entropy.tsv"
CONTEXT="${HOME}/Noderer_context/msb145136-sup-0003-TableS2.txt"
TABLE="${HOME}/uORF_stats/uorfs_stats_table_${CODON}.txt"
ERROR="${HOME}/uORF_stats/uorfs_stats_table_${CODON}.err"
LEETISIDS="${HOME}/ribosome_profiling/globalmapping_uORFs_paper_Leeetal2012/TIS_lee_aTIS_EnsemblIDs"
FRITISIDS="${HOME}/ribosome_profiling/novelhuman_uORFs_paper_Fritschetal2012/UCSCBrowser_identifiedTISs/TIS_fritsch_aTIS_EnsemblIDs"

GERPSTART="${HOME}/GERP_results/GERPrate_start_stop_avg_${CODON}"
GERPELEMENT="${HOME}/GERP_results/uORF_${CODON}_GERPelement_filtered_overlap"

cd ${HOME}/python_scripts/

python uorf_gtf_to_GERP.py ${CODON} ~/uORFs/uORFs_${CODON}.gtf

python uorf_gtf_to_GERP.py ${CODON} ${HOME}/uORFs/uORFs_${CODON}.gtf
python uorfs_GERPrate.py ${CODON} ~/nmd/data/bases ${HOME}/GERP_results/GERP_cache_files -v ${HOME}/uORFs/uORFs_${CODON}.gerp
 python uorfs_GERPratetoCodonAverage.py ${HOME}/GERP_results/GERPrate_start_stop_${CODON} > $GERPSTART #Check.
 bash ${HOME}/bash_scripts/GERPelements.bash ${CODON} #Check.

 cd ~/ribosome_profiling/novelhuman_uORFs_paper_Fritschetal2012/UCSCBrowser_identifiedTISs

 python ~/python_scripts/intersectTIS.py 0 ~/uORFs/uORFs_${CODON}.gtf hub_7251_puroClassifiedTIS_fromUCSCbrowser RefSeq_Ensembl_match.txt -m > fritsch_intersect_${CODON}_uORFs.txt;
 python ~/python_scripts/filestrip.py -fri fritsch_intersect_${CODON}_uORFs.txt -${CODON} > fritsch_${CODON}_uORF_IDs;

 cd ~/ribosome_profiling/globalmapping_uORFs_paper_Leeetal2012

 python ~/python_scripts/lee_to_BED.py ~/gencode.v19.annotation.gtf ~/nmd/data/genome TIS_lee_HEK293 RefSeq_Ensembl_Lee_match.txt;
 python ~/python_scripts/intersectTIS.py 0 ~/uORFs/uORFs_${CODON}.gtf lee_uORF2.bed RefSeq_Ensembl_Lee_match.txt -m > lee_intersect_${CODON}_uORFs.txt ;
 python ~/python_scripts/filestrip.py -fri lee_intersect_${CODON}_uORFs.txt -${CODON} > lee_${CODON}_uORF_IDs;

cd ~/python_scripts
s
cd ${HOME}/
python python_scripts/uorfstats_experiment_allcodons_partition.py uORFs/uORFs_${CODON}.fa five_primes/five_primes.fa CDS/CDSs.fa $FRITSCHUORFIDS $LEEUORFIDS $LEETISIDS $FRITISIDS $GERPSTART $GERPELEMENT $PHYLOCSF $ATGSEQ $CTGSEQ $TTGSEQ $GTGSEQ $SNPS $TISSUE $CONTEXT > $TABLE 2> $ERROR
