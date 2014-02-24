infile=/data/cyau/Ludwig/data/pileup/COLO320.txt.gz
outfile=/data/cyau/temp/COLO320_oncoseq.txt

infile=/data/cyau/wgs500/cll/AS_CLL_003D.txt
outfile=/data/cyau/temp/COLO320_oncoseq.txt

infile=/data/cyau/oncoseq/FH_BLC_4070T.txt
outfile=/data/cyau/temp/FH_BLC_4070T.txt

snpfile=/data/cyau/oncoseq/illumina_hg19.txt
snpfile=/data/cyau/refseq/snps.bed

perl -T process_pileup.pl --infile=$infile --outfile=$outfile --snpfile=$snpfile --minreads=20 --minaltread=2


#perl -T process_pileup.pl --infile=/data/suzaku/yau/Ludwig/data/pileup/COLO320.txt --outfile=/data/suzaku/yau/temp/COLO320_oncoseq.txt --snpfile=/data/suzaku/yau/oncoseq/test-folder/illumina_hg19.txt 

#exit;

#infile=/data/suzaku/yau/cg/data/cancer/HCC1187-H-200-37-ASM/masterVarBeta-HCC1187-H-200-37-ASM-T1-N1.tsv
#outfile=/data/suzaku/yau/temp/masterVarBeta-HCC1187-H-200-37-ASM-T1-N1.tsv
#snpfile=illumina_hg19.txt

#perl process_mastervar.pl --infile=$infile --outfile=$outfile --snpfile=$snpfile 

#infile=/data/suzaku/yau/coombes/mastervar/masterVarBeta-GS000008109-ASM-T1-N1.tsv.bz2
#infile=/data/suzaku/yau/cg/data/cancer/HCC1187-H-200-37-ASM/masterVarBeta-HCC1187-H-200-37-ASM-T1-N1.tsv
#outfile=/data/suzaku/yau/temp/test
#snpfile=illumina_hg19.txt

#perl process_mastervar.pl --infile=$infile --outfile=$outfile

#infile=/data/suzaku/yau/Enric/cg/masterVarBeta/masterVarBeta-GS000017810-ASM.tsv.bz2
#outfile=/data/suzaku/yau/temp/cgtest.txt
#normalfile=/data/suzaku/yau/Enric/cg/masterVarBeta/masterVarBeta-GS000017809-ASM.tsv.bz2
#perl process_mastervar_standard_versus_normal.pl --infile=$infile --outfile=$outfile --normalfile=$normalfile
