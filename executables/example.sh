#!/bin/csh

set MCRDIR=/data/cyau/MCR/r2013b/v82/

set outdir=/data/cyau/temp/
set hgtable=/home/cyau/projects/oncosnpseq/config/hgTables_b37.txt
set tumourstatestable=/home/cyau/projects/oncosnpseq/config/tumourStates.txt

set seqtype=cg
set samplename=snps-100
set infile=/home/cyau/projects/oncosnpseq/example_data/snps-100.txt

echo $MCRDIR 

./run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:22]" --n_train 30000 --maxploidy 4.5 --minploidy 1.5 --normalcontamination --tumourstatestable $tumourstatestable  --maxnormalcontamination 0.6 --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $infile --outdir $outdir
