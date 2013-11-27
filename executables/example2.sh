#!/bin/csh

set MCRDIR=/data/suzaku/yau/MCR/r2012a/v717/

set outdir=/data/suzaku/yau/oncoseq/temp/
set hgtable=/data/suzaku/yau/oncoseq/code/v1.02/config/hgTables_b37.txt
set tumourstatestable=/data/suzaku/yau/oncoseq/code/v1.02/config/tumourStates.txt

set seqtype=cg
set samplename=cellmix-1
set infile=/data/suzaku/yau/Enric/cg/somatic/masterVarBeta-GS000017809-ASM.tsv

./run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:22]" --n_train 30000 --maxploidy 4.5 --minploidy 1.5 --normalcontamination --tumourstatestable $tumourstatestable --maxnormalcontamination 0.6 --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $infile --outdir $outdir
