#!/bin/csh

set MCRDIR=/data/cyau/MCR/r2013b/v82/

set hgtable=../config/hgTables_b37.txt
set tumourstatestable=../config/tumourStatesOncoSNPSEQ.txt

set seqtype=cg
set indir=/data/cyau/Enric/output/oncosnpseq_input/
set outdir=/data/cyau/Enric/output/oncosnpseq_output/


foreach file(`find $indir -type f -name "*"`)
	echo $file
	set samplename=`basename $file`
	echo $samplename
	rm $outdir/$samplename.out
	nohup ../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:23]" --lambda1 "[ 1000 500 100 50 30 10]" --readerror 0.01 --seqerror 0.001 --n_train 30000 --maxploidy 4.5 --minploidy 1.5 --normalcontamination --tumourheterogeneity --fast --tumourstatestable $tumourstatestable  --maxnormalcontamination 0.6 --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $file --outdir $outdir > $outdir/$samplename.out & 
end

#
