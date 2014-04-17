#!/bin/csh

set MCRDIR=/data/cyau/MCR/r2013b/v82/

set hgtable=../config/hgTables_b37.txt

set seqtype=illumina
set indir=/data/cyau/wgs500/pileup/snp/
set outdir=/data/cyau/wgs500/output/v2.0/

set readerror=0.01
set seqerror=0.001
set lambda3=0.05
set ntrain=30000

mkdir -p $outdir

set tumourstatestable=../config/tumourStatesOncoSNPSEQ.txt


foreach file(`find $indir -type f -name "FH_BLC_*.txt*"`)
	echo $file
	set FILE1=`basename $file`
	set SAMPLE=`echo $FILE1 | sed -E  's/(.*)T.txt(.*)/\1/'`
	echo $SAMPLE
	set FILE2="$SAMPLE"N.txt.gz
	echo $FILE1
	echo $FILE2
	../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:50]" --chr_range "[1:23]" --readerror $readerror --seqerror $seqerror --lambda3 $lambda3 --n_train $ntrain --maxploidy 5.5 --minploidy 1.5 --normalcontamination --tumourheterogeneity --fast --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $SAMPLE --infile $indir/$FILE1 --normalfile $indir/$FILE2  --outdir $outdir > $outdir/$SAMPLE.out 
end

