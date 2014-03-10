#!/bin/csh

set MCRDIR=/data/cyau/MCR/r2013b/v82/

set hgtable=../config/hgTables_b37.txt

set seqtype=illumina
set indir=/data/cyau/wgs500/pileup/snp/
set outdir=/data/cyau/wgs500/output/v2.0/

set readerror=0.05
set seqerror=0.05
set ntrain=50000

mkdir -p $outdir

set tumourstatestable=../config/tumourStates.txt


foreach file(`find $indir -type f -name "FH40*T.txt*"`)
	echo $file
	set FILE1=`basename $file`
	set SAMPLE=`echo $FILE1 | sed -E  's/(.*)T.txt(.*)/\1/'`
	echo $SAMPLE
	set FILE2="$SAMPLE"N.txt.gz
	echo $FILE1
	echo $FILE2
	../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:22]" --readerror $readerror --seqerror $seqerror --n_train $ntrain --maxploidy 4.5 --minploidy 1.5 --normalcontamination --tumourheterogeneity --tumourstatestable $tumourstatestable  --maxnormalcontamination 0.6 --hgtable $hgtable --seqtype $seqtype --samplename $SAMPLE --infile $indir/$FILE1 --normalfile $indir/$FILE2  --outdir $outdir > $outdir/$SAMPLE.out
end

