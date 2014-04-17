#!/bin/csh

set MCRDIR=/data/cyau/MCR/r2013b/v82/

set hgtable=../config/hgTables_b37.txt
set tumourstatestable=../config/tumourStatesCLL.txt

set seqtype=illumina
set indir=/data/cyau/oncosnpseq/data/sim/cll/
set outdir=/data/cyau/oncosnpseq/output/sim/

set readerror=0.01
set seqerror=0.001
set ntrain=30000
set maxploidy=2.5
set minploidy=1.5

foreach file(`find $indir -type f -name "$1*.oncoseq.gz"`)
	echo $file
	set samplename=`basename $file`
	echo $samplename
	rm $outdir/$samplename.out
	../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:30]" --chr_range "[1:18]" --lambda1 "[ 1000 500 250 100 50 30 10 ]" --maxploidy $maxploidy --minploidy $minploidy --tumourheterogeneity --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $file --outdir $outdir > $outdir/$samplename.out
end

