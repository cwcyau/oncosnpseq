#!/bin/csh

set MCRDIR=/data/cyau/MCR/r2013b/v82/

set hgtable=../config/hgTables_b37.txt

set seqtype=illumina
set indir=/data/cyau/wgs500/pileup/snp/
set outdir=/data/cyau/wgs500/output/v2.0/

set readerror=0.01
set seqerror=0.001
set lambda3=0.1
set ntrain=30000

mkdir -p $outdir

set tumourstatestable=../config/tumourStatesCLL.txt

set normalfile=$indir/AS_CLL_003GL.txt.gz
foreach file(`find $indir -type f -name "AS_CLL_003*"`)
	echo $file
	set samplename=`basename $file`
	echo $samplename
	rm $outdir/$samplename.out
	../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[20:40]" --chr_range "[1:23]" --readerror $readerror --seqerror $seqerror --lambda3 $lambda3 --n_train $ntrain --maxploidy 3.5 --minploidy 1.5 --tumourheterogeneity  --fast --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $file --normalfile $normalfile  --outdir $outdir > $outdir/$samplename.out
end


set normalfile=$indir/AS_CLL_006GL.txt.gz
foreach file(`find $indir -type f -name "AS_CLL_006*"`)
	echo $file
	set samplename=`basename $file`
	echo $samplename
	rm $outdir/$samplename.out
	../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[20:40]" --chr_range "[1:23]" --readerror $readerror --seqerror $seqerror --lambda3 $lambda3 --n_train $ntrain --maxploidy 3.5 --minploidy 1.5 --tumourheterogeneity  --fast --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $file --normalfile $normalfile  --outdir $outdir > $outdir/$samplename.out
end


set normalfile=$indir/AS_CLL_077GL.txt.gz
foreach file(`find $indir -type f -name "AS_CLL_077*"`)
	echo $file
	set samplename=`basename $file`
	echo $samplename
	rm $outdir/$samplename.out
	../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[20:40]" --chr_range "[1:23]" --readerror $readerror --seqerror $seqerror --lambda3 $lambda3 --n_train $ntrain --maxploidy 3.5 --minploidy 1.5 --tumourheterogeneity  --fast --tumourstatestable $tumourstatestable  --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $file --normalfile $normalfile  --outdir $outdir > $outdir/$samplename.out
end
