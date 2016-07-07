#!/bin/csh

set MCRDIR=/data/cyau/MCR/r2013b/v82/

set hgtable=../config/hgTables_b37.txt
set tumourstatestable=../config/tumourStatesSPS.txt

set seqtype=illumina
set indir=/data/cyau/downloads/allele_counts/
set outdir=/data/cyau/downloads/allele_counts/out/

set readerror=0.01
set seqerror=0.1
set ntrain=30000
set maxploidy=3.5
set minploidy=1.5


#Tumour	Matched Normal
#SPS_0001	SPS_0002

#Tumour	Matched Normal
#SPS_0004	SPS_0005

#Tumour	Matched Normal
#SPS_0006	SPS_0007

#Tumour	Matched Normal
#SPS_0008	SPS_0009

#Unmatched Tumour 
#SPS_0003

set samplename=AW_SPS_0001
set tumourfile=$indir/AW_SPS_0001.1M_markers_AlleleCounts.txt
set normalfile=$indir/AW_SPS_0002.1M_markers_AlleleCounts.txt

../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:24]" --lambda1 "[ 1000 500 250 100 50 30 10 ]" --maxploidy $maxploidy --minploidy $minploidy --normalcontamination --maxnormalcontamination 0.4 --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $tumourfile  --outdir $outdir > $outdir/$samplename.out &

set samplename=AW_SPS_0004
set tumourfile=$indir/AW_SPS_0004.1M_markers_AlleleCounts.txt
set normalfile=$indir/AW_SPS_0005.1M_markers_AlleleCounts.txt

../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:24]" --lambda1 "[ 1000 500 250 100 50 30 10 ]" --maxploidy $maxploidy --minploidy $minploidy --normalcontamination --maxnormalcontamination 0.4 --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $tumourfile  --outdir $outdir > $outdir/$samplename.out &

set samplename=AW_SPS_0006
set tumourfile=$indir/AW_SPS_0006.1M_markers_AlleleCounts.txt
set normalfile=$indir/AW_SPS_0007.1M_markers_AlleleCounts.txt

../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:24]" --lambda1 "[ 1000 500 250 100 50 30 10 ]" --maxploidy $maxploidy --minploidy $minploidy --normalcontamination --maxnormalcontamination 0.4 --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $tumourfile --outdir $outdir > $outdir/$samplename.out &

set samplename=AW_SPS_0008
set tumourfile=$indir/AW_SPS_0008.1M_markers_AlleleCounts.txt
set normalfile=$indir/AW_SPS_0009.1M_markers_AlleleCounts.txt

../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:24]" --lambda1 "[ 1000 500 250 100 50 30 10 ]" --maxploidy $maxploidy --minploidy $minploidy --normalcontamination --maxnormalcontamination 0.4 --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $tumourfile  --outdir $outdir > $outdir/$samplename.out &

set samplename=AW_SPS_0003
set tumourfile=$indir/AW_SPS_0003.1M_markers_AlleleCounts.txt

../executables/run_oncoseq.sh $MCRDIR --read_depth_range "[10:40]" --chr_range "[1:24]" --lambda1 "[ 1000 500 250 100 50 30 10 ]" --maxploidy $maxploidy --minploidy $minploidy --normalcontamination --maxnormalcontamination 0.4 --tumourstatestable $tumourstatestable --hgtable $hgtable --seqtype $seqtype --samplename $samplename --infile $tumourfile  --outdir $outdir > $outdir/$samplename.out &





