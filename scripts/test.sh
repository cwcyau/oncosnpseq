infile=/data/cyau/wgs500/pileup/AS_CLL_003D.bam.pileup.gz
outfile=/data/cyau/temp/test.out

snpfile=illumina_hg19.txt

perl -T process_pileup.pl --infile=$infile --outfile=$outfile --snpfile=$snpfile --minreads=0 --minaltread=0

