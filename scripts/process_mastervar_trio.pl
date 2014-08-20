use strict;
use Getopt::Long;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error) ;
use IO::File ;

my $infile = "";
my $outfile = "";
my $result = GetOptions ( 	"infile=s" => \$infile,    # string
                       		"outfile=s"   => \$outfile );
$infile =~ m/^([a-zA-Z0-9\._\/\-]+)$/ or die "Bad data in infile argument";	
$outfile =~ m/^([a-zA-Z0-9\._\/\-]+)$/ or die "Bad data in outfile argument";	

my($chrInd, $startInd, $endInd, $allele1SeqInd, $allele2SeqInd, $allele1VarQualInd, $allele2VarQualInd, $allele1ReadCountT1Ind, $allele2ReadCountT1Ind, $allele1ReadCountT2Ind, $allele2ReadCountT2Ind, $allele1ReadCountInd, $allele2ReadCountInd );

my($tumour1count,$tumour1total,$tumour2count,$tumour2total,$normal1count,$normal2total);

my $INFILEHANDLE;
if ( $infile =~ /.bz2/ ) {
	$INFILEHANDLE = new IO::Uncompress::Bunzip2 $infile;
} else {
	$INFILEHANDLE = new IO::File "< $infile";	
}

open(OUTFILE_T1, ">", "$outfile.T1") or die "Cannot write to: $outfile.T1\n";
open(OUTFILE_T2, ">", "$outfile.T2") or die "Cannot write to: $outfile.T2\n";

	print OUTFILE_T1 "Chr\tPosition\tAllele1VarQual\tAllele2VarQual\tTumourAllele1Count\tTumourTotalCount\tNormalAllele1Count\tNormalTotalCount\n";
	print OUTFILE_T2 "Chr\tPosition\tAllele1VarQual\tAllele2VarQual\tTumourAllele1Count\tTumourTotalCount\tNormalAllele1Count\tNormalTotalCount\n";

	my $lineCount = 0;
	my $snpCount = 0;
	while ( my $line = <$INFILEHANDLE> ) {

		$lineCount++;			
		if ( $lineCount%100000 == 1 ) {
			print "Lines processed: $lineCount. SNPs found: $snpCount\n";
		}

		chomp($line);
		$line =~ s/\r|\n//;

		if ( $line =~ /#/ ) {
			next;
		}
		if ( $line =~ />/ ) {
			my @linedat = split(/\t/, $line);
			my $indexCount = 0;
			foreach my $cell ( @linedat ) {
				if ( $cell =~ /chromosome/ ) {
					$chrInd = $indexCount;
				}
				if ( $cell =~ /begin/ ) {
					$startInd = $indexCount;
				}
				if ( $cell =~ /end/ ) {
					$endInd = $indexCount;	
				}
				if ( $cell =~ /allele1Seq/ ) {
					$allele1SeqInd = $indexCount;
				}
				if ( $cell =~ /allele2Seq/ ) {
					$allele2SeqInd = $indexCount;
				}						
				if ( ( $cell =~ /allele1VarFilter/ ) | ( $cell =~ /allele1VarQuality/ ) ) {
					$allele1VarQualInd = $indexCount;
				}
				if ( ( $cell =~ /allele2VarFilter/ ) | ( $cell =~ /allele2VarQuality/ ) ) {
					$allele2VarQualInd = $indexCount;
				}				
				if ( $cell =~ /^referenceAlleleReadCount-T1$/ ) {
					$allele1ReadCountT1Ind = $indexCount;	
				}
				if ( $cell =~ /^totalReadCount-T1$/ ) {
					$allele2ReadCountT1Ind = $indexCount;
				}				
				if ( $cell =~ /^referenceAlleleReadCount-T2$/ ) {
					$allele1ReadCountT2Ind = $indexCount;	
				}
				if ( $cell =~ /^totalReadCount-T2$/ ) {
					$allele2ReadCountT2Ind = $indexCount;
				}				
				if ( $cell =~ /^referenceAlleleReadCount$/ ) {
					$allele1ReadCountInd = $indexCount;
				}
				if ( $cell =~ /^totalReadCount$/ ) {
					$allele2ReadCountInd = $indexCount;
				}
												
				$indexCount++;
			}	
			
			next;
			
		}
		if ( $line !~ /dbsnp/ ) {
			next;
		}

		my @linedat = split(/\t/, $line);
	
		my $chr = $linedat[$chrInd];
		my $startPos = $linedat[$startInd];
		my $endPos = $linedat[$endInd];
		my $allele1Seq = $linedat[$allele1SeqInd];
		my $allele2Seq = $linedat[$allele2SeqInd];		
		my $allele1VarQual = $linedat[$allele1VarQualInd];
		my $allele2VarQual = $linedat[$allele2VarQualInd];	
		my $allele1ReadCountT1 = $linedat[$allele1ReadCountT1Ind];
		my $allele2ReadCountT1 = $linedat[$allele2ReadCountT1Ind];	
		my $allele1ReadCountT2 = $linedat[$allele1ReadCountT2Ind];
		my $allele2ReadCountT2 = $linedat[$allele2ReadCountT2Ind];					
		my $allele1ReadCount = $linedat[$allele1ReadCountInd];				
		my $allele2ReadCount = $linedat[$allele2ReadCountInd];

		my $tumour1count = $allele1ReadCountT1;
		my $tumour1total = $allele2ReadCountT1;
		
		my $tumour2count = $allele1ReadCountT2;
		my $tumour2total = $allele2ReadCountT2;
						
		my $normalCount = $allele1ReadCount;
		my $normalTotal = $allele2ReadCount;	
				
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/M/25/;

		my $allele1Qual = 0;
		my $allele2Qual = 0;

		if ( $allele1VarQual =~ /VQHIGH/ ) {
			$allele1Qual = 1;
		}
		if ( $allele1VarQual =~ /PASS/ ) {
			$allele1Qual = 1;
		}
		if ( $allele1VarQual =~ /VQLOW/ ) {
			$allele1Qual = 0;
		}
		if ( $allele2VarQual =~ /VQHIGH/ ) {
			$allele2Qual = 1;
		}
		if ( $allele2VarQual =~ /PASS/ ) {
			$allele2Qual = 1;
		}
		if ( $allele2VarQual =~ /VQLOW/ ) {
			$allele2Qual = 0;
		}
								
	
		print OUTFILE_T1 "$chr\t$startPos\t$allele1Qual\t$allele2Qual\t$tumour1count\t$tumour1total\t$normalCount\t$normalTotal\n";

		print OUTFILE_T2 "$chr\t$startPos\t$allele1Qual\t$allele2Qual\t$tumour2count\t$tumour2total\t$normalCount\t$normalTotal\n";
	
		$snpCount++;

	}

close(OUTFILE_T1);
close(OUTFILE_T2);



