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

my($chrInd, $startInd, $endInd, $allele1SeqInd, $allele2SeqInd, $allele1VarQualInd, $allele2VarQualInd, $allele1ReadCountN1Ind, $allele2ReadCountN1Ind, $allele1ReadCountInd, $allele2ReadCountInd );

my($tumour1count,$tumour2total,$normal1count,$normal2total);

my $INFILEHANDLE;
if ( $infile =~ /.bz2/ ) {
	$INFILEHANDLE = new IO::Uncompress::Bunzip2 $infile;
} else {
	$INFILEHANDLE = new IO::File "< $infile";	
}

open(OUTFILE, ">", $outfile) or die "Cannot write to: $outfile\n";

	print OUTFILE "Chr\tPosition\tAllele1VarQual\tAllele2VarQual\tTumourAllele1Count\tTumourTotalCount\tNormalAllele1Count\tNormalTotalCount\n";

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
				if ( $cell =~ /allele1VarFilter/ ) {
					$allele1VarQualInd = $indexCount;
				}
				if ( $cell =~ /allele2VarFilter/ ) {
					$allele2VarQualInd = $indexCount;
				}				
				if ( $cell =~ /^allele1ReadCount-N1$/ ) {
					$allele1ReadCountN1Ind = $indexCount;	
				}
				if ( $cell =~ /^allele2ReadCount-N1$/ ) {
					$allele2ReadCountN1Ind = $indexCount;
				}				
				if ( $cell =~ /^allele1ReadCount$/ ) {
					$allele1ReadCountInd = $indexCount;
				}
				if ( $cell =~ /^allele2ReadCount$/ ) {
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
		my $allele1ReadCountN1 = $linedat[$allele1ReadCountN1Ind];
		my $allele2ReadCountN1 = $linedat[$allele2ReadCountN1Ind];					
		my $allele1ReadCount = $linedat[$allele1ReadCountInd];				
		my $allele2ReadCount = $linedat[$allele2ReadCountInd];

		if ( $allele1Seq ne $allele2Seq ) {
			$tumour1count = $allele1ReadCount;
			$tumour2total = $allele1ReadCount + $allele2ReadCount;					
			$normal1count = $allele1ReadCountN1;
			$normal2total = $allele1ReadCountN1 + $allele2ReadCountN1;					
		} else {
			$tumour1count = $allele1ReadCount;
			$tumour2total = $allele1ReadCount + $allele2ReadCount;					
			$normal1count = $allele1ReadCountN1;
			$normal2total = $allele1ReadCountN1 + $allele2ReadCountN1;					
		}
				
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
								
	
		print OUTFILE "$chr\t$startPos\t$allele1Qual\t$allele2Qual\t$tumour1count\t$tumour2total\t$normal1count\t$normal2total\n";
	
		$snpCount++;

	}

close(OUTFILE);


