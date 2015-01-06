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

my($chrInd, $startInd, $endInd, $allele1SeqInd, $allele2SeqInd, $allele1VarQualInd, $allele2VarQualInd, $refAlleleReadCountInd, $totalReadCountInd );

my($refAlleleReadCount,$totalReadCount,$normal1count,$normal2total);

print "Reading file: $infile\n";
my $INFILEHANDLE;
if ( $infile =~ /.bz2/ ) {
	print "Found bzip2 file ... decompressing ...\n";
	$INFILEHANDLE = new IO::Uncompress::Bunzip2 $infile or die "Cannot read $infile\n";
} else {
	$INFILEHANDLE = new IO::File "< $infile" or die "Cannot read $infile\n";	
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
				if ( $cell =~ /^referenceAlleleReadCount$/ ) {
					$refAlleleReadCountInd = $indexCount;
				}
				if ( $cell =~ /^totalReadCount$/ ) {
					$totalReadCountInd = $indexCount;
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
		my $refAlleleReadCount = $linedat[$refAlleleReadCountInd];				
		my $totalReadCount = $linedat[$totalReadCountInd];
				
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/M/25/;

		$allele1VarQual =~ s/PASS/1/;
		$allele1VarQual =~ s/VQLOW/0/;
		$allele1VarQual =~ s/AMBIGUOUS/0/;
							
		$allele2VarQual =~ s/PASS/1/;
		$allele2VarQual =~ s/VQLOW/0/;
		$allele2VarQual =~ s/AMBIGUOUS/0/;			
	
		if ( $allele1VarQual != 1 ) {
			next;
		}

		if ( $allele2VarQual != 1 ) {
			next;
		}
	
		print OUTFILE "$chr\t$startPos\t1\t1\t$refAlleleReadCount\t$totalReadCount\t0\t0\n";
	
		$snpCount++;

	}

close(OUTFILE);


