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

my($chrInd, $startInd, $endInd, $refSeqInd, $allele1SeqInd, $allele2SeqInd, $allele1VarQualInd, $allele2VarQualInd, $allele1ReadCountN1Ind, $allele2ReadCountN1Ind, $allele1ReadCountInd, $allele2ReadCountInd, $somaticCategoryInd, $somaticRankInd );

my($tumour1count,$tumour2total,$normal1count,$normal2total);

my $INFILEHANDLE;
if ( $infile =~ /.bz2/ ) {
	$INFILEHANDLE = new IO::Uncompress::Bunzip2 $infile;
} else {
	$INFILEHANDLE = new IO::File "< $infile";	
}

my %subTable = ();
foreach my $refSeq ( "A", "C", "G", "T") {
	foreach my $alleleSeq ( "A", "C", "G", "T") {
		$subTable{$refSeq}{$alleleSeq} = 0;
	}
}

open(OUTFILE, ">", $outfile) or die "Cannot write to: $outfile\n";

	print OUTFILE "Chr\tPosition\tSomaticRank\tRefAllele\tAltAllele\tTumourAltCount\tTumourRefCount\tNormalAltCount\tNormalRefCount\n";

	my $lineCount = 0;
	my $snpCount = 0;
	while ( my $line = <$INFILEHANDLE> ) {

		$lineCount++;			
		if ( $lineCount%100000 == 1 ) {
			print "Lines processed: $lineCount. SNPs found: $snpCount\n";
			
			foreach my $refSeq ( "A", "C", "G", "T") {
				foreach my $alleleSeq ( "A", "C", "G", "T") {
					print "($refSeq,$alleleSeq): $subTable{$refSeq}{$alleleSeq}\n";
				}
			}			
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
			
				print "$indexCount: $cell\n";
			
				if ( $cell =~ /chromosome/ ) {
					$chrInd = $indexCount;
				}
				if ( $cell =~ /begin/ ) {
					$startInd = $indexCount;
				}
				if ( $cell =~ /end/ ) {
					$endInd = $indexCount;	
				}
				if ( $cell =~ /^reference$/ ) {
					$refSeqInd = $indexCount;
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
				if ( $cell =~ /^somaticCategory$/ ) {
					$somaticCategoryInd = $indexCount;
				}
				if ( $cell =~ /^somaticRank$/ ) {
					$somaticRankInd = $indexCount;
				}						
												
				$indexCount++;
			}	
			
			next;
			
		}
		
		if ( $line =~ /dbsnp/ ) {
			next;
		}
		if ( $line =~ /sub/ ) {
			next;
		}
		if ( $line =~ /ins/ ) {
			next;
		}
		if ( $line =~ /del/ ) {
			next;
		}
		
		my @linedat = split(/\t/, $line);
	
		my $chr = $linedat[$chrInd];
		my $startPos = $linedat[$startInd];
		my $endPos = $linedat[$endInd];
		my $refSeq = $linedat[$refSeqInd];
		my $allele1Seq = $linedat[$allele1SeqInd];
		my $allele2Seq = $linedat[$allele2SeqInd];		
		my $allele1VarQual = $linedat[$allele1VarQualInd];
		my $allele2VarQual = $linedat[$allele2VarQualInd];	
		my $allele1ReadCountN1 = $linedat[$allele1ReadCountN1Ind];
		my $allele2ReadCountN1 = $linedat[$allele2ReadCountN1Ind];					
		my $allele1ReadCount = $linedat[$allele1ReadCountInd];				
		my $allele2ReadCount = $linedat[$allele2ReadCountInd];
		my $somaticCategory = $linedat[$somaticCategoryInd];
		my $somaticRank = $linedat[$somaticRankInd];
						
		if ( $somaticCategory ne "snp" ) {
			next;
		}
		if ( $allele1VarQual !~ /PASS/ ) {
			next;
		};
		if ( $allele2VarQual !~ /PASS/ ) {
			next;
		}		
		
		my $tumourAltcount = 0;
		my $tumourReftotal = 0;
		my $normalAltcount = 0;
		my $normalReftotal = 0;
		
		if ( $refSeq eq $allele1Seq  ) {
			$tumourAltcount = $allele2ReadCount;
			$tumourReftotal = $allele1ReadCount;					
			$normalAltcount = $allele2ReadCountN1;
			$normalReftotal = $allele1ReadCountN1;					
		} else {
			$tumourAltcount = $allele1ReadCount;
			$tumourReftotal = $allele2ReadCount;					
			$normalAltcount = $allele1ReadCountN1;
			$normalReftotal = $allele2ReadCountN1;					
		}
				
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/M/25/;

		my $alleleAltSeq;

		if ( $refSeq eq $allele1Seq ) {
			$subTable{$refSeq}{$allele2Seq}++;
			$alleleAltSeq = $allele2Seq;
		}
		if ( $refSeq eq $allele2Seq ) {
			$subTable{$refSeq}{$allele1Seq}++;
			$alleleAltSeq = $allele1Seq;			
		}
		

		#print "$line\n";
		print OUTFILE "$chr\t$startPos\t$somaticRank\t$refSeq\t$alleleAltSeq\t$tumourAltcount\t$tumourReftotal\t$normalAltcount\t$normalReftotal\n";
	
		$snpCount++;

	}

close(OUTFILE);


