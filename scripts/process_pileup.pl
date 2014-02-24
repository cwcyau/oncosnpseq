use strict;

use Getopt::Long;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::File ;

my $infile = "";
my $outfile = "";
my $minReads = 0;
my $minaltread = 0;
my $snpfile = "";

my $result = GetOptions ( 	"infile=s" => \$infile,    # string
                       		"outfile=s"   => \$outfile,      # string
                       		"snpfile=s"   => \$snpfile,      # string
							"minreads=i"  => \$minReads, # numeric
							"minaltread=i" => \$minaltread);  # numeric		

$infile =~ m/^([a-zA-Z0-9\._\-\/]+)$/ or die "Bad data in infile argument";	
$outfile =~ m/^([a-zA-Z0-9\._\-\/]+)$/ or die "Bad data in outfile argument";	
$snpfile =~ m/^([a-zA-Z0-9\._\-\/]+)$/ or die "Bad data in snpfile argument";	
$minReads =~ m/^([0-9]+)$/ or die "Bad data in minreads argument";	
$minaltread =~ m/^([0-9]+)$/ or die "Bad data in minaltread argument";	

print "Reading file: $snpfile\n";
my %snpTable = ();
open(SNPFILE, $snpfile) or die "Cannot read $snpfile\n";

	<SNPFILE>; 
	
	while ( my $line = <SNPFILE> ) {
	
		my @linedat = split(/\t/, $line);
		
		my $chr = $linedat[1];
		my $chrStart = $linedat[2]; 
		my $chrEnd = $linedat[3];
		my $rsId = $linedat[4];

		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/XY/25/;		
		$chr =~ s/MT/26/;					
				
		my $id = "$chr:$chrEnd";
		$snpTable{$id} = 1;
	
	}

close(SNPFILE);

open(OUTFILE, ">", $outfile) or die "Cannot write to $outfile\n";

	print "Reading file: $infile\n";
	my $INFILEHANDLE;
	if ( $infile =~ /.bz2/ ) {
		print "Found bzip2 file ... decompressing ...\n";
		$INFILEHANDLE = new IO::Uncompress::Bunzip2 $infile or die "Cannot read $infile\n";
	} else {
		if ( $infile =~ /.gz/ ) {
		print "Found gzip file ... decompressing ...\n";
			$INFILEHANDLE = new IO::Uncompress::Gunzip $infile or die "Cannot read $infile\n";
		} 
		else {
			$INFILEHANDLE = new IO::File "< $infile" or die "Cannot read $infile\n";	
		}
	}

	print OUTFILE "Chr,Position,ReadRef,AltA,AltC,AltG,AltT\n";

	<$INFILEHANDLE>;

	my $linecount = 0;
	my $linesused = 0;
	my $lastread = 0;

	while ( my $line = <$INFILEHANDLE> ) {

		$linecount++;

		if ( $linecount % 100000 == 1 ) {
			print "Lines processed: $linecount. Lines used: $linesused\n";
		}

		chomp($line);

		$line =~ s/\r|\n//;

		my @linedat = split(/\t/, $line);

		my $chr = $linedat[0];			
		my $pos = $linedat[1];
		my $refAllele = $linedat[2];
		my $readBases = $linedat[3];
		
		my $id = "$chr:$pos";
		
		if ( !defined($snpTable{$id}) ) {
			next;
		}
		
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/XY/25/;		
		$chr =~ s/MT/26/;					
		
		if ( ( $readBases > $minReads ) & ( $chr < 25 ) ) {

			$_ = $linedat[4];
			my $refcount1 = tr/\./\./;
			my $refcount2 = tr/\,/\,/;
			my $refcount = $refcount1 + $refcount2;

			$_ = $linedat[4];
			my $Acount1 = tr/A/A/;
			my $Acount2 = tr/a/a/;		
			my $Acount = $Acount1 + $Acount2;

			$_ = $linedat[4];
			my $Ccount1 = tr/C/C/;
			my $Ccount2 = tr/c/c/;		
			my $Ccount = $Ccount1 + $Ccount2;

			$_ = $linedat[4];
			my $Gcount1 = tr/G/G/;
			my $Gcount2 = tr/g/g/;		
			my $Gcount = $Gcount1 + $Gcount2;

			$_ = $linedat[4];
			my $Tcount1 = tr/T/T/;
			my $Tcount2 = tr/t/t/;		
			my $Tcount = $Tcount1 + $Tcount2;	
			
			my $minaltread = int($readBases+0.5);
			if ( $minaltread < 1 ) {
				$minaltread = 1;
			}
					
			print OUTFILE "$chr,$pos,$refcount,$Acount,$Ccount,$Gcount,$Tcount\n";
			$linesused++;

		} 
		


	}

close(OUTFILE);

print "Lines processed: $linecount. Lines used: $linesused\n";



