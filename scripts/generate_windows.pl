use POSIX;

$snpfile = "illumina_hg19.txt";
$outfile = "windows.nochr.bed";
$hgtables = "../config/hgTables_b37.txt";

print "Reading chromosomal information ...\n";
%hginfo = ();
open(HGFILE, $hgtables) or die "Cannot read $hgtables\n";

	<HGFILE>;

	while ( $line = <HGFILE> ) {

		chomp($line);

		( $chrNo, $startPos, $endPos, $armNo ) = split(/\t/, $line);

		$hginfo{$chrNo}{$armNo}->{startPos} = $startPos;
		$hginfo{$chrNo}{$armNo}->{endPos} = $endPos;

		print "($chrNo, $armNo): $hginfo{$chrNo}{$armNo}->{startPos}, $hginfo{$chrNo}{$armNo}->{endPos}\n";
	}

close(HGFILE);

print "Reading SNPs ...\n";
%snptable = ();
open(SNPFILE, $snpfile) or die "Cannot open $SNPFILE\n";

	<SNPFILE>;

	while ( $line = <SNPFILE> ) {

		chomp($line);

		( $bin, $chrNo, $startPos, $endPos, $name, $tmp ) = split(/\t/, $line);

		$chrNo =~ s/chr//;
		$chrNo =~ s/X/23/;
		$chrNo =~ s/Y/24/;

		if ( ( $endPos >= $hginfo{$chrNo}{1}->{startPos} ) & ( $endPos <= $hginfo{$chrNo}{1}->{endPos} ) ) {
			$snptable{$chrNo}{1}{$endPos} = $name; 
		}
		if ( ( $endPos >= $hginfo{$chrNo}{2}->{startPos} ) & ( $endPos <= $hginfo{$chrNo}{2}->{endPos} ) ) {
			$snptable{$chrNo}{2}{$endPos} = $name; 
		}

	}

close(SNPFILE);

print "Generating windows ...\n";
open(OUTFILE, ">", $outfile) or die "Cannot write to $outfile\n";

	for ( $chrNo = 1; $chrNo < 25; $chrNo++ ) {
		
		for ( $armNo = 1; $armNo < 3; $armNo++ ) {

			@pos = keys %{ $snptable{$chrNo}{$armNo} };

			@sort_by_pos = sort { $pos[$a] <=> $pos[$b] } 0 .. $#pos;

			@names = values %{ $snptable{$chrNo}{$armNo} };

			@pos = @pos[@sort_by_pos];
			@names = @names[@sort_by_pos];

			$nSnps = @pos;

			print "$chrNo\t$armNo\t$nSnps\n";

			if ( $nSnps > 5 ) {

				$startChrPos = $hginfo{$chrNo}{$armNo}->{startPos};
				$endChrPos = $hginfo{$chrNo}{$armNo}->{endPos};

				$windowStart = $startChrPos;
				$windowEnd = ceil( 0.5*($pos[0]+$pos[1]) );
				$name = $names[0];

#				print "$chrNo\t$windowStart\t$windowEnd\t$name\n";
				print OUTFILE "$chrNo\t$windowStart\t$windowEnd\t$name\n";

				for ( $i = 1; $i < $nSnps-1; $i++ ) {
					$windowStart = ceil( 0.5*($pos[$i-1]+$pos[$i]) )+1;
					$windowEnd = ceil( 0.5*($pos[$i]+$pos[$i+1]) );
					$name = $names[$i];
#					print "chr$chrNo\t$windowStart\t$windowEnd\t$name\n";
					print OUTFILE "$chrNo\t$windowStart\t$windowEnd\t$name\n";
				}

				$windowStart = ceil( 0.5*($pos[$nSnps-2]+$pos[$nSnps-1]) )+1;
				$windowEnd = $endChrPos;
				$name = $names[$nSnps-1];

#				print "chr$chrNo\t$windowStart\t$windowEnd\t$name\n";
				print OUTFILE "$chrNo\t$windowStart\t$windowEnd\t$name\n";

			}

		}		

	}

close(OUTFILE);

print "All done.\n";
