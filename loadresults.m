clear all;
close all;
clc;

options.infile = '/data/cyau/ahmed/oncoseq/prechemo-all.oncoseq';
options.chrRange = 1:22;
options.hgtables = 'config/hgTables_b37.txt';
options.gcdir = '/data/cyau/software/gc/b37/';

cnvfile = '/data/cyau/ahmed/oncoseq/output/prechemo-all.cnvs';
outfile = '/data/cyau/ahmed/oncoseq/output/prechemo-all.summary.txt';

R = 5;
ploidyno = 1;

% load data
[chr, arm, pos, k, d, dd] = loaddataCG(options);

chrloc = [];
pos_chr = [];
for chrNo = 1 : 22
	chrloc{chrNo} = find( chr == chrNo );
	pos_chr{chrNo} = pos(chrloc{chrNo});
end

% load data files

N = length(chr);

cn_all = zeros(N, R, 'int8');
loh_all = zeros(N, R, 'int8');
majorcn_all = zeros(N, R, 'int8');
minorcn_all = zeros(N, R, 'int8');

		
[ cnachr, startPos, endPos, copyNumber, loh, rank, loglik_cna, nProbes, normalFraction, tumourState, ploidyNo, majorCN, minorCN ] = ...
	textread( cnvfile, '%n %n %n %n %n %n  %n  %n  %n  %n  %n  %n  %n', 'headerlines', 1);
	
nCNA = length(cnachr);

for i = 1 : nCNA

	if ploidyNo(i) == ploidyno

		chrNo = cnachr(i);

		loc = find( pos_chr{chrNo} >= startPos(i) & pos_chr{chrNo} <= endPos(i) );
	
		cn_all( chrloc{chrNo}(loc), rank(i):R ) = copyNumber(i);
		loh_all( chrloc{chrNo}(loc), rank(i):R ) = loh(i);
		majorcn_all( chrloc{chrNo}(loc), rank(i):R ) = majorCN(i);
		minorcn_all( chrloc{chrNo}(loc), rank(i):R ) = minorCN(i);
	
	end

end		


fid = fopen(outfile, 'wt');
fprintf(fid, '# Ploidy Number: %d\n', ploidyno); 
fprintf(fid, 'Chromosome\tPosition\tTotalCopyNumber\tMajorCopyNumber\tMinorCopyNumber\tLOH\tRank\n');
for i = 1 : N
	for rankno = 1 : R
		fprintf(fid, '%2.0f\t%15.0f\t%d\t%d\t%d\t%d\t%d\n', chr(i), pos(i), cn_all(i, rankno), majorcn_all(i, rankno), minorcn_all(i, rankno), loh_all(i, rankno), rankno);
	end
end
fclose(fid);



