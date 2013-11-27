function output2file(ploidyno, seg, segall, params, options)

outfile = options.outfile;
n_lev = length(options.lambda_1_range);

disp(['Writing CNAs to file: ' outfile ]);
if ploidyno == 1
	fid = fopen(outfile, 'wt');
	fprintf(fid, 'Chromosome\tStartPosition\tEndPosition\tCopyNumber\tLOH\tRank\tLoglik\tnProbes\tNormalFraction\tTumourState\tPloidyNo\tMajorCopyNumber\tMinorCopyNumber\n');
else
	fid = fopen(outfile, 'at');
end

for lev = 1 : n_lev
	nseg = length(seg{lev});
	for ni = 1 : nseg
		chrNo = seg{lev}{ni}.chromosome;
		startPos = seg{lev}{ni}.startPos;
		endPos = seg{lev}{ni}.endPos;
		cn = seg{lev}{ni}.cn;
		loh = seg{lev}{ni}.loh;
		u = seg{lev}{ni}.u;
		loglik = seg{lev}{ni}.loglik;
		nprobes = seg{lev}{ni}.nprobes;
		ts = seg{lev}{ni}.ts; 
		majorcn = seg{lev}{ni}.majorcn; 
		minorcn = seg{lev}{ni}.minorcn; 
		fprintf(fid, '%d\t%d\t%d\t%d\t%d\t%d\t%15.0f\t%d\t%1.1f\t%2.0f\t%2.0f\t%2.0f\t%2.0f\n', chrNo, startPos, endPos, cn, loh, lev, loglik, nprobes, u, ts, ploidyno, majorcn, minorcn);
	end
end
fclose(fid);

disp(['Writing all events file:' options.outfile_all]);
if ploidyno == 1
	fid = fopen(options.outfile_all, 'wt');
	fprintf(fid, 'Chromosome\tStartPosition\tEndPosition\tCopyNumber\tLOH\tRank\tLoglik\tnProbes\tNormalFraction\tTumourState\tPloidyNo\tMajorCopyNumber\tMinorCopyNumber\n');
else
	fid = fopen(options.outfile_all, 'at');
end

nseg = length(segall);
for ni = 1 : nseg
	chrNo = segall{ni}.chromosome;
	startPos = segall{ni}.startPos;
	endPos = segall{ni}.endPos;
	cn = segall{ni}.cn;
	loh = segall{ni}.loh;
	u = segall{ni}.u;
	loglik = segall{ni}.loglik;
	nprobes = segall{ni}.nprobes;
	ts = segall{ni}.ts; 
	majorcn = segall{ni}.majorcn; 
	minorcn = segall{ni}.minorcn; 
	fprintf(fid, '%d\t%d\t%d\t%d\t%d\t%d\t%15.0f\t%d\t%1.1f\t%2.0f\t%2.0f\t%2.0f\t%2.0f\n', chrNo, startPos, endPos, cn, loh, lev, loglik, nprobes, u, ts, ploidyno, majorcn, minorcn);

end
fclose(fid);

