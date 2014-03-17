function outputQC(params, options)

disp(['Writing QC metrics to: ' options.outfile_qc]);

outfile_qc = options.outfile_qc;
nploidyvals = length(params.ploidyvals);

% write QC file
fid = fopen(outfile_qc, 'wt');

fprintf(fid, 'HaploidCoverage\tNormalFraction\tAverageCopyNumber\tLoglikelihood\tReadVariance\tPloidyNo\n');
for i = 1 : nploidyvals
	fprintf(fid, '%d\t%1.1f\t%1.1f\t%f\t%f\t%2.0f\n', params.ploidyvals(i), params.ploidynormal(i), params.ploidycn(i), params.ploidycost(i), params.lambda_s, i);
end

fclose(fid);

