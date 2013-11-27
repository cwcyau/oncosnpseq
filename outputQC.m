function outputQC(params, options)

outfile_qc = options.outfile_qc;
nploidyvals = length(params.ploidyvals);

% write QC file
fid = fopen(outfile_qc, 'wt');

fprintf(fid, 'Haploid Read Depth\tNormal Fraction\tAverage Copy Number\tLog-likelihood\tRead Variance\tPloidyNo\n');
for i = 1 : nploidyvals
	fprintf(fid, '%d\t%1.1f\t%1.1f\t%f\t%f\t%2.0f\n', params.ploidyvals(i), params.ploidynormal(i), params.ploidycn(i), params.ploidycost(i), params.lambda_s, i);
end

fclose(fid);

