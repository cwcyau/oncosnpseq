function [x, seg, u, seg_all] = segment(chr, arm, pos, k, d, dd, log_pr_gg, loglik, params, options)

lambda_1_range = options.lambda_1_range;
lambda_2 = options.lambda_2;
chrRange = options.chrRange;
tumourState = options.tumourState;

n_lev = length(lambda_1_range);
[S, N] = size(loglik);

for lev = 1 : n_lev
	x{lev} = zeros(1, N);
	u{lev} = zeros(1, N);
end

for chrNo = chrRange
	for armNo = 1 : 2
		chrloc = find( chr == chrNo & arm == armNo );
		n_chr = length(chrloc);
		if n_chr > 0
			x{1}(chrloc) = optimPath( loglik(:, chrloc), lambda_1_range(1) );
			u{1}(chrloc) = calclikelihoodU(k(chrloc), d(chrloc), dd(chrloc), log_pr_gg(:, chrloc), x{1}(chrloc), params, options);
		end
	end
end

for lev = 2 : n_lev

	lambda_1 = lambda_1_range(lev);

	for chrNo = chrRange
		for armNo = 1 : 2
			chrloc = find( chr == chrNo & arm == armNo );
			n_chr = length(chrloc);
			if n_chr > 0
				x{lev}(chrloc) = optimMultiPath( loglik(:, chrloc), lambda_1, x{lev-1}(chrloc)-1, lambda_2);				
				u{lev}(chrloc) = calclikelihoodU(k(chrloc), d(chrloc), dd(chrloc), log_pr_gg(:, chrloc), x{lev}(chrloc), params, options);
			end
		end
	end

end

seg{1} = findsegments(chr, arm, pos, x{1}, u{1}, loglik, [], options, params);
for lev = 2 : n_lev
	xd = x{lev} ~= x{lev-1};
	if sum(xd) == 0
		seg{lev} = [];
	else
		seg{lev} = findmultisegments(chr, arm, pos, xd, x{lev}, x{lev-1}, u{lev}, loglik, [], options, params);
	end
end

seg_all = findsegments(chr, arm, pos, x{n_lev}, u{n_lev}, loglik, [], options, params);

