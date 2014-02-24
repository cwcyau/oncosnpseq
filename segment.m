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

[arrayind, log_nu, log_transMat] = gentransmat(-lambda_1_range(1), -lambda_1_range(1), options.tumourState, params.u_range);

for chrNo = chrRange
	for armNo = 1 : 2
		chrloc = find( chr == chrNo & arm == armNo );
		n_chr = length(chrloc);
		if n_chr > 0
			vpath = viterbimex(log_nu, loglik(:, chrloc), log_transMat); 
			x{1}(chrloc) = arrayind(vpath, 1);
			u{1}(chrloc) = params.u0 + (1-params.u0)*params.u_range(arrayind(vpath, 2));
		end
	end
end

for lev = 2 : n_lev

	lambda_1 = lambda_1_range(lev);

	[arrayind, log_nu, log_transMat] = gentransmat(-lambda_1, -lambda_1, options.tumourState, params.u_range);

	for chrNo = chrRange
		for armNo = 1 : 2
			chrloc = find( chr == chrNo & arm == armNo );
			n_chr = length(chrloc);
			if n_chr > 0
				vpath = viterbimex(log_nu, loglik(:, chrloc), log_transMat); 
				x{lev}(chrloc) = arrayind(vpath, 1);				
				u{lev}(chrloc) = params.u0 + (1-params.u0)*params.u_range(arrayind(vpath, 2));
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

