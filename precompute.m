% precompute

dd_unique = 0:dd_max;
n_dd_unique = length(dd_unique);

ut_range = 0:0.1:1.0;
n_ut = length(ut_range);


loglik_r_mat = zeros(S, n_ut, n_dd_unique);
for si = 1 : S


	% normal/tumour copy number
	cn_t = tumourState(si, 4);
	cn_n = normalState(si, 4);

	for ut_i = 1 : n_ut
	
		ut = ut_range(ut_i);
		un = 1 - ut;
		
		m = un*read_depth*cn_n + ut*cn_t*read_depth;
		tmp(1, :) = log(1-seq_error) + logtpdf(dd_unique, m, lambda_s^2, nu);
		tmp(2, :) = log(seq_error) - log(dd_max);
	
		loglik_r_mat(si, ut_i, :) = logsumexp(tmp, 1);
	
	end

end



