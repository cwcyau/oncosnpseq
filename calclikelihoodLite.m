%
% OncoSEQ
%
% Written by: Christopher Yau
% 
% 19 March 2012
%
% Copyright (C) 2012 Imperial College London
%
function log_pr_s = calclikelihoodLite(k, d, dd, log_pr_gg, params, options)

N = length(d);
S = params.S;
K = params.K;
G = params.G;
U = params.U;

lambda_s = params.lambda_s;
nu = params.nu;
u = params.u_range;
p_u = params.p_u;
u0 = params.u0;
read_error = params.read_error;
read_depth = params.read_depth;
seq_error = params.seq_error;

tumourState = options.tumourState;
normalState = options.normalState;
chrRange = options.chrRange;
genotypeState = options.genotypeState;
pg = options.pg;

binomcoeff = gammaln(d+1) - gammaln(k+1) - gammaln(d-k+1);
dd_max = max(dd);

dloc = find( d > 0 );
binomcoeff_dloc = binomcoeff(dloc); 
k_dloc = k(dloc);
d_minus_k_dloc = d(dloc)-k(dloc);

loglik_r = -Inf*ones(2, N);
loglik_r(2, :) = log(seq_error) - log(dd_max);

loglik_b = zeros(2, N);
loglik_b(2, d==0) = log(seq_error);
loglik_b(2, dloc) = log(seq_error) - log(d(dloc));

log_pr_s = -Inf*ones(S*U, N);
log_pr_g = zeros(G, N);

for si = 1 : S
	
	% normal/tumour copy number
	cn_t = tumourState(si, 4);
	cn_n = normalState(si, 4);
		
	for ui = 1 : U

		% normal/tumour content
		ut = (1-u0)*(1-u(ui));
		un = u0 + (1-u0)*u(ui);
		
		% compute log-likelihood of observed read count
		m = un*read_depth*cn_n + ut*cn_t*read_depth;
		loglik_r(1, :) = log(1-seq_error) + logtpdf(dd, m, lambda_s^2, nu);
		loglik_r_sum = logsumexp(loglik_r, 1);
		
		% if this is a pure homozygous deletion
		if si == 1 & un > 0.01 
		
			ind = (si-1)*U + ui;
			log_pr_s(ind, :) = log(p_u(si, ui)) + loglik_r_sum + -1e9*(k > 0);
			
	 	else
		% otherwise
		
			% for each genotype, compute log-likehood of observed allele fraction
			for gi = 1 : G
	
				t_b = tumourState(si, gi);
				n_b = normalState(si, gi);
				gn = genotypeState(si, gi);	
		
				% compute expected allele fraction
				pp = ( un*n_b + ut*t_b )/( un*cn_n + ut*cn_t );
			
				% adjust for read errors 
				pp = pp*(1-read_error) + (1-pp)*read_error;

				% use binomial model to calculate alternate allele count probability
				loglik_b(1, dloc) = log(1-seq_error) + binomcoeff_dloc + k_dloc.*log(pp) + d_minus_k_dloc.*log(1-pp);
				loglik_b_sum = logsumexp(loglik_b, 1);

				% compute joint probability
				if options.paired
					log_pr_g(gi, :) = log_pr_gg(gn, :) + loglik_r_sum + loglik_b_sum;
				else
					log_pr_g(gi, :) = log(pg(si, gi)) + loglik_r_sum + loglik_b_sum;
				end

			end

			ind = (si-1)*U + ui;
			log_pr_s(ind, :) = log(p_u(si, ui)) + logsumexp(log_pr_g, 1);
			
		end
			
	end
	
end


