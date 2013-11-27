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
u = params.u;
p_u = params.p_u;
u0 = params.u0;
read_error = params.read_error;
read_depth = params.read_depth;
seq_error = params.seq_error;

tumourState = options.tumourState;
normalState = options.normalState;
chrRange = options.chrRange;

tumourState = options.tumourState;
normalState = options.normalState;
genotypeState = options.genotypeState;
pg = options.pg;

u0 = params.u0;
u = params.u;
p_u = params.p_u;

binomcoeff = gammaln(d+1) - gammaln(k+1) - gammaln(d-k+1);

log_pr_s = -Inf*ones(S, N);
log_pr_g = zeros(G, N);
log_pr_u = zeros(U, N);

d_max = max(d);

loglik_r = -Inf*ones(U, N, 2);
loglik_r(:, :, 2) = log(seq_error) - log(d_max);

dloc = find( d > 0 );
loglik_b = -Inf*ones(U, N, 2);
loglik_b(:, :, 2) = log(seq_error) - log(d_max);

for si = 1 : S
	
	cn_t = tumourState(si, 4);
	cn_n = normalState(si, 4);
		
	for gi = 1 : G
	
		t_b = tumourState(si, gi);
		n_b = normalState(si, gi);
		gn = genotypeState(si, gi);	
		
		ut = (1-u0)*(1-u);
		un = u0 + (1-u0)*u;
			
		if si > 1
			m = un*read_depth*cn_n + ut*cn_t*read_depth;
		else
			m = un*read_depth*cn_n;
		end
		
		loglik_r(:, :, 1) = log(1-seq_error) + logtpdf(d, m, lambda_s^2, nu);
		loglik_r_sum = logsumexp(loglik_r, 3);	
		
		if si > 1
			pp = ( un*n_b + ut*t_b )./( un*cn_n + ut*cn_t );
		else
			pp = (1/2)*ones(U, 1);
		end

		pp = pp*(1-read_error) + (1-pp)*read_error;
		
		loglik_b(:, dloc, 1) = log(1-seq_error) + repmat(binomcoeff(dloc), [U 1]);
		loglik_b(:, dloc, 1) = loglik_b(:, dloc, 1) + bsxfun(@times, k(dloc), log(pp)) + bsxfun(@times, d(dloc)-k(dloc), log(1-pp) );
		loglik_b_sum = logsumexp(loglik_b, 3);			
		
		log_pr_u = bsxfun(@plus, loglik_r_sum + loglik_b_sum, log(p_u(:, si)));	
		
%		for ui = 1 : U
%			
%			ut = (1-u0)*(1-u(ui));
%			un = u0 + (1-u0)*u(ui);
%		
%			if si > 1
%				m = un*read_depth*cn_n + ut*cn_t*read_depth;
%			else
%				m = un*read_depth*cn_n;
%			end
%		
%			loglik_r(1, :) = log(1-seq_error) + logtpdf(d, m, lambda_s^2, nu);
%			loglik_r_sum = logsumexp(loglik_r, 1);
%			
%			if si > 1
%				pp = ( un*n_b + ut*t_b )/( un*cn_n + ut*cn_t );
%			else
%				pp = 1/2;
%			end

%			pp = pp*(1-read_error) + (1-pp)*read_error;
%			
%			loglik_b(1, dloc) = log(1-seq_error) + binomcoeff(dloc) + k(dloc).*log(pp) + (d(dloc)-k(dloc)).*log(1-pp);
%			loglik_b_sum = logsumexp(loglik_b, 1);

%			log_pr_u(ui, :) = log(p_u(si, ui)) + loglik_r_sum + loglik_b_sum;
%		
%		end
		
		if options.paired
			log_pr_g(gi, :) = log_pr_gg(gn, :) + logsumexp(log_pr_u, 1);
		else
			log_pr_g(gi, :) = log(pg(si, gi)) + logsumexp(log_pr_u, 1);
		end
		
	end
	
	log_pr_s(si, :) = logsumexp(log_pr_g, 1);
	
end

