%
% OncoSEQ
%
% Written by: Christopher Yau
% 
% 19 March 2012
%
% Copyright (C) 2012 Imperial College London
%
function u_vec = calclikelihoodU(k, d, dd, log_pr_gg, x, params, options)

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

dd_max = max(dd);

log_pr_u = -Inf*ones(U, N);

for ui = 1 : U

	ut = (1-u0)*(1-u(ui));
	un = u0 + (1-u0)*u(ui);
	
	for si = 1 : S
		
		loc = find( x == si & d > 0 );
		n_loc = length(loc);
		if n_loc == 0
			continue;
		end
	
		cn_t = tumourState(si, 4);
		cn_n = normalState(si, 4);
	
		if si > 1
			m = un*read_depth*cn_n + ut*cn_t*read_depth;
		else
			m = un*read_depth*cn_n;
		end
		
		loglik_r = -Inf*ones(2, N);
		loglik_r(1, loc) = log(1-seq_error) + logtpdf(dd(loc), m, lambda_s^2, nu);
		loglik_r(2, loc) = log(seq_error) - log(dd_max);
		loglik_r_sum = logsumexp(loglik_r, 1);
		
		log_pr_g = -Inf + zeros(G, N);
		for gi = 1 : G
			
			t_b = tumourState(si, gi);
			n_b = normalState(si, gi);
			gn = genotypeState(si, gi);

			if si > 1
				pp = ( un*n_b + ut*t_b )/( un*cn_n + ut*cn_t );
			else
				pp = 1/2;
			end
			pp = pp*(1-read_error) + (1-pp)*read_error;
			
			loglik_b = -Inf*ones(2, N);
			loglik_b(1, loc) = log(1-seq_error) + binomcoeff(loc) + k(loc).*log(pp) + (d(loc)-k(loc)).*log(1-pp);
			loglik_b(2, loc) = log(seq_error) - log(d(loc));
			loglik_b_sum = logsumexp(loglik_b(:, loc), 1);
	
			if options.paired
				log_pr_g(gi, loc) = log_pr_gg(gn, loc) + loglik_b_sum;
			else
				log_pr_g(gi, loc) = log(pg(si, gi)) + loglik_b_sum;
			end
			
		end
		
		log_pr_u(ui, loc) = log(p_u(ui, si)) + loglik_r_sum(loc) + logsumexp(log_pr_g(:, loc), 1);
	
	end
	
end

u_vec = zeros(1, N);

startInd = 1;
for i = 2 : N
	if x(i) ~= x(i-1) | i == N
	
		seglen = length(startInd:i-1);
	
		pu = sum(log_pr_u(:, startInd:i-1), 2);
		
		[ mxval, mxloc ] = max(pu);
		u_vec(startInd:i-1) = u0 + (1-u0)*u(mxloc);
		
		startInd = i;
	end
end

