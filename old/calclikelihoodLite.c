#include "math.h"


double logtpdf( double x, double m, double S, double v ) {

	double z = x-m;
	double d = (1/S)*z*z;
	double y = -0.5*log(S) - 0.5*(v+1)*log(1 + d/v);
	return y;
	
};

double logsumexp(double nums[], int ct) {
  
  	double max_exp = nums[0], sum = 0.0;
  	int i;

  	for (i = 1 ; i < ct ; i++)
    	if (nums[i] > max_exp)
      		max_exp = nums[i];

  	for (i = 0; i < ct ; i++)
    	sum += exp(nums[i] - max_exp);

  	return log(sum) + max_exp;
}


double* calc_log_obslik(double k, double d, double* log_pr_gg, double d_max, double u0, double read_depth, double* u, double seq_error, double read_error, double lambda_s, double nu, double* tumourState, double* normalState, double* genotypeState, int S, int U, int G) { 

	int si, ui, gi;
	
	double cn_t, cn_n, ut, un, m;
	
	double log_pr_s[S*U];
	
	double loglik_r[2], loglik_b[2];
	double loglik_r_sum, loglik_b_sum;		
	
	double binomcoeff;
	double pp;
	double log_pr_g[G];
	
	double t_b, n_b;
	int gn;

	for ( si = 0; si < S; si++ ) {
	
		cn_t = tumourState[4*S + si];
		cn_n = normalState[4*S + si];
	
		for ( ui = 0; ui < U; ui++ ) {
		
			ut = (1-u0)*(1-u[ui]);
			un = u0 + (1-u0)*u[ui];

			if ( si > 1 ) {
				m = un*read_depth*cn_n + ut*cn_t*read_depth;
			} else {
				m = un*read_depth*cn_n;
			}
			
			loglik_r[0] = log(1-seq_error) + logtpdf(d, m, lambda_s*lambda_s, nu);
			loglik_r[1] = log(seq_error) - log(d_max);
			loglik_r_sum = logsumexp( loglik_r, 2);

			binomcoeff = lgamma(d+1) - lgamma(k+1) - lgamma(d-k+1);
	
			for ( gi = 0; gi < G; gi++ ) {
				
				t_b = tumourState[ (gi-1)*S + si ];
				n_b = normalState[ (gi-1)*S + si ];
				gn = (int )genotypeState[ (gi-1)*S + si ];	
		
				if ( si > 1 ) {
					pp = ( un*n_b + ut*t_b )/( un*cn_n + ut*cn_t );
				} 
				else {
					pp = 0.5;
				}
				pp = pp*(1-read_error) + (1-pp)*read_error;
			
		
				loglik_b[0] = log(1-seq_error) + binomcoeff + k*log(pp) + (d-k)*log(1-pp);
				loglik_b[1] = log(seq_error) - log(1);
				loglik_b_sum = logsumexp(loglik_b, 2);

				log_pr_g[gi] = log_pr_gg[gn] + loglik_r_sum + loglik_b_sum;
			}
	
			log_pr_s[(si-1)*U + ui] = logsumexp(log_pr_g, G);
			if ( ui > 1 ) {
				log_pr_s[(si-1)*U + ui] -= 1;
			}
		
		}
	
	}

}

int main () {
	return 0;
}
