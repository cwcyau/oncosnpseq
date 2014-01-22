#include "mex.h"
#include "matrix.h"
#include "math.h"
    
/* If you are using a compiler that equates NaN to zero, you must
* compile this example using the flag -DNAN_EQUALS_ZERO. For 
* example:
*
*     mex -DNAN_EQUALS_ZERO findnz.c  
*
* This will correctly define the IsNonZero macro for your
compiler. */

#if NAN_EQUALS_ZERO
    #define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
    #define IsNonZero(d) ((d) != 0.0)
#endif

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


double calc_log_obslik(double k, double d, double* log_pr_gg, double ut, double un, double cn_t, double cn_n, double d_max, double read_depth, double seq_error, double read_error, double m, double lambda_s, double nu, double* pp, double* t_b, double* n_b, int* gn, int G) { 

	int gi;
	
	double log_pr_s;
	
	double loglik_r[2], loglik_b[2];
	double loglik_r_sum, loglik_b_sum;		
	
	double binomcoeff;
	double log_pr_g[G];


	loglik_r[0] = log(1-seq_error) + logtpdf(d, m, lambda_s*lambda_s, nu);
	loglik_r[1] = log(seq_error) - log(d_max);
	loglik_r_sum = logsumexp( loglik_r, 2);

	binomcoeff = lgamma(d+1) - lgamma(k+1) - lgamma(d-k+1);

	for ( gi = 0; gi < G; gi++ ) {

		loglik_b[0] = log(1-seq_error) + binomcoeff + k*log(pp[gi]) + (d-k)*log(1-pp[gi]);
		loglik_b[1] = log(seq_error) - log(1);
		loglik_b_sum = logsumexp(loglik_b, 2);

		log_pr_g[gi] = log_pr_gg[gn[gi]] + loglik_r_sum + loglik_b_sum;
	}

	log_pr_s = logsumexp(log_pr_g, G);

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* Declare variables. */ 
    int i, j, t, g;
	mxArray *tmpVec;
	double tmpMax;
	int tmpLoc, ind;
    
    /* Get the data. */
	double* logNu = (double *)mxGetPr(prhs[0]);
	double* logtransMat = (double *)mxGetPr(prhs[1]);
	int* arrayind = (int *)mxGetPr(prhs[1]);
    double* k = (double *)mxGetPr(prhs[5]);
    double* d = (double *)mxGetPr(prhs[6]);
	double* log_pr_gg = (double *)mxGetPr(prhs[6]);
	double u0 = (double )mxGetScalar(prhs[6]);
	double read_depth = (double )mxGetScalar(prhs[6]);
    double* u = (double *)mxGetPr(prhs[6]);
	double seq_error = (double )mxGetScalar(prhs[6]);
	double read_error = (double )mxGetScalar(prhs[6]);
	double lambda_s = (double )mxGetScalar(prhs[6]);
	double nu = (double )mxGetScalar(prhs[6]);
	double* tumourState = (double *)mxGetPr(prhs[6]); 
	double* normalState = (double *)mxGetPr(prhs[6]);
	double* genotypeState = (double *)mxGetPr(prhs[6]);
    int S = (int )mxGetScalar(prhs[3]);
    int U = (int )mxGetScalar(prhs[4]);
    int G = (int )mxGetScalar(prhs[4]);    
	int T = (int )mxGetScalar(prhs[4]);

	int si, sj, ui, uj;

	double logobslik;
    
	plhs[0] = mxCreateDoubleMatrix(1, T, mxREAL);
    double* path = mxGetPr(plhs[0]);	

	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* loglik = mxGetPr(plhs[1]);	

	double phi[S];
	double phi_old[S]; 

	mxArray* delta_mat = mxCreateNumericMatrix(S, T, mxINT32_CLASS, mxREAL);
    int* delta = (int *) mxGetPr(delta_mat);

	double tmp[S];
	double p;

	double ut, un, cn_t, cn_n, m;

	double log_pr_gg_t[G], t_b[G], n_b[G], pp[G];
	int gn[G];

	double dmax = d[0];
	for ( t = 1; t < T; t++ ) {
		if ( d[t] > dmax ) {
			dmax = d[t];
		}
	}

    for ( i = 0; i < S; i++ ) {
		phi[i] = 0;
		phi_old[i] = 0;
	}
    for ( i = 0; i < S; i++ ) {

		si = arrayind[i];
		ui = arrayind[S+i];


		cn_t = tumourState[4*S + si];
		cn_n = normalState[4*S + si];

		ut = (1-u0)*(1-u[ui]);
		un = u0 + (1-u0)*u[ui];

		if ( si > 1 ) {
			m = un*read_depth*cn_n + ut*cn_t*read_depth;
		} else {
			m = un*read_depth*cn_n;
		}

		for ( g = 0; g < G; g++ ) {
			log_pr_gg_t[g] = log_pr_gg[G*t + g];
			t_b[g] = tumourState[ g*S + si ];
			n_b[g] = normalState[ g*S + si ];
			gn[g] = (int )genotypeState[ g*S + si ];	
			if ( si > 1 ) {
				pp[g] = ( un*n_b[g] + ut*t_b[g] )/( un*cn_n + ut*cn_t );
			} 
			else {
				pp[g] = 0.5;
			}
			pp[g] = pp[g]*(1-read_error) + (1-pp[g])*read_error;
		}
		logobslik = calc_log_obslik(k[t], d[t], log_pr_gg_t, ut, un, cn_t, cn_n, dmax, read_depth, seq_error, read_error, m, lambda_s, nu, pp, t_b, n_b, gn, G);
		if ( ui > 1 ) {
			logobslik -= 1;
		}

        phi[i] = logNu[i] + logobslik;
    }

    for ( t = 1; t < T; t++ ) {

		for ( i = 0; i < S; i++ ) {
			phi_old[i] = phi[i];
		}	

		for ( j = 0; j < S; j++ ) {

			sj = arrayind[j];
			uj = arrayind[S+j];

			cn_t = tumourState[4*S + sj];
			cn_n = normalState[4*S + sj];

			ut = (1-u0)*(1-u[uj]);
			un = u0 + (1-u0)*u[uj];

			if ( sj > 1 ) {
				m = un*read_depth*cn_n + ut*cn_t*read_depth;
			} else {
				m = un*read_depth*cn_n;
			}

			for ( g = 0; g < G; g++ ) {
				log_pr_gg_t[g] = log_pr_gg[G*t + g];
				t_b[g] = tumourState[ g*S + sj ];
				n_b[g] = normalState[ g*S + sj ];
				gn[g] = (int )genotypeState[ g*S + sj ];	
				if ( sj > 1 ) {
					pp[g] = ( un*n_b[g] + ut*t_b[g] )/( un*cn_n + ut*cn_t );
				} 
				else {
					pp[g] = 0.5;
				}
				pp[g] = pp[g]*(1-read_error) + (1-pp[g])*read_error;
			}
			logobslik = calc_log_obslik(k[t], d[t], log_pr_gg_t, ut, un, cn_t, cn_n, dmax, read_depth, seq_error, read_error, m, lambda_s, nu, pp, t_b, n_b, gn, G);
			if ( uj > 1 ) {
				logobslik -= 1;
			}

			for ( i = 0; i < S; i++ ) {

				si = arrayind[i];
				ui = arrayind[S+i];

				p = logobslik + logtransMat[j*S + i];
				tmp[i] = phi_old[i] + p;
			}

			tmpMax = tmp[0];
			tmpLoc = 0;
			for ( i = 1; i < S; i++ ) {
				if ( tmp[i] > tmpMax ) {
					tmpMax = tmp[i];
					tmpLoc = i;
				}
			}

			phi[j] = tmpMax;
			delta[S*t + j] = tmpLoc;

		}
		
    }

	tmpMax = phi[0];
	tmpLoc = 0;
	for ( i = 1; i < S; i++ ) {
		if ( phi[i] > tmpMax ) {
			tmpMax = phi[i];
			tmpLoc = i;
		}
	}
	loglik[0] = tmpMax;

	path[T-1] = tmpLoc;
	for ( t = T-2; t > -1; t-- ) {
		ind = (t+1)*S + path[t+1];
		path[t] = delta[ind];
	}

}


