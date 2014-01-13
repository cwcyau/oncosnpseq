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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* Declare variables. */ 
    int i, j, k, m;
	mxArray *tmpVec;
	double tmpMax;
	int tmpLoc, ind;
    
    /* Get the data. */
	double* logNu = (double *)mxGetPr(prhs[0]);
	double* logtransMat = (double *)mxGetPr(prhs[1]);
    double* logobsMatrix = (double *)mxGetPr(prhs[2]);
    int nStates = (int )mxGetScalar(prhs[3]);
    int nSnps = (int )mxGetScalar(prhs[4]);
    
	plhs[0] = mxCreateDoubleMatrix(1, nSnps, mxREAL);
    double* path = mxGetPr(plhs[0]);	

	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* loglik = mxGetPr(plhs[1]);	

	double phi[nStates];
	double phi_old[nStates]; 

	mxArray* delta_mat = mxCreateNumericMatrix(nStates, nSnps, mxINT32_CLASS, mxREAL);
    int* delta = (int *) mxGetPr(delta_mat);

	double tmp[nStates];
	double p;

    for ( i = 0; i < nStates; i++ ) {
		phi[i] = 0;
		phi_old[i] = 0;
	}
    for ( i = 0; i < nStates; i++ ) {
        phi[i] = logNu[i] + logobsMatrix[i];
    }

    for ( k = 1; k < nSnps; k++ ) {

		for ( i = 0; i < nStates; i++ ) {
			phi_old[i] = phi[i];
		}	

		for ( j = 0; j < nStates; j++ ) {
		
			for ( i = 0; i < nStates; i++ ) {
				p = logobsMatrix[nStates*k + j] + logtransMat[j*nStates + i];
				tmp[i] = phi_old[i] + p;
			}

			tmpMax = tmp[0];
			tmpLoc = 0;
			for ( m = 1; m < nStates; m++ ) {
				if ( tmp[m] > tmpMax ) {
					tmpMax = tmp[m];
					tmpLoc = m;
				}
			}

			phi[j] = tmpMax;
			delta[nStates*k + j] = tmpLoc;

		}
		
    }

	tmpMax = phi[0];
	tmpLoc = 0;
	for ( m = 1; m < nStates; m++ ) {
		if ( phi[m] > tmpMax ) {
			tmpMax = phi[m];
			tmpLoc = m;
		}
	}
	loglik[0] = tmpMax;

	path[nSnps-1] = tmpLoc;
	for ( k = nSnps-2; k > -1; k-- ) {
		ind = (k+1)*nStates + path[k+1];
		path[k] = delta[ind];
	}

}


