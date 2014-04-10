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
    int nStates, nSnps;
    double *obsMatrix, *cost; 
	mxArray *tmpVec;
	double *tmp, *delta, *phi, *path;
	double tmpMax;
	int tmpLoc, ind;
    
    /* Get the data. */
    obsMatrix           = (double *)mxGetPr(prhs[0]);
    nStates             = (int )mxGetScalar(prhs[1]);
    nSnps               = (int )mxGetScalar(prhs[2]);
	cost 				= (double *)mxGetPr(prhs[3]);
    
	plhs[0] = mxCreateDoubleMatrix(1, nSnps, mxREAL);
    path = mxGetPr(plhs[0]);	
	
    plhs[1] = mxCreateDoubleMatrix(nStates, nSnps, mxREAL);
    phi = mxGetPr(plhs[1]);
	
    plhs[2] = mxCreateDoubleMatrix(nStates, nSnps, mxREAL);
    delta = mxGetPr(plhs[2]);

	tmpVec = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    tmp = mxGetPr(tmpVec);

	for ( k = 0; k < nSnps; k++ ) {
	    for ( i = 0; i < nStates; i++ ) {
			phi[nStates*k + i] = 0;
		}
	}
    for ( i = 0; i < nStates; i++ ) {
        phi[i] = obsMatrix[i];
    }

    for ( k = 1; k < nSnps; k++ ) {

		for ( j = 0; j < nStates; j++ ) {
		
			for ( i = 0; i < nStates; i++ ) {
				tmp[i] = phi[nStates*(k-1) + i] + obsMatrix[nStates*k + j] - cost[i*nStates + j];
			}

			tmpMax = tmp[0];
			tmpLoc = 0;
			for ( m = 1; m < nStates; m++ ) {
				if ( tmp[m] > tmpMax ) {
					tmpMax = tmp[m];
					tmpLoc = m;
				}
			}

			phi[nStates*k + j] = tmpMax;
			delta[nStates*k + j] = tmpLoc;

		}
		
    }

	tmpMax = phi[(nSnps-1)*nStates + 0];
	tmpLoc = 0;
	for ( m = 1; m < nStates; m++ ) {
		if ( phi[(nSnps-1)*nStates + m] > tmpMax ) {
			tmpMax = phi[(nSnps-1)*nStates + m];
			tmpLoc = m;
		}
	}

	path[nSnps-1] = tmpLoc;
	for ( k = nSnps-2; k > -1; k-- ) {
		ind = (k+1)*nStates + path[k+1];
		path[k] = delta[ind];
	}

}



