/*
 * logsumexp2_mexutil.c
 *
 * Code generation for function 'logsumexp2_mexutil'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "logsumexp2_mexutil.h"

/* Function Definitions */
void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "error", TRUE, location);
}

/* End of code generation (logsumexp2_mexutil.c) */
