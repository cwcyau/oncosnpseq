/*
 * exp.c
 *
 * Code generation for function 'exp'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "exp.h"

/* Function Definitions */
void b_exp(emxArray_real_T *x)
{
  int32_T i2;
  int32_T k;
  i2 = x->size[1];
  for (k = 0; k < i2; k++) {
    x->data[k] = muDoubleScalarExp(x->data[k]);
  }
}

/* End of code generation (exp.c) */
