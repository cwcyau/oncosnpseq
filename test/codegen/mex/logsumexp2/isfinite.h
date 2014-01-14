/*
 * isfinite.h
 *
 * Code generation for function 'isfinite'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

#ifndef __ISFINITE_H__
#define __ISFINITE_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "logsumexp2_types.h"

/* Function Declarations */
extern void b_isfinite(const emlrtStack *sp, const emxArray_real_T *x, emxArray_boolean_T *b);
#endif
/* End of code generation (isfinite.h) */
