/*
 * isfinite.c
 *
 * Code generation for function 'isfinite'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "isfinite.h"
#include "logsumexp2_emxutil.h"

/* Variable Definitions */
static emlrtRTEInfo f_emlrtRTEI = { 1, 14, "isfinite",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/isfinite.m" };

/* Function Definitions */
void b_isfinite(const emlrtStack *sp, const emxArray_real_T *x,
                emxArray_boolean_T *b)
{
  int32_T b_b;
  int32_T loop_ub;
  emxArray_boolean_T *r6;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  b_b = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = x->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)b, b_b, (int32_T)sizeof(boolean_T),
                    &f_emlrtRTEI);
  loop_ub = x->size[0] * x->size[1];
  for (b_b = 0; b_b < loop_ub; b_b++) {
    b->data[b_b] = muDoubleScalarIsInf(x->data[b_b]);
  }

  emxInit_boolean_T(sp, &r6, 2, &f_emlrtRTEI, TRUE);
  b_b = r6->size[0] * r6->size[1];
  r6->size[0] = 1;
  r6->size[1] = x->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)r6, b_b, (int32_T)sizeof(boolean_T),
                    &f_emlrtRTEI);
  loop_ub = x->size[0] * x->size[1];
  for (b_b = 0; b_b < loop_ub; b_b++) {
    r6->data[b_b] = muDoubleScalarIsNaN(x->data[b_b]);
  }

  b_b = b->size[0] * b->size[1];
  b->size[0] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)b, b_b, (int32_T)sizeof(boolean_T),
                    &f_emlrtRTEI);
  b_b = b->size[0];
  loop_ub = b->size[1];
  loop_ub *= b_b;
  for (b_b = 0; b_b < loop_ub; b_b++) {
    b->data[b_b] = ((!b->data[b_b]) && (!r6->data[b_b]));
  }

  emxFree_boolean_T(&r6);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (isfinite.c) */
