/*
 * log.c
 *
 * Code generation for function 'log'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "log.h"
#include "eml_error.h"

/* Variable Definitions */
static emlrtRSInfo lb_emlrtRSI = { 14, "log",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elfun/log.m" };

/* Function Definitions */
void b_log(const emlrtStack *sp, emxArray_real_T *x)
{
  int32_T k;
  int32_T i3;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  for (k = 0; k < x->size[1]; k++) {
    if (x->data[k] < 0.0) {
      st.site = &lb_emlrtRSI;
      eml_error(&st);
    }
  }

  i3 = x->size[1];
  for (k = 0; k < i3; k++) {
    x->data[k] = muDoubleScalarLog(x->data[k]);
  }
}

/* End of code generation (log.c) */
