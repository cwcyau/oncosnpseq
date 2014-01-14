/*
 * eml_error.c
 *
 * Code generation for function 'eml_error'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "eml_error.h"

/* Variable Definitions */
static emlrtRTEInfo i_emlrtRTEI = { 20, 5, "eml_error",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_error.m" };

/* Function Definitions */
void eml_error(const emlrtStack *sp)
{
  static char_T cv4[3][1] = { { 'l' }, { 'o' }, { 'g' } };

  emlrtErrorWithMessageIdR2012b(sp, &i_emlrtRTEI,
    "Coder:toolbox:ElFunDomainError", 3, 4, 3, cv4);
}

/* End of code generation (eml_error.c) */
