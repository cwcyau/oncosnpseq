/*
 * logsumexp2_api.c
 *
 * Code generation for function 'logsumexp2_api'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "logsumexp2_api.h"
#include "logsumexp2_emxutil.h"

/* Variable Definitions */
static emlrtRTEInfo h_emlrtRTEI = { 1, 1, "logsumexp2_api", "" };

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *a, const
  char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(emxArray_real_T *u);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  c_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv3[2];
  boolean_T bv0[2];
  int32_T i1;
  static const boolean_T bv1[2] = { FALSE, TRUE };

  int32_T iv4[2];
  for (i1 = 0; i1 < 2; i1++) {
    iv3[i1] = 1 + 9999999 * i1;
    bv0[i1] = bv1[i1];
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", FALSE, 2U, iv3, bv0, iv4);
  ret->size[0] = iv4[0];
  ret->size[1] = iv4[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *a, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(sp, emlrtAlias(a), &thisId, y);
  emlrtDestroyArray(&a);
}

static const mxArray *emlrt_marshallOut(emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv2[2] = { 0, 0 };

  const mxArray *m2;
  y = NULL;
  m2 = mxCreateNumericArray(2, (int32_T *)&iv2, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m2, (void *)u->data);
  mxSetDimensions((mxArray *)m2, u->size, 2);
  emlrtAssign(&y, m2);
  return y;
}

void logsumexp2_api(emlrtStack *sp, const mxArray * const prhs[1], const mxArray
                    *plhs[1])
{
  emxArray_real_T *a;
  emxArray_real_T *z;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &a, 2, &h_emlrtRTEI, TRUE);
  emxInit_real_T(sp, &z, 2, &h_emlrtRTEI, TRUE);

  /* Marshall function inputs */
  emlrt_marshallIn(sp, emlrtAlias(prhs[0]), "a", a);

  /* Invoke the target function */
  logsumexp2(sp, a, z);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(z);
  z->canFreeData = FALSE;
  emxFree_real_T(&z);
  a->canFreeData = FALSE;
  emxFree_real_T(&a);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (logsumexp2_api.c) */
