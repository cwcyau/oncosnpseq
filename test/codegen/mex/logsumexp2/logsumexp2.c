/*
 * logsumexp2.c
 *
 * Code generation for function 'logsumexp2'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "exp.h"
#include "logsumexp2_emxutil.h"
#include "isfinite.h"
#include "log.h"
#include "bsxfun.h"
#include "logsumexp2_mexutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 5, "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m" };

static emlrtRSInfo b_emlrtRSI = { 10, "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m" };

static emlrtRSInfo c_emlrtRSI = { 12, "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m" };

static emlrtRSInfo d_emlrtRSI = { 14, "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m" };

static emlrtRSInfo e_emlrtRSI = { 18, "max",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/datafun/max.m" };

static emlrtRSInfo f_emlrtRSI = { 15, "eml_min_or_max",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo rb_emlrtRSI = { 11, "eml_li_find",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRSInfo sb_emlrtRSI = { 29, "eml_li_find",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtMCInfo g_emlrtMCI = { 14, 5, "eml_li_find",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRTEInfo emlrtRTEI = { 83, 1, "eml_min_or_max",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRTEInfo b_emlrtRTEI = { 1, 14, "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m" };

static emlrtRTEInfo g_emlrtRTEI = { 17, 9, "eml_li_find",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtECInfo emlrtECI = { 2, 12, 5, "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m" };

static emlrtECInfo b_emlrtECI = { -1, 14, 1, "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m" };

static emlrtBCInfo emlrtBCI = { -1, -1, 14, 22, "amax", "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m", 0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 14, 1, "z", "logsumexp2",
  "/Net/fs1/home/cyau/projects/oncosnpseq/logsumexp2.m", 0 };

static emlrtDCInfo emlrtDCI = { 17, 37, "eml_li_find",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_li_find.m", 4 };

static emlrtRSInfo ub_emlrtRSI = { 14, "eml_li_find",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

/* Function Declarations */
static int32_T compute_nones(const emxArray_boolean_T *x, int32_T n);
static void eml_li_find(const emlrtStack *sp, const emxArray_boolean_T *x,
  emxArray_int32_T *y);

/* Function Definitions */
static int32_T compute_nones(const emxArray_boolean_T *x, int32_T n)
{
  int32_T k;
  int32_T i;
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  return k;
}

static void eml_li_find(const emlrtStack *sp, const emxArray_boolean_T *x,
  emxArray_int32_T *y)
{
  int32_T k;
  const mxArray *b_y;
  const mxArray *m1;
  int32_T i;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &rb_emlrtRSI;
  k = compute_nones(x, x->size[1]);
  if (k <= x->size[1]) {
  } else {
    b_y = NULL;
    m1 = mxCreateString("Assertion failed.");
    emlrtAssign(&b_y, m1);
    st.site = &ub_emlrtRSI;
    error(&st, b_y, &g_emlrtMCI);
  }

  emlrtNonNegativeCheckFastR2012b(k, &emlrtDCI, sp);
  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = k;
  emxEnsureCapacity(sp, (emxArray__common *)y, i, (int32_T)sizeof(int32_T),
                    &g_emlrtRTEI);
  k = 0;
  for (i = 1; i <= x->size[1]; i++) {
    if (x->data[i - 1]) {
      y->data[k] = i;
      st.site = &sb_emlrtRSI;
      k++;
    }
  }
}

void logsumexp2(const emlrtStack *sp, const emxArray_real_T *a, emxArray_real_T *
                z)
{
  int32_T outsz[2];
  int32_T iy;
  emxArray_real_T *extremum;
  int32_T ix;
  int32_T i;
  int32_T b_z[2];
  emxArray_boolean_T *r0;
  emxArray_boolean_T *r1;
  emxArray_int32_T *r2;
  emxArray_boolean_T *r3;
  emxArray_real_T *r4;
  emxArray_int32_T *r5;
  int32_T i0;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  st.site = &emlrtRSI;
  b_st.site = &e_emlrtRSI;
  c_st.site = &f_emlrtRSI;
  for (iy = 0; iy < 2; iy++) {
    outsz[iy] = a->size[iy];
  }

  emxInit_real_T(&c_st, &extremum, 2, &b_emlrtRTEI, TRUE);
  iy = extremum->size[0] * extremum->size[1];
  extremum->size[0] = 1;
  extremum->size[1] = outsz[1];
  emxEnsureCapacity(&c_st, (emxArray__common *)extremum, iy, (int32_T)sizeof
                    (real_T), &emlrtRTEI);
  ix = -1;
  iy = -1;
  for (i = 1; i <= a->size[1]; i++) {
    ix++;
    iy++;
    extremum->data[iy] = a->data[ix];
  }

  st.site = &b_emlrtRSI;
  bsxfun(&st, a, extremum, z);
  st.site = &c_emlrtRSI;
  b_exp(z);
  st.site = &c_emlrtRSI;
  st.site = &c_emlrtRSI;
  b_log(&st, z);
  for (iy = 0; iy < 2; iy++) {
    outsz[iy] = extremum->size[iy];
  }

  for (iy = 0; iy < 2; iy++) {
    b_z[iy] = z->size[iy];
  }

  emlrtSizeEqCheck2DFastR2012b(outsz, b_z, &emlrtECI, sp);
  iy = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = extremum->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)z, iy, (int32_T)sizeof(real_T),
                    &b_emlrtRTEI);
  i = extremum->size[0] * extremum->size[1];
  for (iy = 0; iy < i; iy++) {
    z->data[iy] += extremum->data[iy];
  }

  emxInit_boolean_T(sp, &r0, 2, &b_emlrtRTEI, TRUE);
  emxInit_boolean_T(sp, &r1, 2, &b_emlrtRTEI, TRUE);
  st.site = &d_emlrtRSI;
  b_isfinite(&st, extremum, r0);
  iy = r1->size[0] * r1->size[1];
  r1->size[0] = 1;
  r1->size[1] = r0->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)r1, iy, (int32_T)sizeof(boolean_T),
                    &b_emlrtRTEI);
  i = r0->size[0] * r0->size[1];
  for (iy = 0; iy < i; iy++) {
    r1->data[iy] = !r0->data[iy];
  }

  emxInit_int32_T(sp, &r2, 2, &b_emlrtRTEI, TRUE);
  emxInit_boolean_T(sp, &r3, 2, &b_emlrtRTEI, TRUE);
  st.site = &d_emlrtRSI;
  eml_li_find(&st, r1, r2);
  st.site = &d_emlrtRSI;
  b_isfinite(&st, extremum, r0);
  iy = r3->size[0] * r3->size[1];
  r3->size[0] = 1;
  r3->size[1] = r0->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)r3, iy, (int32_T)sizeof(boolean_T),
                    &b_emlrtRTEI);
  i = r0->size[0] * r0->size[1];
  emxFree_boolean_T(&r1);
  for (iy = 0; iy < i; iy++) {
    r3->data[iy] = !r0->data[iy];
  }

  emxFree_boolean_T(&r0);
  emxInit_real_T(sp, &r4, 2, &b_emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r5, 2, &b_emlrtRTEI, TRUE);
  st.site = &d_emlrtRSI;
  eml_li_find(&st, r3, r5);
  iy = r4->size[0] * r4->size[1];
  r4->size[0] = 1;
  r4->size[1] = r5->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)r4, iy, (int32_T)sizeof(real_T),
                    &b_emlrtRTEI);
  i = r5->size[0] * r5->size[1];
  emxFree_boolean_T(&r3);
  for (iy = 0; iy < i; iy++) {
    i0 = extremum->size[1];
    ix = r5->data[iy];
    r4->data[iy] = extremum->data[emlrtDynamicBoundsCheckFastR2012b(ix, 1, i0,
      &emlrtBCI, sp) - 1];
  }

  emxFree_int32_T(&r5);
  emxFree_real_T(&extremum);
  iy = r2->size[1];
  i0 = r4->size[1];
  emlrtSizeEqCheck1DFastR2012b(iy, i0, &b_emlrtECI, sp);
  ix = z->size[1];
  i = r4->size[0] * r4->size[1];
  for (iy = 0; iy < i; iy++) {
    i0 = r2->data[iy];
    z->data[emlrtDynamicBoundsCheckFastR2012b(i0, 1, ix, &b_emlrtBCI, sp) - 1] =
      r4->data[iy];
  }

  emxFree_real_T(&r4);
  emxFree_int32_T(&r2);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (logsumexp2.c) */
