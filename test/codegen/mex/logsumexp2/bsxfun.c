/*
 * bsxfun.c
 *
 * Code generation for function 'bsxfun'
 *
 * C source code generated on: Tue Jan 14 10:28:31 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "logsumexp2.h"
#include "bsxfun.h"
#include "logsumexp2_emxutil.h"
#include "logsumexp2_mexutil.h"

/* Variable Definitions */
static emlrtRSInfo n_emlrtRSI = { 21, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo o_emlrtRSI = { 23, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo p_emlrtRSI = { 75, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo q_emlrtRSI = { 78, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo r_emlrtRSI = { 82, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo s_emlrtRSI = { 88, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo t_emlrtRSI = { 96, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo u_emlrtRSI = { 97, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo v_emlrtRSI = { 101, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo w_emlrtRSI = { 102, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo x_emlrtRSI = { 104, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo y_emlrtRSI = { 107, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo ab_emlrtRSI = { 108, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo bb_emlrtRSI = { 109, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo cb_emlrtRSI = { 110, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo c_emlrtMCI = { 22, 5, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo d_emlrtMCI = { 21, 15, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo e_emlrtMCI = { 24, 5, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo f_emlrtMCI = { 23, 15, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRTEInfo c_emlrtRTEI = { 41, 1, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo vb_emlrtRSI = { 24, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo wb_emlrtRSI = { 22, "bsxfun",
  "/opt/well/matlab/matlab-2013b/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

/* Function Declarations */
static boolean_T bsxfun_compatible(const emxArray_real_T *a, const
  emxArray_real_T *b);
static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location);
static boolean_T no_dynamic_expansion(const emxArray_real_T *a, const
  emxArray_real_T *b);

/* Function Definitions */
static boolean_T bsxfun_compatible(const emxArray_real_T *a, const
  emxArray_real_T *b)
{
  boolean_T p;
  if ((a->size[1] != 1) && (b->size[1] != 1) && (a->size[1] != b->size[1])) {
    p = FALSE;
  } else {
    p = TRUE;
  }

  return p;
}

static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m3;
  pArray = b;
  return emlrtCallMATLABR2012b(sp, 1, &m3, 1, &pArray, "message", TRUE, location);
}

static boolean_T no_dynamic_expansion(const emxArray_real_T *a, const
  emxArray_real_T *b)
{
  boolean_T p;
  if (a->size[1] != b->size[1]) {
    p = FALSE;
  } else {
    p = TRUE;
  }

  return p;
}

void bsxfun(const emlrtStack *sp, const emxArray_real_T *a, const
            emxArray_real_T *b, emxArray_real_T *c)
{
  const mxArray *y;
  static const int32_T iv0[2] = { 1, 44 };

  const mxArray *m0;
  char_T cv0[44];
  int32_T i;
  static const char_T cv1[44] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'b', 's', 'x', 'f', 'u', 'n', '_', 'a', 'r', 'r', 'a',
    'y', 'D', 'i', 'm', 'e', 'n', 's', 'i', 'o', 'n', 's', 'M', 'u', 's', 't',
    'M', 'a', 't', 'c', 'h' };

  const mxArray *b_y;
  static const int32_T iv1[2] = { 1, 37 };

  char_T cv2[37];
  static const char_T cv3[37] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'b', 's', 'x', 'f', 'u', 'n', '_', 'd', 'y', 'n',
    'a', 'm', 'i', 'c', 'E', 'x', 'p', 'a', 'n', 's', 'i', 'o', 'n' };

  int32_T asub;
  int32_T bsub;
  int32_T ak;
  int32_T bk;
  int32_T ck;
  int32_T exitg1;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &n_emlrtRSI;
  b_st.prev = sp;
  b_st.tls = sp->tls;
  if (bsxfun_compatible(a, b)) {
  } else {
    y = NULL;
    m0 = mxCreateCharArray(2, iv0);
    for (i = 0; i < 44; i++) {
      cv0[i] = cv1[i];
    }

    emlrtInitCharArrayR2013a(sp, 44, m0, cv0);
    emlrtAssign(&y, m0);
    st.site = &n_emlrtRSI;
    b_st.site = &wb_emlrtRSI;
    error(&st, message(&b_st, y, &c_emlrtMCI), &d_emlrtMCI);
  }

  st.site = &o_emlrtRSI;
  if (no_dynamic_expansion(a, b)) {
  } else {
    b_y = NULL;
    m0 = mxCreateCharArray(2, iv1);
    for (i = 0; i < 37; i++) {
      cv2[i] = cv3[i];
    }

    emlrtInitCharArrayR2013a(sp, 37, m0, cv2);
    emlrtAssign(&b_y, m0);
    st.site = &o_emlrtRSI;
    b_st.site = &vb_emlrtRSI;
    error(&st, message(&b_st, b_y, &e_emlrtMCI), &f_emlrtMCI);
  }

  i = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = a->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)c, i, (int32_T)sizeof(real_T),
                    &c_emlrtRTEI);
  i = a->size[1];
  if (i == 0) {
  } else {
    asub = 1;
    bsub = 1;
    ak = 0;
    bk = 0;
    st.site = &p_emlrtRSI;
    ck = 0;
    do {
      exitg1 = 0;
      i = a->size[1];
      if (ck <= i - 1) {
        st.site = &q_emlrtRSI;
        st.site = &r_emlrtRSI;
        st.site = &s_emlrtRSI;
        c->data[ck] = a->data[ak] - b->data[bk];
        st.site = &t_emlrtRSI;
        if (asub < a->size[1]) {
          st.site = &u_emlrtRSI;
          ak++;
          st.site = &v_emlrtRSI;
          bk++;
          st.site = &w_emlrtRSI;
          bsub++;
          st.site = &x_emlrtRSI;
          asub++;
        } else {
          st.site = &y_emlrtRSI;
          if (bsub < b->size[1]) {
            st.site = &ab_emlrtRSI;
            st.site = &bb_emlrtRSI;
            bk++;
            st.site = &cb_emlrtRSI;
            bsub++;
          } else {
            asub = 1;
            bsub = 1;
          }
        }

        ck++;
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
}

/* End of code generation (bsxfun.c) */
