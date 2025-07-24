#include "mex.h"
#include <cblas.h>

/*
 * matvec_blas_mex.c
 *
 * Compute y = A * v efficiently by skipping zero rows/cols
 * based on a logical mask, using BLAS for the inner mat-vec on the
 * compacted active block.
 *
 * Syntax:
 *   y = matvec_blas_mex(A, v, mask)
 *
 * Inputs (assumed valid):
 *   A    - n-by-n double matrix (column-major as MATLAB stores)
 *   v    - double vector of length n
 *   mask - logical vector of length n; true = include, false = skip
 *
 * Output:
 *   y    - double vector of length n, with y(i)=0 if mask(i)==false
 *
 * Strategy:
 * 1) Count active indices (m), allocate contiguous B[m*m] on heap.
 * 2) Extract B from A into column-major buffer.
 * 3) Extract v_active into contiguous buffer.
 * 4) Call cblas_dgemv to compute y_active = B * v_active.
 * 5) Scatter y_active back into y.
 *
 * No per-element loops in the inner product.
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A    = mxGetDoubles(prhs[0]);  /* n x n */
    double *v    = mxGetDoubles(prhs[1]);  /* n x 1 */
    mxLogical *mask = mxGetLogicals(prhs[2]); /* n x 1 */

    mwSize n = mxGetM(prhs[0]);

    /* Count active entries */
    mwSize i, j, m = 0;
    for (i = 0; i < n; ++i) if (mask[i]) ++m;

    /* Prepare output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *y = mxGetDoubles(plhs[0]);
    for (i = 0; i < n; ++i) y[i] = 0.0;

    if (m == 0) return;

    /* Allocate buffers */
    double *B      = (double*)mxMalloc(m * m * sizeof(double));
    double *v_act  = (double*)mxMalloc(m * sizeof(double));
    double *y_act  = (double*)mxMalloc(m * sizeof(double));
    mwSize *idx    = (mwSize*)mxMalloc(m * sizeof(mwSize));

    /* Build active idx list */
    mwSize p = 0;
    for (i = 0; i < n; ++i) if (mask[i]) idx[p++] = i;

    /* Extract B and v_act */
    for (j = 0; j < m; ++j) {
        mwSize col = idx[j];
        v_act[j] = v[col];
        for (i = 0; i < m; ++i) {
            mwSize row = idx[i];
            B[i + j*m] = A[row + col*n];
        }
    }

    /* Compute y_act = B * v_act */
    cblas_dgemv(CblasColMajor, CblasNoTrans,
                (int)m, (int)m,
                1.0, B, (int)m,
                v_act, 1,
                0.0, y_act, 1);

    /* Scatter result */
    for (j = 0; j < m; ++j) {
        y[idx[j]] = y_act[j];
    }

    /* Cleanup */
    mxFree(B);
    mxFree(v_act);
    mxFree(y_act);
    mxFree(idx);
}

