#ifndef __DDA_UTILS_H
#define __DDA_UTILS_H

#include "misc/cppwrap.h"

#define KSP_DIMS 2 // 2 k-space dimensions
#define MAX_T 100
#define MAPS2_DIM (TE_DIM+1)
#define MAPS2_FLAG (1u << MAPS2_DIM)
#include "misc/mri.h" // DIMS

typedef struct sparsew_s{
    long *kw[MAX_T][MAX_T][2];
    long D;
    long dims[DIMS];
    double w00tt[MAX_T];
    // list of nonzeros of w
    double *w[MAX_T][MAX_T][2];
    long len_w[MAX_T][MAX_T][2];
} SparseW;

extern void compute_dd_sparse_pat(const int D, const long dims[], double *p, const long *samples[], const long N[]);

extern void compute_dd_fft(const int D, const long dimsp[], int *p, const int *pat);

extern void computeDeltaJ(const int D, const long dims[], double *deltaJ, const double *w, long *samples[], const long Nsamps[]);

extern void circRev2d(const long dims[], _Complex float *dst, const _Complex float *src);
extern void computeDeltaJFftMethod(const int kDims, const long w_dims[], _Complex float *deltaJ, const _Complex float *w, const _Complex float *pat);

extern void check_sns_dims(const long sns_dims[]);
extern void get_w_dims(long dims[], const long sns_dims[]);

// Convert sensitivity maps to w, _Complex output, but real part = 0
extern void buildW2(_Complex float *w, const long sns_dims[], 
            const _Complex float *sns_maps);


// Convert sensitivity maps to w, real output
extern void buildW(double *w, const long sns_dims[], 
             const _Complex float *sns_maps);

extern void exactBestCandidate(const int D, const long dims[], double *cost, double *deltaJ, 
                    int *mask, const double *w, 
                    long *samples[], const long maxSamps[], 
                    const long totSamps);

extern void approxBestCandidate(const int D, double *cost, double *deltaJ, 
                                 int *mask, const SparseW* wsp,
                                 const long maxSamps[], 
                                 const long totSamps);

extern void sparsify_w_to_k(const long D, const long dims[], SparseW *wsp, double *wmd, const long k);

extern void sparsify_w(const long D, const long dims[], SparseW *wsp, double *wmd, const double T);

extern void print_wsp(const SparseW *wsp);

#include "misc/cppwrap.h"

#endif	// __DDA_UTILS_H
