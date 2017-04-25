

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

void compute_dd_sparse_pat(const int D, const long dims[], double *p, const long *samples[], const long N[]);

void compute_dd_fft(const int D, const long dimsp[], int *p, const int *pat);

void computeDeltaJ(const int D, const long dims[], double *deltaJ, const double *w, long *samples[], const long Nsamps[]);
void computeDeltaJ2(const int D, const long dims[], double *deltaJ, const double *w, long *samples[], const long Nsamps[], const int w_type);

void circRev2d(const long dims[], complex float *dst, const complex float *src);
void computeDeltaJFftMethod(const int kDims, const long w_dims[], complex float *deltaJ, const complex float *w, const complex float *pat);

void check_sns_dims(const long sns_dims[]);
void get_w_dims(long dims[], const long sns_dims[]);

// Convert sensitivity maps to w, complex output, but real part = 0
void buildW2(complex float *w, const long sns_dims[], 
            const complex float *sns_maps);


// Convert sensitivity maps to w, real output
void buildW(double *w, const long sns_dims[], 
             const complex float *sns_maps);

void exactBestCandidate(const int D, const long dims[], double *cost, double *deltaJ, 
                    int *mask, const double *w, 
                    long *samples[], const long maxSamps[], 
                    const long totSamps);

void approxBestCandidate(const int D, double *cost, double *deltaJ, 
                                 int *mask, const SparseW* wsp,
                                 const long maxSamps[], 
                                 const long totSamps);

void sparsify_w_to_k(const long D, const long dims[], SparseW *wsp, double *wmd, const long k);

void sparsify_w(const long D, const long dims[], SparseW *wsp, double *wmd, const double T);

void print_wsp(const SparseW *wsp);


