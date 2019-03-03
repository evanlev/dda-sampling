
#include "mdarray.h"
#include "misc.h"

#define MAX_DIMS 16
#define MAX_T 100

typedef struct sparsew_s{
    long *kw[MAX_T][MAX_T][2];
    long D;
    long dims[MAX_DIMS];
    double w00tt[MAX_T];
    // list of nonzeros of w
    double *w[MAX_T][MAX_T][2];
    long len_w[MAX_T][MAX_T][2];
} SparseW;

void computeDiffDist(const int D, MDArray<double> *p, const long *samples[], const long N[]);

void computeDeltaJ(const int D, MDArray<double> *deltaJ, const MDArray<double> *w, long *samples[], const long Nsamps[]);
void computeDeltaJ2(const int D, double *deltaJ, const MDArray<double> *w, long *samples[], const long Nsamps[], const int w_type);


void sparsifyWToK(const long D, SparseW *wsp, MDArray<double> *wmd, const long k);
void sparsifyW(const long D, SparseW *wsp, MDArray<double> *wmd, const double T);
void printWsp(const SparseW *wsp);

void exactBestCandidate(const int D, double *cost, double *deltaJ, int *mask, const MDArray<double> *w, long *samples[], const long maxSamps[], const long totSamps);

void approxBestCandidate(const int D, 
                         double *cost, 
                         double *deltaJ, 
                         int *mask, 
                         const SparseW* wsp,
                         const long maxSamps[], 
                         const long totSamps);

