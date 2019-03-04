#ifndef DDA_UTILS_H
#define DDA_UTILS_H 1u

#include "mdarray.h"
#include "misc.h"
#include "config.h"

#define MAX_DIMS (kPhaseEncodeDims + 2u)
#define MAX_T 100

struct SparseKernel {
    long *kw[MAX_T][MAX_T];
    long D;
    long dims[MAX_DIMS];
    double w00tt[MAX_T];

    // list of nonzeros of w
    double *w[MAX_T][MAX_T];
    long len_w[MAX_T][MAX_T];

    void Print() const;
};

void computeDiffDist(const int D, const long *samples[], const long N[], MDArray<4, double> &p);

MDArray<3, double> computeDeltaJ(const MDArray<4, double> &w, const vector<vector<long> > &samples);

SparseKernel sparsifyWToK(const MDArray<4, double> &wmd, const long k);
SparseKernel sparsifyW(const MDArray<4, double> &wmd, const double T);

void exactBestCandidate(const int D, double &cost, MDArray<3, double> &deltaJ, 
    int *mask, const MDArray<4, double> &w, long *samples[], const long maxSamps[], const long totSamps);

void approxBestCandidate(const int D, 
                         double &cost, 
                         MDArray<3, double> &deltaJ, 
                         int *mask, 
                         const SparseKernel& wsp,
                         const long maxSamps[], 
                         const long totSamps);
#endif // #ifndef DDA_UTILS_H
