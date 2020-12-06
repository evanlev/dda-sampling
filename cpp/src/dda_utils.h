#ifndef DDA_UTILS_H
#define DDA_UTILS_H 1u

#include "mdarray.h"
#include "misc.h"
#include "config.h"

#define MAX_T 100

// TODO clean up
struct SparseKernel
{
    long *kw[MAX_T][MAX_T];
    long D;
    long dims[kPhaseEncodeDims + 2u];
    double w00tt[MAX_T];

    // list of nonzeros of w
    double *w[MAX_T][MAX_T];
    long len_w[MAX_T][MAX_T];

    void Print() const;
};

// TODO test this function
MDArray<4, double> computeDiffDist(const vector<vector<long> > &samples, const long dims[]);

MDArray<3, double> computeDeltaJ(const MDArray<4, double> &w, const MDArray<3, int> &mask);

SparseKernel sparsifyWToK(const MDArray<4, double> &wmd, const long k);

SparseKernel sparsifyW(const MDArray<4, double> &wmd, const double T);

void exactBestCandidate(const MDArray<4, double> &kernel, const vector<long> &maxSamps, const long totSamps, 
    MDArray<3, int> &mask, double &cost, MDArray<3, double> &deltaJ);

void approxBestCandidate(const SparseKernel& sparse_kernel,
                         const vector<long> &maxSamps, 
                         const long totSamps,
                         MDArray<3, int> &mask,
                         double &cost);
#endif // #ifndef DDA_UTILS_H
