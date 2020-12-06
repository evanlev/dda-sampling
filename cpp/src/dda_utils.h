#ifndef DDA_UTILS_H
#define DDA_UTILS_H 1u

#include "mdarray.h"
#include "misc.h"
#include "config.h"

#include <array>

static constexpr size_t kMaxPatterns = 100;

// TODO clean up. This is used more like a class
struct SparseKernel
{
    SparseKernel(const MDArray<kPhaseEncodeDims + 2, double> &kernel, double T);
    SparseKernel(const MDArray<kPhaseEncodeDims + 2, double> &kernel, long K);

    using WType = std::vector<std::vector<std::vector<std::pair<std::array<long, kPhaseEncodeDims>, double>>>>;

    const std::array<long, kPhaseEncodeDims + 2>& Dims() const
    {
        return m_dims;
    }

    double GetW00tt(unsigned int t) const
    {
        return m_w00tt[t];
    }

    const WType& GetW() const
    {
        return m_wKeyValuePairs;
    }

    void Print() const;

private:
    // Stores key-value pairs for w. m_w[t1][t2][j] contains the location and value in the first and second
    // elements.
    WType m_wKeyValuePairs;

    // holds w(0, 0, t, t)
    std::vector<double> m_w00tt;

    std::array<long, kPhaseEncodeDims + 2> m_dims;
};

MDArray<kPhaseEncodeDims + 2, double> computeDiffDist(const vector<vector<long>> &samples, const long dims[]);

MDArray<kPhaseEncodeDims + 1, double> computeDeltaJ(const MDArray<4, double> &w, const MDArray<kPhaseEncodeDims + 1, int> &mask);

void exactBestCandidate(const MDArray<4, double> &kernel, const vector<long> &maxSamps, const long totSamps, 
    MDArray<3, int> &mask, double &cost, MDArray<3, double> &deltaJ);

void approxBestCandidate(const SparseKernel& sparse_kernel,
                         const vector<long> &maxSamps, 
                         const long totSamps,
                         MDArray<3, int> &mask,
                         double &cost);
#endif // #ifndef DDA_UTILS_H
