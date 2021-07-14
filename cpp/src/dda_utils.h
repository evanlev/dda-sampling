#ifndef DDA_UTILS_H
#define DDA_UTILS_H 1u

#include "mdarray.h"
#include "misc.h"
#include "config.h"

#include <array>
#include <vector>

class SparseKernel
{
public:
    SparseKernel(const MDArray<kPhaseEncodeDims + 2, double> &kernel, double T);

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

template <size_t PhaseEncodeDims>
MDArray<PhaseEncodeDims + 2, double> computeDiffDist(const std::vector<std::vector<long>> &samples, const long dims[]);

template <size_t PhaseEncodeDims>
MDArray<PhaseEncodeDims + 1, double> ComputeDeltaJ(const MDArray<PhaseEncodeDims + 2, double> &w,
                                                   const MDArray<PhaseEncodeDims + 1, int> &mask);

template <size_t PhaseEncodeDims>
void SampleExactBestCandidate(const MDArray<PhaseEncodeDims + 2, double> &kernel,
                              const std::vector<long> &maxSamplesPerFrame,
                              long totSamps,
                              MDArray<PhaseEncodeDims + 1, int> &mask,
                              double &cost,
                              MDArray<PhaseEncodeDims + 1, double> &deltaJ);

template <size_t PhaseEncodeDims>
void SampleApproxBestCandidate(const SparseKernel& sparse_kernel,
                               const std::vector<long> &maxSamps, 
                               long totSamps,
                               MDArray<PhaseEncodeDims + 1, int> &mask,
                               double &cost);
#endif // #ifndef DDA_UTILS_H
