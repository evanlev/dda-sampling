#include "misc.h"
#include "misc.hpp"
#include "multind.h"
#include "mdarray.h"
#include "debug.h"
#include "dda_utils.h"
#include "sampleHeap.h"
#include "config.h"

#include <algorithm>
#include <vector>
#include <random>
#include <limits>

constexpr static double Inf = std::numeric_limits<double>::max();

void SparseKernel::Print() const
{
    debug_printf(DP_INFO, "sparse w: \n");
    for (long t2 = 0; t2 < m_wKeyValuePairs.size(); t2++)
    {
        for (long t1 = 0; t1 < m_wKeyValuePairs[t2].size(); t1++)
        {
            for (const auto &[location, value] : m_wKeyValuePairs[t2][t1] )
            {
                debug_printf(DP_INFO, "w(");
                for (int i = 0; i < kPhaseEncodeDims; i++)
                {
                    debug_printf(DP_INFO, "%d, ", location[i]);
                }
                debug_printf(DP_INFO, "%d, %d) = %f\n", t2, t1, value);
            }
        }
    }
}

SparseKernel::SparseKernel(const MDArray<kPhaseEncodeDims + 2, double> &kernel, double T)
    : m_dims(kernel.Dims())
{
    const long kSize = md_calc_size(kPhaseEncodeDims, m_dims.data());
    const long nt = m_dims[kPhaseEncodeDims];

    debug_printf(DP_DEBUG3, "sparsifying w, nt: %d\n", nt);
    const std::array<long, kPhaseEncodeDims + 2> strs = md_calc_strides(kernel.Dims(), 1);

    // set sparseKernel.w00tt
    m_w00tt.resize(nt);
    for (long t = 0; t < nt; t++)
    {
        m_w00tt[t] = kernel[t * strs[kPhaseEncodeDims] + t * strs[kPhaseEncodeDims + 1]];
    }

    // Allocate w k-v pairs
    m_wKeyValuePairs.resize(nt);
    for (auto& kv : m_wKeyValuePairs)
    {
        kv.resize(nt);
    }

    // Assign k-v pairs for w
    long supp = 0;
    for (long t2 = 0; t2 < nt; t2++)
    {
        for (long t1 = 0; t1 < nt; t1++)
        {
            for (long k = 0; k < kSize; k++)
            {
                std::array<long, kPhaseEncodeDims> kSub;
                ind2sub(kPhaseEncodeDims, kernel.Dims().data(), kSub.data(), k);
                long kt1t2r_ind = k + t1 * strs[kPhaseEncodeDims] + t2 * strs[kPhaseEncodeDims + 1];
                if (kernel[kt1t2r_ind] > T)
                {
                    // append (k) to w_sp[t1][t2]
                    supp++;
                    m_wKeyValuePairs[t1][t2].emplace_back(kSub, kernel[kt1t2r_ind]);
                }
            }
        }
    }
    debug_printf(DP_INFO, "w support: %d/%d = %f\n", supp, nt * nt * kSize, (float)supp / (float)(nt * nt * kSize));
}

void printPat(const MDArray<kPhaseEncodeDims+1, int> &mask)
{
    static int maskPrintCount = 1;
    debug_printf(DP_INFO, "mask(:,%d) = [", maskPrintCount++);
    for (long i = 0; i < mask.Length(); i++)
    {
        debug_printf(DP_INFO, "%d ", mask[i]);
    }
    debug_printf(DP_INFO, "];\n");
}

/*
 * Best candidate sampling but use a heap to store J
 */
template <size_t PhaseEncodeDims>
void SampleApproxBestCandidate(const SparseKernel &wsp,
                               const vector<long> &maxSamplesPerFrame,
                               long totSamps,
                               MDArray<PhaseEncodeDims + 1, int> &mask,
                               double &cost)
{
    const long nt = wsp.Dims()[PhaseEncodeDims];
    const long kSize = md_calc_size(PhaseEncodeDims, wsp.Dims().data());

    debug_printf(DP_INFO, "Approximate best candidate sampling, nt: %d\n", nt);

    std::array<long, PhaseEncodeDims + 2> strs = md_calc_strides(wsp.Dims(), 1);

    MDArray<PhaseEncodeDims+1, double> deltaJ(mask.Dims());

    // Initialize deltaJ
    for (long t = 0; t < nt; t++)
    {
        double currentDeltaJ = wsp.GetW00tt(t);
        md_fill(PhaseEncodeDims, wsp.Dims().data(), &deltaJ[t * kSize], &currentDeltaJ, sizeof(double));
    }

    // Clear out the mask
    mask.Clear();

    // Current number of samples
    std::vector<long> nSamplesCurrent(nt, 0);

    // Optional constraint on max samples in each pattern
    int numFull = 0;
    std::vector<bool> frameIsFull(nt);
    for (long t = 0; t < nt; t++)
    {
        frameIsFull[t] = maxSamplesPerFrame[t] == 0;
        if (frameIsFull[t])
        {
            numFull++;
        }
    }

    // Heapify delta J, later can use the fact that deltaJ is constant to speed this up
    debug_printf(DP_INFO, "Construct heap...\n");
    SampleHeap heap(deltaJ);

    // Best candidate selection loop
    debug_printf(DP_INFO, "Main loop...\n");
    cost = 0;
    for (long i = 0; i < totSamps; i++)
    {
        // TODO: pop or have an option if we can reacquire
        const Sample &snew = heap.getArr(0);

        std::array<long, PhaseEncodeDims + 1> ktNewSub;
        ind2sub(mask.Dims(), ktNewSub.data(), snew.getKTIndex());
        const long tNew = ktNewSub[PhaseEncodeDims];

        debug_printf(DP_DEBUG4, "Sampling (%d %d)\n", 1 + ktNewSub[0], 1 + ktNewSub[1]);

        mask[snew.getKTIndex()]++;
        nSamplesCurrent[tNew]++;

        // --- Update total cost
        cost += deltaJ[snew.getKTIndex()];

        const SparseKernel::WType& weighting = wsp.GetW();

        // --- Update insert cost (deltaJ), even for frames that are full
        for (long t = 0; t < nt; t++)
        {
            // w( Delta k, t, t' )
            for (const auto& [location, value] : weighting[t][tNew])
            {
                // tmp = k' + (Delta k)_j
                long tmp[PhaseEncodeDims + 1];
                add(PhaseEncodeDims, tmp, &ktNewSub[0], location.data());
                mod(PhaseEncodeDims, tmp, tmp, wsp.Dims().data());
                tmp[PhaseEncodeDims] = t;
                long i1 = sub2ind(strs, tmp);
                heap.increaseKey(heap.getKt2idx(i1), value);
            }

            // w( Delta k, t', t)
            for (const auto& [location, value] : weighting[tNew][t] )
            {
                // tmp = k' - (Delta k)_j
                long tmp[PhaseEncodeDims + 1];
                sub(PhaseEncodeDims, tmp, &ktNewSub[0], location.data());
                mod(PhaseEncodeDims, tmp, tmp, wsp.Dims().data());
                tmp[PhaseEncodeDims] = t;
                long i1 = sub2ind(strs, tmp);
                // deltaJ[k' + (Delta k)_j] += w[t'][t][j]
                heap.increaseKey(heap.getKt2idx(i1), value);
            }
        } // end t loop

        // --- Check if this pattern is full
        if (nSamplesCurrent[tNew] == maxSamplesPerFrame[tNew])
        {
            frameIsFull[tNew] = true;
            numFull++;
            debug_printf(DP_DEBUG3, "%d is full\n", tNew);
            // TODO: do something about this
            //md_fill(D, wsp.Dims(), &deltaJ[tNew*kSize], &Inf, sizeof(double));
        }

        if (numFull == nt)
        {
            debug_printf(DP_DEBUG3, "All frames full, breaking at %d/%d\n", i + 1, totSamps);
            break;
        }

        if (totSamps > 10 && totSamps % (totSamps / 10) == 0)
        {
            debug_printf(DP_DEBUG1, "\n%d%", (int)(100.0 * (float)i / (float)totSamps));
        }
    } // end i loop

    debug_printf(DP_INFO, "Done!\n");
}

// Functor to get the minimum element of an array breaking ties randomly. This does bias toward choosing elements
struct GetArgmin
{
    explicit GetArgmin(const int N)
    {
        for (long i = 0; i < N; i++)
        {
            m_perm.push_back(i);
        }
        std::random_shuffle(m_perm.begin(), m_perm.end());
    }

    double operator()(const MDArray<kPhaseEncodeDims+1, double> &array) const
    {
        assert(array.Length() > 0);
        long argmin = 0;
        double min = array[argmin];
        for (long j = 0; j < array.Length(); j++)
        {
            long i = m_perm[j];
            // long i = j;
            if (array[i] < min || (array[i] == min && (rand() % 2 < 1)))
            {
                argmin = i;
                min = array[i];
            }
        }
        return argmin;
    }
private:
    std::vector<long> m_perm;
};

/*
 * Best candidate sampling
 * compute deltaJ, cost change associated with adding a sample 
 *
 * INPUTS:
 *   w        = weighting cost function
 *   maxSamplesPerFrame = max samples per frame
 *   totSamps = total number of samples
 *
 * OUTPUTS:
 *   cost    = final cost
 *   deltaJ   = final insertion cost
 *   mask    = final sampling pattern
 *   samples = sample list
 */
template <size_t PhaseEncodeDims>
void SampleExactBestCandidate(const MDArray<PhaseEncodeDims + 2, double> &kernel,
                              const std::vector<long> &maxSamplesPerFrame,
                              long totSamps,
                              MDArray<PhaseEncodeDims + 1, int> &mask,
                              double &cost,
                              MDArray<PhaseEncodeDims + 1, double> &deltaJ)
{
    const long nt = kernel.Dims()[PhaseEncodeDims];
    const long kSize = md_calc_size(PhaseEncodeDims, kernel.Dims().data());
    const long cSize = nt * kSize;

    // Initialize deltaJ
    for (long t = 0; t < nt; t++)
    {
        double currentDeltaJ = kernel[t * kernel.Strides()[PhaseEncodeDims] + t * kernel.Strides()[PhaseEncodeDims + 1]];
        md_fill(PhaseEncodeDims, kernel.Dims().data(), &deltaJ[t * kSize], &currentDeltaJ, sizeof(double));
    }

    // Clear out the mask
    mask.Clear();

    // Current number of samples
    std::vector<long> nSamplesCurrent(nt, 0);

    // Optional constraint on max samples in each pattern
    long numFull = 0;
    std::vector<bool> frameIsFull(nt);
    for (long t = 0; t < nt; t++)
    {
        frameIsFull[t] = maxSamplesPerFrame[t] == 0;
        if (frameIsFull[t])
            numFull++;
    }

    GetArgmin argmin(deltaJ.Length());

    // Best candidate selection loop
    cost = 0;
    for (long i = 0; i < totSamps; i++)
    {
        // k_new = argmin deltaJ(k,t)
        long ktNewInd = argmin(deltaJ);
        long ktNewSub[PhaseEncodeDims + 1];
        ind2sub(mask.Dims(), ktNewSub, ktNewInd);
        long tNew = ktNewSub[PhaseEncodeDims];
        myAssert(!frameIsFull[tNew], "argmin should not have selected a full frame!");

        // --- Add sample
        // Increment mask
        mask[ktNewInd]++;
        nSamplesCurrent[tNew]++;

        // --- Update total cost
        cost += deltaJ[ktNewInd];

        // --- Update insert cost (deltaJ), even for frames that are full
        for (long sj_ind = 0; sj_ind < cSize; sj_ind++)
        {
            // diff = (ktNewSub - sj, tNew, tj)
            long sj[PhaseEncodeDims + 1];
            ind2sub(mask.Dims(), sj, sj_ind);
            const long tj = sj[PhaseEncodeDims];

            long diff[PhaseEncodeDims + 2];
            sub(PhaseEncodeDims, diff, ktNewSub, sj);
            mod(PhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[PhaseEncodeDims + 0] = tNew;
            diff[PhaseEncodeDims + 1] = tj;

            // deltaJ[sj_ind] += w[snew - sj, tNew, tj]
            deltaJ[sj_ind] += kernel[sub2ind(kernel.Strides(), diff)];

            // diff = (sj - ktNewSub, tj, tNew)
            sub(PhaseEncodeDims, diff, sj, ktNewSub);
            mod(PhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[PhaseEncodeDims + 0] = tj;
            diff[PhaseEncodeDims + 1] = tNew;

            // deltaJ[sj_ind] += w[sj - snew, tj, tNew]
            deltaJ[sj_ind] += kernel[sub2ind(kernel.Strides(), diff)];
        } // end sj_ind loop

        // --- Check if this pattern is full
        if (nSamplesCurrent[tNew] == maxSamplesPerFrame[tNew])
        {
            frameIsFull[tNew] = true;
            numFull++;
            debug_printf(DP_DEBUG3, "%d is full\n", tNew);
            md_fill(PhaseEncodeDims, kernel.Dims().data(), &deltaJ[tNew * kSize], &Inf, sizeof(double));
        }

        // --- Break if all patterns are full
        if (numFull == nt)
        {
            debug_printf(DP_DEBUG3, "All frames full, breaking at %d/%d\n", i, totSamps);
            break;
        }

        // --- Progress
        if (totSamps % (totSamps / 10) == 0)
        {
            debug_printf(DP_DEBUG1, "\r%d%", (int)(100.0 * (float)i / (float)totSamps));
        }
    } // end i loop
}

template <size_t PhaseEncodeDims>
MDArray<PhaseEncodeDims+1, double> ComputeDeltaJ(const MDArray<PhaseEncodeDims + 2, double> &kernel,
                                                 const MDArray<PhaseEncodeDims+1, int> &mask)
{

    MDArray<PhaseEncodeDims+1, double> deltaJ(mask.Dims());

    debug_printf(DP_DEBUG1, "Computing first term in DeltaJ...\n");
    const long nframes = kernel.Dims()[PhaseEncodeDims];
    const long kSize = md_calc_size(PhaseEncodeDims, kernel.Dims().data());
    const long cSize = nframes * kSize;
    for (long t = 0; t < nframes; t++)
    {
        double currentDeltaJ = kernel[t * kernel.Strides()[PhaseEncodeDims] + t * kernel.Strides()[PhaseEncodeDims + 1]];
        md_fill(PhaseEncodeDims, kernel.Dims().data(), &deltaJ[t * kSize], &currentDeltaJ, sizeof(double));
    }

    std::array<long, PhaseEncodeDims + 1> ktNewSub;
    std::array<long, PhaseEncodeDims + 1> ktSampleSub;
    for (long ktNewInd = 0; ktNewInd < cSize; ktNewInd++)
    {
        for (long ktSampleInd = 0; ktSampleInd < cSize; ktSampleInd++)
        {
            if (mask[ktSampleInd] == 0)
            {
                continue;
            }

            ind2sub(mask.Dims(), ktNewSub.data(), ktNewInd);
            ind2sub(mask.Dims(), ktSampleSub.data(), ktSampleInd);

            // diff = (k_new - si, tNew, t)
            long diff[PhaseEncodeDims + 2];
            sub(PhaseEncodeDims, diff, ktNewSub.data(), ktSampleSub.data());
            mod(PhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[PhaseEncodeDims + 0] = ktNewSub[PhaseEncodeDims];
            diff[PhaseEncodeDims + 1] = ktSampleSub[PhaseEncodeDims];

            // deltaJ[k_new] += w[diff]
            deltaJ[ktNewInd] += static_cast<double>(mask[ktSampleInd]) *
                                 kernel[sub2ind(kernel.Strides(), diff)];

            // diff = (si - ktNewSub, t, tNew);
            sub(PhaseEncodeDims, diff, ktSampleSub.data(), ktNewSub.data());
            mod(PhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[PhaseEncodeDims + 0] = ktSampleSub[PhaseEncodeDims];
            diff[PhaseEncodeDims + 1] = ktNewSub[PhaseEncodeDims];

            // deltaJ[k_new] += w[diff]
            deltaJ[ktNewInd] += static_cast<double>(mask[ktSampleInd]) *
                                 kernel[sub2ind(kernel.Strides(), diff)];

        } // end i loop
        debug_printf(DP_DEBUG1, "\r%d%", 100.0 * (float)ktNewInd / (float)cSize);
    } // end k new loop
    debug_printf(DP_DEBUG1, "Done computing first term in DeltaJ.\n");

    return deltaJ;
}

/* 
 * Note this is slower than the method using the Fourier transform or the method that 
 * exploits the sparsity of the mask.
 */
template <size_t PhaseEncodeDims>
MDArray<PhaseEncodeDims + 2, double> computeDiffDist(const MDArray<PhaseEncodeDims + 1, int> &mask)
{
    const std::array<long, PhaseEncodeDims + 1>& maskDims = mask.Dims();
    const std::array<long, PhaseEncodeDims + 2> pDims = {maskDims[0], maskDims[1], maskDims[2], maskDims[2]};
    MDArray<PhaseEncodeDims + 2, double> p(pDims);
    p.Clear();

    long diff[PhaseEncodeDims + 2];
    std::array<long, PhaseEncodeDims + 1> si, sj;
    for (long i = 0; i < mask.Length(); i++)
    {
        for (long j = 0; j < mask.Length(); j++)
        {
            if (mask[i] == 0 || mask[j] == 0)
            {
                continue;
            }

            ind2sub(mask.Dims(), si.data(), i);
            ind2sub(mask.Dims(), sj.data(), j);

            sub(PhaseEncodeDims, diff, si.data(), sj.data());
            mod(PhaseEncodeDims, diff, diff, p.Dims().data());
            diff[PhaseEncodeDims + 0] = si[PhaseEncodeDims];
            diff[PhaseEncodeDims + 1] = sj[PhaseEncodeDims];
#ifdef DEBUG
            assert_in_bounds(PhaseEncodeDims, diff, p.Dims().data(), "diff out of bounds of mask");
#endif
            p[sub2ind(p.Strides(), diff)]++;
        }
    }
    return p;
}

// Explicit instantation.
template
MDArray<kPhaseEncodeDims + 1, double> ComputeDeltaJ<kPhaseEncodeDims>(const MDArray<kPhaseEncodeDims+2, double> &kernel,
                                                                     const MDArray<kPhaseEncodeDims+1, int> &mask);

template
void SampleApproxBestCandidate<kPhaseEncodeDims>(const SparseKernel &wsp,
                                                          const std::vector<long> &maxSamplesPerFrame,
                                                          long totSamps,
                                                          MDArray<kPhaseEncodeDims + 1, int> &mask,
                                                          double &cost);

template
MDArray<kPhaseEncodeDims+2, double> computeDiffDist<kPhaseEncodeDims>(const MDArray<kPhaseEncodeDims + 1, int> &mask);

template
void SampleExactBestCandidate<kPhaseEncodeDims>(const MDArray<kPhaseEncodeDims + 2, double> &kernel,
                                                const std::vector<long> &maxSamplesPerFrame,
                                                long totSamps,
                                                MDArray<kPhaseEncodeDims + 1, int> &mask,
                                                double &cost,
                                                MDArray<kPhaseEncodeDims+1, double> &deltaJ);
