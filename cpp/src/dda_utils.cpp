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
    int nt = m_dims[kPhaseEncodeDims];
    for (int t2 = 0; t2 < nt; t2++)
    {
        for (int t1 = 0; t1 < nt; t1++)
        {
            for (int j = 0; j < m_wKeyValuePairs[t1][t2].size(); j++)
            {
                debug_printf(DP_INFO, "w(");
                for (int i = 0; i < kPhaseEncodeDims; i++)
                {
                    debug_printf(DP_INFO, "%d, ", m_wKeyValuePairs[t1][t2][j * kPhaseEncodeDims + i]);
                }
                debug_printf(DP_INFO, "%d, %d) = %f\n", t1, t2, m_wKeyValuePairs[t1][t2][j]);
            }
        }
    }
}

SparseKernel::SparseKernel(const MDArray<kPhaseEncodeDims + 2, double> &kernel, double T)
    : m_dims(kernel.Dims())
{
    long kSize = md_calc_size(kPhaseEncodeDims, m_dims.data());
    long nt = m_dims[kPhaseEncodeDims];

    debug_printf(DP_DEBUG3, "sparsifying w, nt: %d\n", nt);
    std::array<long, kPhaseEncodeDims + 2> strs;
    md_calc_strides(kPhaseEncodeDims + 2, strs.data(), kernel.Dims().data(), 1);

    // set sparseKernel.w00tt
    m_w00tt.resize(nt);
    for (long t = 0; t < nt; t++)
    {
        m_w00tt[t] = kernel[t * strs[kPhaseEncodeDims] + t * strs[kPhaseEncodeDims + 1]];
    }

    // Allocate w k-v pairs
    m_wKeyValuePairs.resize(nt);
    for (long t2 = 0; t2 < nt; t2++)
    {
        m_wKeyValuePairs[t2].resize(nt);
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
                if (kernel.m_data[kt1t2r_ind] > T)
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
void approxBestCandidate(const SparseKernel &wsp,
                         const vector<long> &maxSamplesPerFrame,
                         const long totSamps, MDArray<kPhaseEncodeDims+1, int> &mask,
                         double &cost)
{
    //debug_level = DP_ALL;

    const long nt = wsp.Dims()[kPhaseEncodeDims];
    const long kSize = md_calc_size(kPhaseEncodeDims, wsp.Dims().data());

    debug_printf(DP_INFO, "Approximate best candidate sampling, nt: %d\n", nt);

    long strs[kPhaseEncodeDims + 2];
    md_calc_strides(kPhaseEncodeDims + 2, strs, wsp.Dims().data(), 1);

    // Total cost
    cost = 0;
    MDArray<kPhaseEncodeDims+1, double> deltaJ(mask.Dims());

    // Initialize deltaJ
    for (long t = 0; t < nt; t++)
    {
        double currentDeltaJ = wsp.GetW00tt(t);
        md_fill(kPhaseEncodeDims, wsp.Dims().data(), &deltaJ[t * kSize], &currentDeltaJ, sizeof(double));
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
    for (long i = 0; i < totSamps; i++)
    {
        // TODO: pop or have an option if we can reacquire
        const Sample &snew = heap.getArr(0);

        std::array<long, kPhaseEncodeDims + 1> ktNewSub;
        ind2sub(kPhaseEncodeDims + 1, wsp.Dims().data(), ktNewSub.data(), snew.getKTIndex());
        const long tNew = ktNewSub[kPhaseEncodeDims];

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
            for (long j = 0; j < weighting[t][tNew].size(); j++)
            {
                // tmp = k' + (Delta k)_j
                const auto& [location, value] = weighting[tNew][t][j];
                long tmp[kPhaseEncodeDims + 1];
                add(kPhaseEncodeDims, tmp, &ktNewSub[0], location.data());
                mod(kPhaseEncodeDims, tmp, tmp, wsp.Dims().data());
                tmp[kPhaseEncodeDims] = t;
                long i1 = sub2ind(kPhaseEncodeDims + 1, strs, tmp);
                heap.increaseKey(heap.getKt2idx(i1), value);
            }

            // w( Delta k, t', t)
            for (long j = 0; j < weighting[tNew][t].size(); j++)
            {
                // tmp = k' - (Delta k)_j
                const auto& [location, value] = weighting[tNew][t][j];
                long tmp[kPhaseEncodeDims + 1];
                sub(kPhaseEncodeDims, tmp, &ktNewSub[0], location.data());
                mod(kPhaseEncodeDims, tmp, tmp, wsp.Dims().data());
                tmp[kPhaseEncodeDims] = t;
                long i1 = sub2ind(kPhaseEncodeDims + 1, strs, tmp);
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

        // --- Break if all patterns are full
        if (numFull == nt)
        {
            debug_printf(DP_DEBUG3, "All frames full, breaking at %d/%d\n", i + 1, totSamps);
            break;
        }

        // --- Progress
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
    std::vector<long> m_perm;

    GetArgmin(const int N)
    {
        for (long i = 0; i < N; i++)
            m_perm.push_back(i);
        std::random_shuffle(m_perm.begin(), m_perm.end());
    }

    double operator()(const MDArray<kPhaseEncodeDims+1, double> &array) const
    {
        assert(array.Length() > 0);
        long argmin = 0;
        double min = array[argmin];
        for (long j = 0; j < array.Length(); j++)
        {
            //long i = m_perm[j];
            long i = j;
            if (array[i] < min || (array[i] == min && (rand() % 2 < 1)))
            {
                argmin = i;
                min = array[i];
            }
        }
        return argmin;
    }
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
void exactBestCandidate(const MDArray<kPhaseEncodeDims + 2, double> &kernel,
                        const std::vector<long> &maxSamplesPerFrame,
                        const long totSamps, MDArray<kPhaseEncodeDims + 1, int> &mask, double &cost,
                        MDArray<kPhaseEncodeDims+1, double> &deltaJ)
{
    const long nt = kernel.Dims()[kPhaseEncodeDims];
    const long kSize = md_calc_size(kPhaseEncodeDims, kernel.Dims().data());
    const long cSize = nt * kSize;

    // Total cost
    cost = 0;

    // Initialize deltaJ
    for (long t = 0; t < nt; t++)
    {
        double currentDeltaJ = kernel.m_data[t * kernel.m_strs[kPhaseEncodeDims] + t * kernel.m_strs[kPhaseEncodeDims + 1]];
        md_fill(kPhaseEncodeDims, kernel.Dims().data(), &deltaJ.m_data[t * kSize], &currentDeltaJ, sizeof(double));
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
    for (long i = 0; i < totSamps; i++)
    {
        // k_new = argmin deltaJ(k,t)
        long ktNew_ind = argmin(deltaJ);
        long ktNewSub[kPhaseEncodeDims + 1];
        ind2sub(kPhaseEncodeDims + 1, mask.Dims().data(), ktNewSub, ktNew_ind);
        long tNew = ktNewSub[kPhaseEncodeDims];
        myAssert(!frameIsFull[tNew], "argmin should not have selected a full frame!");

        // --- Add sample
        // Increment mask
        mask[ktNew_ind]++;
        nSamplesCurrent[tNew]++;

        // --- Update total cost
        cost += deltaJ[ktNew_ind];

        // --- Update insert cost (deltaJ), even for frames that are full
        for (long sj_ind = 0; sj_ind < cSize; sj_ind++)
        {

            // diff = (ktNewSub - sj, tNew, tj)
            long sj[kPhaseEncodeDims + 1];
            ind2sub(kPhaseEncodeDims + 1, mask.Dims().data(), sj, sj_ind);
            const long tj = sj[kPhaseEncodeDims];

            long diff[kPhaseEncodeDims + 2];
            sub(kPhaseEncodeDims, diff, ktNewSub, sj);
            mod(kPhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[kPhaseEncodeDims + 0] = tNew;
            diff[kPhaseEncodeDims + 1] = tj;

            // deltaJ[sj_ind] += w[snew - sj, tNew, tj]
            deltaJ.m_data[sj_ind] += kernel.m_data[sub2ind(kPhaseEncodeDims + 2, kernel.m_strs, diff)];

            // diff = (sj - ktNewSub, tj, tNew)
            sub(kPhaseEncodeDims, diff, sj, ktNewSub);
            mod(kPhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[kPhaseEncodeDims + 0] = tj;
            diff[kPhaseEncodeDims + 1] = tNew;

            // deltaJ[sj_ind] += w[sj - snew, tj, tNew]
            deltaJ.m_data[sj_ind] += kernel.m_data[sub2ind(kPhaseEncodeDims + 2, kernel.m_strs, diff)];
        } // end sj_ind loop

        // --- Check if this pattern is full
        if (nSamplesCurrent[tNew] == maxSamplesPerFrame[tNew])
        {
            frameIsFull[tNew] = true;
            numFull++;
            debug_printf(DP_DEBUG3, "%d is full\n", tNew);
            md_fill(kPhaseEncodeDims, kernel.Dims().data(), &deltaJ[tNew * kSize], &Inf, sizeof(double));
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

MDArray<kPhaseEncodeDims+1, double> computeDeltaJ(const MDArray<4, double> &kernel, const MDArray<kPhaseEncodeDims+1, int> &mask)
{

    MDArray<kPhaseEncodeDims+1, double> deltaJ(mask.Dims());

    debug_printf(DP_DEBUG1, "Computing first term in DeltaJ...\n");
    const long nframes = kernel.Dims()[kPhaseEncodeDims];
    const long kSize = md_calc_size(kPhaseEncodeDims, kernel.Dims().data());
    const long cSize = nframes * kSize;

    for (long t = 0; t < nframes; t++)
    {
        double currentDeltaJ = kernel.m_data[t * kernel.m_strs[kPhaseEncodeDims] + t * kernel.m_strs[kPhaseEncodeDims + 1]];
        md_fill(kPhaseEncodeDims, kernel.Dims().data(), &deltaJ[t * kSize], &currentDeltaJ, sizeof(double));
    }

    for (long ktNew_ind = 0; ktNew_ind < cSize; ktNew_ind++)
    {
        for (long kt_sample_ind = 0; kt_sample_ind < cSize; kt_sample_ind++)
        {
            if (mask[kt_sample_ind] == 0)
            {
                continue;
            }

            long ktNewSub[kPhaseEncodeDims + 1];
            ind2sub(kPhaseEncodeDims + 1, mask.Dims().data(), ktNewSub, ktNew_ind);

            long kt_sample_sub[kPhaseEncodeDims + 1];
            ind2sub(kPhaseEncodeDims + 1, mask.Dims().data(), kt_sample_sub, kt_sample_ind);

            // diff = (k_new - si, tNew, t)
            long diff[kPhaseEncodeDims + 2];
            sub(kPhaseEncodeDims, diff, ktNewSub, kt_sample_sub);
            mod(kPhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[kPhaseEncodeDims + 0] = ktNewSub[kPhaseEncodeDims];
            diff[kPhaseEncodeDims + 1] = kt_sample_sub[kPhaseEncodeDims];

            // deltaJ[k_new] += w[diff]
            deltaJ[ktNew_ind] += static_cast<double>(mask[kt_sample_ind]) *
                                  kernel.m_data[sub2ind(kPhaseEncodeDims + 2, kernel.m_strs, diff)];

            // diff = (si - ktNewSub, t, tNew);
            sub(kPhaseEncodeDims, diff, kt_sample_sub, ktNewSub);
            mod(kPhaseEncodeDims, diff, diff, mask.Dims().data());
            diff[kPhaseEncodeDims + 0] = kt_sample_sub[kPhaseEncodeDims];
            diff[kPhaseEncodeDims + 1] = ktNewSub[kPhaseEncodeDims];

            // deltaJ[k_new] += w[diff]
            deltaJ[ktNew_ind] += static_cast<double>(mask[kt_sample_ind]) *
                                  kernel.m_data[sub2ind(kPhaseEncodeDims + 2, kernel.m_strs, diff)];

        } // end i loop
        debug_printf(DP_DEBUG1, "\r%d%", 100.0 * (float)ktNew_ind / (float)cSize);
    } // end k new loop
    debug_printf(DP_DEBUG1, "Done computing first term in DeltaJ.\n");

    return deltaJ;
}

/* 
 * Note this is slower than the method using the Fourier transform or the method that 
 * exploits the sparsity of the mask.
 */
MDArray<4, double> computeDiffDist(const MDArray<kPhaseEncodeDims + 1, int> &mask)
{
    const std::array<long, kPhaseEncodeDims + 1>& maskDims = mask.Dims();
    const std::array<long, kPhaseEncodeDims + 2> p_dims = {maskDims[0], maskDims[1], maskDims[2], maskDims[2]};
    MDArray<kPhaseEncodeDims + 2, double> p(p_dims);
    p.Clear();

    for (int i = 0; i < mask.Length(); i++)
    {

        for (int j = 0; j < mask.Length(); j++)
        {
            if (mask[i] == 0 || mask[j] == 0)
            {
                continue;
            }

            long diff[kPhaseEncodeDims + 2];

            long si[kPhaseEncodeDims + 1u];
            long sj[kPhaseEncodeDims + 1u];

            ind2sub(kPhaseEncodeDims + 1u, mask.Dims().data(), si, i);
            ind2sub(kPhaseEncodeDims + 1u, mask.Dims().data(), sj, j);

            sub(kPhaseEncodeDims, diff, si, sj);
            mod(kPhaseEncodeDims, diff, diff, p.Dims().data());
            diff[kPhaseEncodeDims + 0] = si[kPhaseEncodeDims];
            diff[kPhaseEncodeDims + 1] = sj[kPhaseEncodeDims];
#ifdef DEBUG
            assert_in_bounds(kPhaseEncodeDims, diff, p.Dims().data(), "diff out of bounds of mask");
#endif
            p.m_data[sub2ind(kPhaseEncodeDims + 2, p.m_strs, diff)]++;
        }
    }
    return p;
}
