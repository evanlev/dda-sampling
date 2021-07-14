#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>

#include "config.h"
#include "mdarray.h"
#include "multind.h"
#include "misc.h"
#include "dda_utils.h"
#include "debug.h"

// Functor for the window function.
struct WindowFunction
{
    const double beta;
    const double T;

    WindowFunction(double beta, double T) : beta(beta), T(T) {}

    double operator()(const double kr) const
    {
        if (abs(kr) <= (1 - beta) / (2 * T))
        {
            return 1;
        }
        else if (abs(kr) <= (1 + beta) / (2 * T) && abs(kr) >= (1 - beta) / (2 * T))
        {
            //return 0.5 * (1 + cos(3.14159 * T / beta * (abs(kr) - (1 - beta) / (2 * T))));
            return 0.5 * (1 + cos(M_PI * T / beta * (abs(kr) - (1 - beta) / (2 * T))));
        }
        else
        {
            return 0;
        }
    }
};

int main(int argc, char *argv[])
{

    /*
     * TODO: this needs to be tested with the following cases:
     * - More than  1 frame
     * - Max samples constraint is active
     * - 
    */

    // Dimensions
    const std::array<long, kPhaseEncodeDims + 1> patDims = {50, 60, 1};
    const long reduction_factor = 12;
    const long nt = patDims[kPhaseEncodeDims];

    const std::array<long, kPhaseEncodeDims + 2> kernelDims = {patDims[0], patDims[1], patDims[2], patDims[2]};
    const long N = md_calc_size(patDims);
    const long totSamps = N / reduction_factor;
    const std::vector<long> max_samples_per_frame_array(nt, 100 * patDims[0] * patDims[1] / reduction_factor);

    // Create kernel
    MDArray<4, double> kernel(kernelDims);
    kernel.Clear();
    WindowFunction window(1, 0.25);

    for (long kx = 0; kx < patDims[0]; kx++)
    {
        for (long ky = 0; ky < patDims[1]; ky++)
        {
            for (long i = 0; i < patDims[2]; i++)
            {
                for (long j = 0; j < patDims[2]; j++)
                {
                    long kx_abs_diff = MIN(kx, -kx + patDims[0]);
                    long ky_abs_diff = MIN(ky, -ky + patDims[1]);
                    std::array<long, 4> sub = {kx, ky, i, j};
                    const double kr = sqrt(pow(kx_abs_diff, 2) + pow(ky_abs_diff, 2));
                    const double kr_max = 4.0;
                    const double time_scaling = (0.1 / (0.1 + 10000 * abs(static_cast<double>(j - i))));
                    kernel[sub2ind(kernel.Strides(), sub)] = time_scaling * window(kr);
                }
            }
        }
    }

    // Faster version using a priority queue. Faster if w is sparse enough
    SparseKernel sparse_kernel(kernel, kernel.kthLargest(16));
    MDArray<kPhaseEncodeDims + 1, int> patApprox(patDims);
    double costApprox;
    SampleApproxBestCandidate<kPhaseEncodeDims>(sparse_kernel, max_samples_per_frame_array, totSamps, patApprox, costApprox);
    const MDArray<kPhaseEncodeDims + 1, double> deltaJApprox = ComputeDeltaJ<kPhaseEncodeDims>(kernel, patApprox);

    // Exact version, slower
    MDArray<kPhaseEncodeDims + 1, int> patExact(patDims);
    double costExact;
    MDArray<kPhaseEncodeDims + 1, double> deltaJExact(patDims);
    SampleExactBestCandidate<kPhaseEncodeDims>(kernel, max_samples_per_frame_array, totSamps, patExact, costExact, deltaJExact);

    printf("Cost with approximate best candidate: %f\n", costApprox);
    printf("Cost with exact best candidate: %f\n", costExact);

    assert(costApprox <= 250.0);
    assert(costExact <= 336.0);

    const MDArray<kPhaseEncodeDims + 1, double> deltaJExact2 = ComputeDeltaJ<kPhaseEncodeDims>(kernel, patExact);

    // Test that exact best candidate maintains the correct cost.
    for (long i = 0; i < deltaJExact.Length(); i++)
    {
        assert(abs(deltaJExact[i] - deltaJExact2[i]) < 1e-5);
    }

    const bool kDumpData = false;
    if (kDumpData)
    {
        kernel.Write("kernel");
        patApprox.Write("patApprox");
        patExact.Write("patExact");
        deltaJExact.Write("deltaJExact");
        deltaJExact2.Write("deltaJExact2");
        deltaJApprox.Write("deltaJApprox");
    }

    printf("TEST PASSED!\n");

    return 0;
}
