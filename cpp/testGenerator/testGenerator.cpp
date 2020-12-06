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
#include "misc.hpp"
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
    const std::array<long, kPhaseEncodeDims + 1> pat_dims = {50, 60, 1};
    const long reduction_factor = 12;
    const long nt = pat_dims[kPhaseEncodeDims];

    const std::array<long, kPhaseEncodeDims + 2> ker_dims = {pat_dims[0], pat_dims[1], pat_dims[2], pat_dims[2]};
    const long N = md_calc_size(pat_dims);
    const long totSamps = N / reduction_factor;
    const vector<long> max_samples_per_frame_array(nt, 100 * pat_dims[0] * pat_dims[1] / reduction_factor);

    // Create kernel
    MDArray<4, double> kernel(ker_dims);
    kernel.Clear();
    WindowFunction window(1, 0.25);

    for (long kx = 0; kx < pat_dims[0]; kx++)
    {
        for (long ky = 0; ky < pat_dims[1]; ky++)
        {
            for (long i = 0; i < pat_dims[2]; i++)
            {
                for (long j = 0; j < pat_dims[2]; j++)
                {
                    long kx_abs_diff = MIN(kx, -kx + pat_dims[0]);
                    long ky_abs_diff = MIN(ky, -ky + pat_dims[1]);
                    long sub[4] = {kx, ky, i, j};
                    const double kr = sqrt(pow(kx_abs_diff, 2) + pow(ky_abs_diff, 2));
                    const double kr_max = 4.0;
                    const double time_scaling = (0.1 / (0.1 + 10000 * abs(static_cast<double>(j - i))));
                    kernel[sub2ind(4, kernel.m_strs, sub)] = time_scaling * window(kr);
                }
            }
        }
    }

    // Faster version using a heap. Faster if w is sparse enough
    SparseKernel sparse_kernel(kernel, kernel.kthLargest(16));
    MDArray<3, int> pat_approx(pat_dims);
    double cost_approx;
    approxBestCandidate(sparse_kernel, max_samples_per_frame_array, totSamps, pat_approx, cost_approx);
    MDArray<3, double> deltaJ_approx = computeDeltaJ(kernel, pat_approx);

    // Exact version, slower
    MDArray<3, int> pat_exact(pat_dims);
    double cost_exact;
    MDArray<3, double> deltaJ_exact(pat_dims);
    exactBestCandidate(kernel, max_samples_per_frame_array, totSamps, pat_exact, cost_exact, deltaJ_exact);

    printf("Cost with approximate best candidate: %f\n", cost_approx);
    printf("Cost with exact best candidate: %f\n", cost_exact);

    assert(cost_approx <= 250.0);
    assert(cost_exact <= 336.0);

    MDArray<3, double> deltaJ_exact2 = computeDeltaJ(kernel, pat_exact);

    // Test that exact best candidate maintains the correct cost.
    for (long i = 0; i < deltaJ_exact.Length(); i++)
    {
        assert(abs(deltaJ_exact[i] - deltaJ_exact2[i]) < 1e-5);
    }

    const bool kDumpData = false;
    if (kDumpData)
    {
        kernel.Write("kernel");
        pat_approx.Write("pat_approx");
        pat_exact.Write("pat_exact");
        deltaJ_exact.Write("dJ_exact");
        deltaJ_exact2.Write("dJ_exact2");
        deltaJ_approx.Write("dJ_approx");
    }

    printf("TEST PASSED!\n");

    return 0;
}
