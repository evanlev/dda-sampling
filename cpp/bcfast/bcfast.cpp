#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>

#include "mdarray.h"
#include "multind.h"
#include "misc.h"
#include "dda_utils.h"
#include "debug.h"
#include "config.h"

// Faster version

// Boost
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace
{
    /* Write pattern. First line contains dims, each line contains a linear index of a sample location
    *
    * filename = file name
    * D        = number of dims in pat_dims
    * pat_dims = dimensions, {yres, zres, #frames}
    * pat      = pattern data
    */
    template <size_t D>
    void writePat(const std::string& filename, const std::array<long, D>& pat_dims, int *pat)
    {
        debug_printf(DP_DEBUG4, "Writing pattern...\n");
        long N = md_calc_size(pat_dims);
        FILE *pFile = fopen(filename.c_str(), "w");

        for (int i = 0; i < D; i++)
        {
            fprintf(pFile, "%ld ", pat_dims[i]);
        }
        fprintf(pFile, "\n");

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < pat[i]; j++)
            {
                fprintf(pFile, "%d\n", i);
            }
        }
        fclose(pFile);
    }
}

/*
 * Use boost to process command line arguments
 */
static bool processCommandLine(int argc, char *argv[], double &T, int &K, int &maxSamplesPerFrame, int &totalSamples, int &exact, std::string &kernel_file, std::string &patfile)
{
    try
    {
        //po::options_description desc("Allowed options");
        po::options_description desc("Options");
        desc.add_options()("help", "produce help message")
            ("max-per-frame", po::value<int>(&maxSamplesPerFrame), "set max samples per phase")
            ("T", po::value<double>(&T), "hard threshold for w")
            ("K", po::value<int>(&K), "set support of thresholded w")
            ("S", po::value<int>(&totalSamples), "total samples")
            ("e", po::value<int>(&exact), "exact best candidate")
            ("w", po::value<std::string>(&kernel_file), "weighting file")
            ("pat", po::value<std::string>(&patfile), "pattern file");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        // S = total samples
        if (!vm.count("S") || !vm.count("w") || !vm.count("pat"))
        {
            std::cout << desc << "\n";
            return false;
        }
        // default max samples per phase
        if (!vm.count("max-per-frame"))
        {
            maxSamplesPerFrame = totalSamples;
        }

    } // end try
    catch (std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    catch (...)
    {
        std::cerr << "Unknown error!"
                  << "\n";
        return false;
    }

    return true;
}

int main(int argc, char *argv[])
{
    // Process arguments
    double T = 0;
    int useExact, maxSamplesPerFrame, totalSamples, K = 0;
    std::string kernel_file, patfile;
    if (!processCommandLine(argc, argv, T, K, maxSamplesPerFrame, totalSamples, useExact, kernel_file, patfile))
    {
        return 0;
    }

    // Read in the kernel w from kernel_file
    MDArray<kPhaseEncodeDims + 2, double> kernel(kernel_file);

    // Output dims
    const long numFrames = kernel.Dims()[kPhaseEncodeDims];
    std::array<long, kPhaseEncodeDims + 1> pat_dims;
    md_select_dims(4, 7u, pat_dims.data(), kernel.Dims().data()); // 111
    long N = md_calc_size(3, pat_dims.data());

    const std::vector<long> maxSamplesPerFrameVector(numFrames, maxSamplesPerFrame);

    //  Build sparse w
    if (K)
    {
        T = kernel.kthLargest(K);
    }

    // Generate pat
    double cost;
    MDArray<kPhaseEncodeDims + 1, int> pat(pat_dims);
    pat.Clear();
    MDArray<kPhaseEncodeDims + 1, double> deltaJ(pat_dims);

    if (useExact)
    {
        // Exact version, slower
        SampleExactBestCandidate<kPhaseEncodeDims>(kernel, maxSamplesPerFrameVector, totalSamples, pat, cost, deltaJ);
    }
    else
    {
        // Faster version using a heap. Faster if w is sparse enough
        SparseKernel sparseKernel(kernel, T);
        sparseKernel.Print();
        SampleApproxBestCandidate<kPhaseEncodeDims>(sparseKernel, maxSamplesPerFrameVector, totalSamples, pat, cost);
    }

    // Write out pattern
    debug_printf(DP_INFO, "Writing pattern\n");
    writePat(patfile, pat_dims, pat.Data());
    debug_printf(DP_INFO, "Done writing pattern\n");

    return 0;
}
