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
#include "misc.hpp"
#include "dda_utils.h"
#include "debug.h"
#include "config.h"

// Faster version

// Boost
#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace std;

/* Write pattern. First line contains dims, each line contains a linear index of a sample location
 *
 * filename = file name
 * D        = number of dims in pat_dims
 * pat_dims = dimensions, {yres, zres, #frames}
 * pat      = pattern data
*/
void writePat(const string filename, const long D, const long pat_dims[], const int *pat){
    //debug_printf(DP_INFO, "Writing pattern...\n");
    long N = md_calc_size(D, pat_dims);
    FILE * pFile = fopen(filename.c_str(), "w");
            
    for( int i = 0 ; i < D ; i++ ){
        fprintf(pFile, "%ld ", pat_dims[i]);
    }
    fprintf(pFile, "\n");

    for( int i = 0 ; i < N ; i++ ){
        for( int j = 0 ; j < pat[i] ; j++ ){
            fprintf(pFile, "%d\n", i);
        }
    }
    fclose(pFile);
}

/*
 * Use boost to process command line arguments
 */
static bool processCommandLine(int argc, char *argv[], double &T, int &K, int &max_samples_per_frame, int &totSamps, int &exact, string &kernel_file, string &patfile){
    try{
        //po::options_description desc("Allowed options");
        po::options_description desc("Options");
        desc.add_options()
            ("help", "produce help message")
            ("max-per-frame", po::value<int>(&max_samples_per_frame), "set max samples per phase")
            ("T", po::value<double>(&T), "hard threshold for w")
            ("K", po::value<int>(&K), "set support of thresholded w")
            ("S", po::value<int>(&totSamps), "total samples")
            ("e", po::value<int>(&exact), "exact best candidate")
            ("w", po::value<string>(&kernel_file), "weighting file")
            ("pat", po::value<string>(&patfile), "pattern file")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        // S = total samples
        if ( !vm.count("S") || !vm.count("w") || !vm.count("pat") ) {
            cout << desc << "\n";
            return false;
        }
        // default max samples per phase
        if( !vm.count("max-per-frame") ){
            max_samples_per_frame = totSamps;
        }

    } // end try
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    catch(...)
    {
        std::cerr << "Unknown error!" << "\n";
        return false;
    }

    return true;
}

int main( int argc, char* argv[] )
{
    // Process arguments
    double T = 0;
    int exact, max_samples_per_frame, totSamps, K = 0;
    string kernel_file, patfile;
    if( !processCommandLine(argc, argv, T, K, max_samples_per_frame, totSamps, exact, kernel_file, patfile) ){
        return 0;
    }

    // Read in the kernel w from kernel_file
    MDArray<4, double> kernel(kernel_file);

    // Output dims
    const long nframes = kernel.dims[kPhaseEncodeDims];
    long pat_dims[4];
    md_select_dims(4, 7u, pat_dims, kernel.dims);
    long N = md_calc_size(3, pat_dims);

    vector<long> max_samples_per_frame_array(nframes, max_samples_per_frame);

    //cout << kernel.toString() << endl;
#if 0
    for( long i = 0 ; i < kernel.len ; ++i ){
        debug_printf(DP_INFO, "w[%d] = %f\n", i, kernel.data[i]);
    }
#endif

    //  Build sparse w
    debug_level = DP_ALL;

    SparseKernel sparse_kernel;
    if( !exact ){
        debug_printf(DP_DEBUG3, "building sparse w, K = %d\n", K);
        if( K ){
            sparse_kernel = sparsifyWToK(kernel, K);
        }else{
            sparse_kernel = sparsifyW(kernel, T);
        }
    }
    sparse_kernel.Print();

    // Generate pat
    double cost;
    MDArray<3, int> pat(pat_dims);
    pat.Clear();
    MDArray<3, double> deltaJ(pat_dims);

    if( exact ){
        // Exact version, slower
        exactBestCandidate(kernel, max_samples_per_frame_array, totSamps, pat, cost, deltaJ);
    }else{
        // Faster version using a heap. Faster if w is sparse enough
        approxBestCandidate(sparse_kernel, max_samples_per_frame_array, totSamps, pat, cost);
    }

    // Write out pattern
    debug_printf(DP_INFO, "Writing pattern\n");
    writePat(patfile, 3, pat_dims, pat.data);
    debug_printf(DP_INFO, "Done writing pattern\n");

    return 0;
}
