#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <boost/program_options.hpp>

#include "debug/debug.h"
#include "dda_utils/dda_utils.h"
#include "misc/misc.h"
#include "misc/multind.h"


static int verbose;
using namespace std;

namespace po = boost::program_options;

/*
 * Use boost to process command line arguments
 */
static bool processCommandLine(int argc, char *argv[], string &wfile, string &patfile, string &deltaJfile, int &w_type){
    try{
        //po::options_description desc("Allowed options");
        po::options_description desc("Options");
        desc.add_options()
            ("help", "produce help message")
            ("w", po::value<string>(&wfile), "weighting file")
            ("pat", po::value<string>(&patfile), "pattern file")
            ("dJ", po::value<string>(&deltaJfile), "output deltaJ file")
            ("w_type", po::value<int>(&w_type), "1 => <w1,p1>, 2 => <w2,p2>")
            ("debug", po::value<int>(&debug_level), "debug level")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        // S = total samples
        if ( !vm.count("w") || !vm.count("pat") || !vm.count("dJ") ) {
            cout << desc << "\n";
            return false;
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
    string wfile, patfile, deltaJfile;
    int w_type = 3;
    if( !processCommandLine(argc, argv, wfile, patfile, deltaJfile, w_type) ){
        return 0;
    }

    // Read in w from wfile
    MDArray<double> *wmd = MDArray<double>::read_array(wfile);

    // Read in pattern
    MDArray<int> *pat = MDArray<int>::read_array(patfile);

    // Create output
    MDArray<double> *deltaJ = new MDArray<double>(pat->D, pat->dims);

    // Convert pattern to sample lists
    long Nt = pat->dims[PE_DIMS];
    long *samples[Nt];
    long Nsamps[Nt];
    find_samples(samples, Nsamps, pat->data, pat->dims, PE_DIMS);
    delete pat;
    
    // Compute Delta J
    if( w_type == 1 || w_type == 2 ){
        assert(wmd->dims[4] == 1);
        computeDeltaJ2(PE_DIMS, deltaJ->data, wmd, samples, Nsamps, w_type);
    }else{
        computeDeltaJ(PE_DIMS, deltaJ, wmd, samples, Nsamps);
    }

    // Write output file
    MDArray<double>::write_array(deltaJfile, deltaJ);

    // Clean up
    delete wmd;
    delete deltaJ;
    for( int t = 0 ; t < Nt ; t++ ){
        if( Nsamps[t] > 0 ){
            delete [] samples[t];
        }
    }
    return 0;
}


