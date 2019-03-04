#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <boost/program_options.hpp>

#include "debug.h"
#include "dda_utils.h"
#include "misc.h"
#include "multind.h"

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
    MDArray<4, double> wmd(wfile);

    // Read in pattern
    MDArray<3, int> pat(patfile);

    // Convert pattern to sample lists
    long Nt = pat.dims[kPhaseEncodeDims];
    vector<vector<long> > samples = find_samples(pat);

    // Compute Delta J
    MDArray<3, double> deltaJ = computeDeltaJ(wmd, samples);

    // Write output file
    deltaJ.Write(deltaJfile);
    
    return 0;
}


