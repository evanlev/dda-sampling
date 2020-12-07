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

using std::string;

namespace po = boost::program_options;

/*
 * Use boost to process command line arguments
 */
static bool processCommandLine(int argc, char *argv[], string &wfile, string &patfile, string &deltaJfile){
    try{
        po::options_description desc("Options");
        desc.add_options()
            ("help", "produce help message")
            ("w", po::value<string>(&wfile), "weighting file")
            ("pat", po::value<string>(&patfile), "pattern file")
            ("dJ", po::value<string>(&deltaJfile), "output deltaJ file")
            ("debug", po::value<int>(&debug_level), "debug level")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        // S = total samples
        if (!vm.count("w") || !vm.count("pat") || !vm.count("dJ"))
        {
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

int main(int argc, char *argv[])
{
    // Process arguments
    string wfile, mask_file, deltaJ_file;
    if (!processCommandLine(argc, argv, wfile, mask_file, deltaJ_file))
    {
        return 0;
    }

    // Read in w from wfile
    MDArray<kPhaseEncodeDims + 2, double> kernel(wfile);

    // Read in pattern
    MDArray<kPhaseEncodeDims + 1, int> mask(mask_file);

    // Compute Delta J
    MDArray<kPhaseEncodeDims + 1, double> deltaJ = ComputeDeltaJ<kPhaseEncodeDims>(kernel, mask);

    // Write output file
    deltaJ.Write(deltaJ_file);

    return 0;
}
