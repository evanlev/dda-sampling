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
#include "sampleHeap.h"

// Faster version

// Boost
using namespace std;

int main( int argc, char* argv[] )
{
    // Dimensions
    const long pat_dims[3] = {386, 192, 1};
    const long N = md_calc_size(3, pat_dims);

    //  Build sparse w
    debug_level = DP_ALL;

    // Random cost
    MDArray<double> deltaJ(3, pat_dims);
    for( int i = 0 ; i < N ; i++ ){
        deltaJ[i] = static_cast<double>(rand() % 1000);
    }

    SampleHeap heap(N, deltaJ);

    for( int i = 0 ; i < 500 ; i++ ){
        heap.Verify();
        int idx = rand() % N;
        heap.increaseKey(idx, rand() % 100);
    }

    

    return 0;
}

