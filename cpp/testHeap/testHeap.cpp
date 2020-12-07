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
#include "sampleHeap.h"

static bool CheckHeap(const SampleHeap &heap)
{
    // Check heap kt2idx array
    for (unsigned long kt_ind = 0; kt_ind < heap.Size(); kt_ind++)
    {
        if (heap.getKt2idx(kt_ind) != -1 && heap.getArr(heap.getKt2idx(kt_ind)).getKTIndex() != kt_ind)
        {
            debug_printf(DP_ERROR, "kt_ind: %lu, kt: %ld, size: %ld\n", kt_ind,
                         heap.getArr(heap.getKt2idx(kt_ind)).getKTIndex(), heap.Size());
            return false;
        }
    }
    return true;
}

int main(int argc, char *argv[])
{
    // Dimensions
    const std::array<long, kPhaseEncodeDims + 1> pat_dims = {7, 9, 2};
    const long N = md_calc_size(pat_dims);

    // Random cost
    MDArray<3, double> deltaJ(pat_dims);
    for (int i = 0; i < deltaJ.Length(); i++)
    {
        deltaJ[i] = static_cast<double>(rand() % 1000);
    }

    SampleHeap heap(deltaJ);

    // Check the inverse index
    for (int i = 0; i < 500; i++)
    {
        assert(CheckHeap(heap));
        int idx = rand() % N;
        heap.increaseKey(idx, rand() % 100);
    }

    // Check that popping the samples returns them in sorted order and maintains the inverse index
    std::vector<double> costs;
    costs.reserve(heap.Size());
    while (heap.Size() > 0)
    {
        costs.push_back(heap.pop().dJ);
        assert(CheckHeap(heap));
    }
    assert(costs.size() == N);
    assert(std::is_sorted(costs.begin(), costs.end()));

    printf("TEST PASSED!\n");

    return 0;
}
