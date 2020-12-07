#include "sampleHeap.h"

#include "debug.h"
#include "misc.h"
#include "misc.hpp"
#include "sample.h"

#include <algorithm> // std::random_shuffle

static int s_deltaJPrintCount = 1;

// Sample heap class
namespace
{
    int GetParentIndex(int nodeIndex)
    {
        return (nodeIndex + 1) / 2 - 1;
    }
    int GetRightChildIndex(int nodeIndex)
    {
        return 2 * (nodeIndex + 1);
    }
    int GetLeftChildIndex(int nodeIndex)
    {
        return 2 * (nodeIndex + 1) - 1;
    }
}

void SampleHeap::swapSamples(int i1, int i2)
{
    std::swap(arr[i1], arr[i2]);
    kt2idx[arr[i1].getKTIndex()] = i1;
    kt2idx[arr[i2].getKTIndex()] = i2;
}

// Percolate down for pop
void SampleHeap::percolateDown(int nodeIndex)
{
    assert(nodeIndex >= 0);
    // Percolate down
    int childIndexR, childIndexL;
    while (1)
    {
        childIndexL = GetLeftChildIndex(nodeIndex);
        childIndexR = GetRightChildIndex(nodeIndex);

        int smallest;
        if (childIndexL < Size() && arr[childIndexL].dJ < arr[nodeIndex].dJ)
        {
            smallest = childIndexL;
        }
        else
        {
            smallest = nodeIndex;
        }
        if (childIndexR < Size() && arr[childIndexR].dJ < arr[smallest].dJ)
        {
            smallest = childIndexR;
        }
        if (smallest == nodeIndex)
        {
            return;
        }
        else
        {
            swapSamples(nodeIndex, smallest);
            nodeIndex = smallest;
        }
    } // end while loop
}

// Percolate up for push
void SampleHeap::percolateUp(int nodeIndex)
{
    assert(nodeIndex >= 0);
    // Percolate up
    int parentIndex = GetParentIndex(nodeIndex);
    while (nodeIndex != 0 && arr[nodeIndex].dJ < arr[parentIndex].dJ)
    {
        swapSamples(nodeIndex, parentIndex);
        nodeIndex = parentIndex;
        parentIndex = GetParentIndex(nodeIndex);
    }
}

int SampleHeap::getKt2idx(const int nodeIndex) const
{
    return kt2idx[nodeIndex];
}

Sample SampleHeap::getArr(int idx) const
{
    return arr[idx];
}

// Increase key
void SampleHeap::increaseKey(int nodeIndex, double delta)
{
    arr[nodeIndex].dJ += delta;
    if (delta > 0)
    {
        percolateDown(nodeIndex);
    }
    else
    {
        percolateUp(nodeIndex);
    }
}

// Push data onto the heap
void SampleHeap::push(Sample data)
{
    arr.push_back(data);
    kt2idx[data.getKTIndex()] = arr.size() - 1;
    percolateUp(arr.size() - 1);
}

// Pop a sample from the heap
Sample SampleHeap::pop()
{
    Sample result = arr[0];

    arr[0] = arr.back();

    kt2idx[arr[0].getKTIndex()] = 0;
    kt2idx[result.getKTIndex()] = -1;

    arr.pop_back();
    percolateDown(0);

    return result;
}

SampleHeap::SampleHeap(const MDArray<kPhaseEncodeDims + 1, double> &deltaJ)
{
    debug_printf(DP_DEBUG3, "Building heap, size %d...\n", deltaJ.Length());
    kt2idx.resize(deltaJ.Length());

    // Push in random order to break ties randomly
    std::vector<int> perm(deltaJ.Length());
    for (int i = 0; i < deltaJ.Length(); ++i)
    {
        perm[i] = i;
    }
    std::random_shuffle(perm.begin(), perm.end());
    for (int kt_ind = 0; kt_ind < deltaJ.Length(); kt_ind++)
    {
        push(Sample(perm[kt_ind], deltaJ[perm[kt_ind]]));
    }
    debug_printf(DP_DEBUG3, "heap size: %d\n", Size());
}

void SampleHeap::Print(const long dims[]) const
{
    for (unsigned long i = 0; i < this->Size(); i++)
    {
        long kt_sub[kPhaseEncodeDims+1];
        ind2sub(kPhaseEncodeDims+1, dims, kt_sub, this->getArr(i).getKTIndex());
        debug_printf(DP_INFO, "deltaJ(%d,%d,%d) = %f;\n",
                     kt_sub[0] + 1, kt_sub[1] + 1, s_deltaJPrintCount,
                     this->getArr(i).dJ);
    }
    s_deltaJPrintCount++;

    debug_printf(DP_INFO, "SampleHeap = [\n");
    for (unsigned long i = 0; i < this->Size(); i++)
    {
        long kt_sub[kPhaseEncodeDims+1];
        ind2sub(kPhaseEncodeDims+1, dims, kt_sub, this->getArr(i).getKTIndex());
        debug_printf(DP_INFO, "%d %d %f; ...\n", 1 + kt_sub[0], 1 + kt_sub[1], this->getArr(i).dJ);
    }
    debug_printf(DP_INFO, "];\n");
}
