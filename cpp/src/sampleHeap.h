#ifndef SAMPLEHEAP_H
#define SAMPLEHEAP_H 1

#include "config.h"
#include "sample.h"
#include "mdarray.h"
#include <vector>

// Heap class allowing increasing arbitrary elements. This is a min-heap.
class SampleHeap {
    private:
        std::vector<Sample> arr;
        std::vector<int> kt2idx; // inverse index: kt2idx[kt] = index in array containing element at k-space index kt
        void copyTo(const int nodeIndex, const Sample &data);

    public:
        SampleHeap(const MDArray<kPhaseEncodeDims + 1, double> &deltaJ);

        unsigned long Size() const
        {
            return static_cast<unsigned long>(arr.size());
        }

        int getKt2idx(const int nodeIndex) const;

        void swapSamples(int i1, int i2);

        void percolateDown(int nodeIndex);

        void percolateUp(int nodeIndex);

        void increaseKey(int nodeIndex, double delta);

        void push(Sample data);

        Sample pop();

        Sample getArr(int idx) const;

        void Print(const long dims[]) const;

        bool Verify() const;
};

#endif

