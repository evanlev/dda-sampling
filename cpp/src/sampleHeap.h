#ifndef SAMPLEHEAP_H
#define SAMPLEHEAP_H

#include "sample.h"
#include "mdarray.h"
#include <vector>

using std::vector;

// Heap class allowing increasing arbitrary elements. This is a min heap
class SampleHeap {
    private:
        std::vector<Sample> arr;
        std::vector<int> kt2idx; // inverse index: kt2idx[kt] = index in array containing element at k-space index kt
        inline void copyTo(const int nodeIndex, const Sample &data);

    public:
        SampleHeap(const MDArray<3, double> &deltaJ);

        unsigned long Size() const {
            return static_cast<unsigned long>(arr.size());
        }

        int getKt2idx(const int nodeIndex) const ;

        static inline int getParentIndex(int nodeIndex);

        static inline int getRightChildIndex(int nodeIndex);

        static inline int getLeftChildIndex(int nodeIndex);

        inline void swapSamples(int i1, int i2);

        inline void percolateDown(int nodeIndex);

        inline void percolateUp(int nodeIndex);

        void increaseKey(int nodeIndex, double delta);

        void push(Sample data);

        Sample pop();

        Sample getArr(int idx) const;

        void Print(const long dims[]) const;

        bool Verify() const;
};

#endif

