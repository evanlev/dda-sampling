#ifndef SAMPLEHEAP_H
#define SAMPLEHEAP_H

#include "sample.h"
#include <vector>

using std::vector;
// Heap class allowing increasing arbitrary elements. This is a min heap
class SampleHeap {
    private:
        unsigned long size; 
        unsigned long ksp_size;
        std::vector<Sample *> arr;
        std::vector<int> kt2idx; // kt2idx[kt] = element of arr with this sample
    public:
        SampleHeap(const long ksp_size, const double *deltaJ);
        ~SampleHeap();

        unsigned long getSize(){
            return size;
        }

        int getKt2idx(int nodeIndex);

        static inline int getParentIndex(int nodeIndex);

        static inline int getRightChildIndex(int nodeIndex);

        static inline int getLeftChildIndex(int nodeIndex);

        inline void swapSamples(int i1, int i2);

        inline void percolateDown(int nodeIndex);

        inline void percolateUp(int nodeIndex);

        void increaseKey(int nodeIndex, double delta);

        inline void copyTo(int nodeIndex, Sample *data);

        void push(Sample *data);

        Sample * pop();

        Sample * getArr(int idx);
};
/*
*/

void printHeap(SampleHeap *heap, const long dims[]);

#endif

