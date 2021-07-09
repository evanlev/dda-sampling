#ifndef __SAMPLEHEAP_H
#define __SAMPLEHEAP_H

#include "sample.h"

// Heap class allowing increasing arbitrary elements. This is a min heap
typedef struct SampleHeap_s {
    unsigned long arr_size;
    unsigned long size; 
    unsigned long ksp_size;
    Sample *arr; // TODO Double pointer?
    int *kt2idx; // kt2idx[kt] = element of arr with this sample
} SampleHeap;

extern void push(SampleHeap *heap, long kt_ind, double cplus);

extern void free_heap(SampleHeap *heap);

extern int validateHeap(SampleHeap *heap);
extern int validateHeapIndexArray(SampleHeap *heap);

extern Sample pop(SampleHeap *heap);

extern void swapSamples(SampleHeap *heap, int i1, int i2);

extern void increaseKey(SampleHeap *heap, int nodeIndex, double delta);


extern void percolateDown(SampleHeap *heap, int nodeIndex);

extern void percolateUp(SampleHeap *heap, int nodeIndex);

extern void printHeap(SampleHeap *heap, const long dims[]);

extern SampleHeap *build_sample_heap(const long _ksp_size, const double *cplus);

#endif // __SAMPLEHEAP_H
