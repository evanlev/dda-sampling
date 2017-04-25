#ifndef SAMPLEHEAP_H
#define SAMPLEHEAP_H

#include "sample.h"

// Heap class allowing increasing arbitrary elements. This is a min heap
typedef struct SampleHeap_s {
    unsigned long arr_size;
    unsigned long size; 
    unsigned long ksp_size;
    Sample *arr; // TODO Double pointer?
    int *kt2idx; // kt2idx[kt] = element of arr with this sample
} SampleHeap;

void push(SampleHeap *heap, long kt_ind, double cplus);

void free_heap(SampleHeap *heap);

int validateHeap(SampleHeap *heap);
int validateHeapIndexArray(SampleHeap *heap);

Sample pop(SampleHeap *heap);

void swapSamples(SampleHeap *heap, int i1, int i2);

void increaseKey(SampleHeap *heap, int nodeIndex, double delta);


void percolateDown(SampleHeap *heap, int nodeIndex);

void percolateUp(SampleHeap *heap, int nodeIndex);

void printHeap(SampleHeap *heap, const long dims[]);

SampleHeap *build_sample_heap(const long _ksp_size, const double *cplus);


#endif

