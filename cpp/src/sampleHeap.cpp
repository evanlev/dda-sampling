#include "debug.h"
#include "misc.h"
#include "misc.hpp"
#include "sampleHeap.h"
#include "sample.h"

#include <algorithm>    // std::random_shuffle

static int deltaJPrintCount = 1;

// Sample heap class
inline int SampleHeap::getParentIndex(int nodeIndex){
    return (nodeIndex + 1)/ 2 - 1;
}
inline int SampleHeap::getRightChildIndex(int nodeIndex){
    return 2*(nodeIndex + 1);
}
inline int SampleHeap::getLeftChildIndex(int nodeIndex){
    return 2*(nodeIndex + 1) - 1;
}
inline void SampleHeap::swapSamples(int i1, int i2){
    //swap(&arr[i1], &arr[i2], sizeof(Sample *));
    std::swap(arr[i1], arr[i2]);
    kt2idx[arr[i1].getIndex()] = i1;
    kt2idx[arr[i2].getIndex()] = i2;
}

// Percolate down for pop
inline void SampleHeap::percolateDown(int nodeIndex){
    assert(nodeIndex >= 0);
    // Percolate down
    int childIndexR, childIndexL;
    while(1) {
        childIndexL = getLeftChildIndex(nodeIndex); 
        childIndexR = getRightChildIndex(nodeIndex); 

        int smallest;
        if( childIndexL < size && arr[childIndexL] < arr[nodeIndex] ){
            smallest = childIndexL;
        }else{
            smallest = nodeIndex;
        }
        if( childIndexR < size && arr[childIndexR] < arr[smallest] ){
            smallest = childIndexR;
        }
        if( smallest == nodeIndex ){
            return;
        }else{
            swapSamples(nodeIndex, smallest);
            nodeIndex = smallest;
        }
    } // end while loop
}

// Percolate up for push
inline void SampleHeap::percolateUp(int nodeIndex){
    assert(nodeIndex >= 0);
    // Percolate up
    int parentIndex = getParentIndex(nodeIndex);
    //debug_printf(DP_INFO, "Node index: %d, parent: %d\n", nodeIndex, parentIndex);
    while( nodeIndex != 0 && arr[nodeIndex] < arr[parentIndex] ){
        //debug_printf(DP_INFO, "Swap...\n");
        swapSamples(nodeIndex, parentIndex);
        //debug_printf(DP_INFO, "Done Swap...\n");
        nodeIndex = parentIndex;
        parentIndex = getParentIndex(nodeIndex);
    }
}

int SampleHeap::getKt2idx(int nodeIndex) const {
    return kt2idx[nodeIndex];
}

Sample SampleHeap::getArr(int idx) const {
    return arr[idx];
}

// Increase key
void SampleHeap::increaseKey(int nodeIndex, double delta){
    arr[nodeIndex].dJ += delta;
    if( delta > 0 ){
        percolateDown(nodeIndex);
    }else{
        percolateUp(nodeIndex);
    }
}
    
inline void SampleHeap::copyTo(const int nodeIndex, const Sample &data){
    if( nodeIndex >= size || nodeIndex < 0 ){
        debug_printf(DP_INFO, "size: %d, nodeIndex: %d\n", size, nodeIndex);
        assert(0);
    }
    arr[nodeIndex] = data;
    kt2idx[data.getIndex()] = nodeIndex;
}

// Push data onto the heap
void SampleHeap::push(Sample data){

    if( size == arr.size() ){
        debug_printf(DP_INFO, "resizing...\n");
        arr.resize(size * 2);
        kt2idx.resize(size * 2);
    }
    int nodeIndex = size++;
    copyTo(nodeIndex, data); // arr[nodeIndex] = data

    percolateUp(nodeIndex);
}

// Pop a sample from the heap
Sample SampleHeap::pop(){
    Sample res = arr[0];
    size--;
    copyTo(0, arr[size]);

    percolateDown(0);

    return res;
}

SampleHeap::~SampleHeap(){
    debug_printf(DP_DEBUG3, "Destructor...\n");
    long size = getSize();
    debug_printf(DP_DEBUG3, "Destructor done\n");
}


SampleHeap::SampleHeap(const long _ksp_size, const MDArray<double> &deltaJ){
    debug_printf(DP_DEBUG3, "Building heap, size %d...\n", _ksp_size);

    // k-space size including time, mostly for debugging
    ksp_size = _ksp_size;
    // Array of samples to store the heap
    arr.resize(ksp_size);
    // Initialize mapping from (k,t) to index in the array
    kt2idx.resize(ksp_size);
    //for( std::vector<long>::iterator iter = kt2idx.begin() ; iter < kt2idx.end() ; ++iter ){
    for( long i = 0 ; i < (long) kt2idx.size() ; ++i ){
        //*iter = -1;
        kt2idx[i] = -1;
    }

    // Heap size = 0
    size = 0;
    // Push in random order to break ties randomly
    vector<int> perm(ksp_size);
    for( int i = 0 ; i < ksp_size ; ++i )
        perm[i] = i;
    std::random_shuffle(perm.begin(), perm.end());
    for( int kt_ind = 0 ; kt_ind < ksp_size ; kt_ind++ ){
        assert(perm[kt_ind] < ksp_size && perm[kt_ind] >= 0);
        push(Sample(perm[kt_ind], deltaJ[perm[kt_ind]]));
    }
    debug_printf(DP_DEBUG3, "heap size: %d\n", getSize());
}

void SampleHeap::Print(const long dims[]) const {
    /*
    debug_printf(DP_INFO, "Indices = [\n");
    for( unsigned long i = 0 ; i < heap.getSize() ; i++ ){
        debug_printf(DP_INFO, "%d -> %d, ... \n", 1+i, 1+heap.kt2idx[i]);
    }
    debug_printf(DP_INFO, "];\n");
    */
    
    for( unsigned long i = 0 ; i < this->getSize() ; i++ ){
        long kt_sub[3];
        ind2sub<long>(3, dims, kt_sub, this->getArr(i).getIndex());
        debug_printf(DP_INFO, "deltaJ(%d,%d,%d) = %f;\n", 
                    kt_sub[0]+1, kt_sub[1]+1, deltaJPrintCount,
                    this->getArr(i).dJ);
    }
    deltaJPrintCount++;

    debug_printf(DP_INFO, "SampleHeap = [\n");
    for( unsigned long i = 0 ; i < this->getSize() ; i++ ){
        long kt_sub[3];
        ind2sub<long>(3, dims, kt_sub, this->getArr(i).getIndex());
        debug_printf(DP_INFO, "%d %d %f; ...\n", 1+kt_sub[0], 1+kt_sub[1], this->getArr(i).dJ);
    }
    debug_printf(DP_INFO, "];\n");
}

bool SampleHeap::Verify() const {
    // Check heap kt2idx array
    for( unsigned long kt_ind = 0 ; kt_ind < this->getSize() ; kt_ind++ ){
        if(this->getKt2idx(kt_ind) != -1 && this->getArr(this->getKt2idx(kt_ind)).getIndex() != kt_ind){
            debug_printf(DP_ERROR, "kt_ind: %lu, kt: %ld, size: %ld\n", kt_ind, 
                this->getArr(this->getKt2idx(kt_ind)).getIndex(), this->getSize());
            return false;
        }
    }
    return true;
}


