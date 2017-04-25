
//#include "../debug/debug.h"
//#include "../misc/misc.h"
//#include "../misc/misc.hpp"
#include "misc/debug.h"
#include "misc/misc.h"
#include <assert.h>
#include "sampleHeap.h"
#include "sample.h"
#include "misc_dd.h"

static int deltaJPrintCount = 1;

// Sample heap class
static inline int getParentIndex(int nodeIndex){
    return (nodeIndex + 1)/ 2 - 1;
}
static inline int getRightChildIndex(int nodeIndex){
    return 2*(nodeIndex + 1);
}
static inline int getLeftChildIndex(int nodeIndex){
    return 2*(nodeIndex + 1) - 1;
}

inline void swapSamples(SampleHeap *heap, int i1, int i2){
    swap(&heap->arr[i1], &heap->arr[i2], sizeof(Sample));
    heap->kt2idx[heap->arr[i1].kt] = i1;
    heap->kt2idx[heap->arr[i2].kt] = i2;
}

// Percolate down for pop
void percolateDown(SampleHeap *heap, int nodeIndex){
    assert(nodeIndex >= 0);
    // Percolate down
    int childIndexR, childIndexL;
    while(1) {
        childIndexL = getLeftChildIndex(nodeIndex); 
        childIndexR = getRightChildIndex(nodeIndex); 

        int smallest;
        if( childIndexL < (signed) heap->size && dJCompare(&heap->arr[childIndexL], &heap->arr[nodeIndex]) == -1 ){ // childL < node
            smallest = childIndexL;
        }else{
            smallest = nodeIndex;
        }
        if( childIndexR < (signed) heap->size && dJCompare(&heap->arr[childIndexR], &heap->arr[smallest]) == -1 ){ // childR < node
            smallest = childIndexR;
        }
        if( smallest == nodeIndex ){
            return;
        }else{
            swapSamples(heap, nodeIndex, smallest);
            nodeIndex = smallest;
        }
    } // end while loop
}

// Percolate up for push
inline void percolateUp(SampleHeap *heap, int nodeIndex){
    assert(nodeIndex >= 0);
    // Percolate up
    int parentIndex = getParentIndex(nodeIndex);
    while( nodeIndex != 0 && dJCompare(&heap->arr[nodeIndex], &heap->arr[parentIndex]) == -1 ){ // node < parent
        swapSamples(heap, nodeIndex, parentIndex);
        nodeIndex = parentIndex;
        parentIndex = getParentIndex(nodeIndex);
    }
}


// Increase key
void increaseKey(SampleHeap *heap, int nodeIndex, double delta){
    heap->arr[nodeIndex].dJ += delta;
    if( delta > 0 ){
        percolateDown(heap, nodeIndex);
    }else{
        percolateUp(heap, nodeIndex);
    }
}
   
// Push data onto the heap
void push(SampleHeap *heap, long kt_ind, double dJ){
    if( heap->size == heap->arr_size ){ /* arr.size() -> arr_size */
        debug_printf(DP_DEBUG3, "resizing...\n");
        heap->arr = realloc(heap->arr, sizeof(Sample) * heap->size * 2);
        heap->kt2idx = realloc(heap->arr, sizeof(Sample) * heap->size * 2);
        heap->arr_size *= 2;
    }
    int nodeIndex = heap->size;
    //debug_printf(DP_INFO, "Array size: %d\n", heap->arr_size);
    //debug_printf(DP_INFO, "Allocating %d for Sample pointer at %d\n", sizeof(Sample *), nodeIndex);
    //heap->arr[nodeIndex] = (Sample *) xmalloc(sizeof(Sample *)); // TODO: free this

    //debug_printf(DP_INFO, "Allocated %d for Sample pointer at %d\n", sizeof(Sample *), nodeIndex);
    heap->arr[nodeIndex].dJ = dJ;
    heap->arr[nodeIndex].kt = kt_ind;
    heap->kt2idx[kt_ind] = nodeIndex;
    heap->size++;

    //debug_printf(DP_INFO, "Percolate up...\n");
    percolateUp(heap, nodeIndex);
}
 
// Pop a sample from the heap
Sample pop(SampleHeap *heap){
    Sample res = heap->arr[0];
    heap->size--;
    swapSamples(heap, 0, heap->size);

    percolateDown(heap, 0);
    /*
    debug_printf(DP_INFO, "Validing heap after pop\n");
    validateHeapIndexArray(heap);
    */

    //assert(res != heap->arr[0]);
    return res;
}

void free_heap(SampleHeap *heap){
    debug_printf(DP_DEBUG3, "Destructor...\n");
    free(heap->arr);
    free(heap);
    debug_printf(DP_DEBUG3, "Destructor done\n");
}
     
SampleHeap *build_sample_heap(const long _ksp_size, const double *deltaJ){
    debug_printf(DP_DEBUG3, "Building heap...\n");

    SampleHeap *heap = (SampleHeap *) xmalloc(sizeof(SampleHeap));

    // k-space size including time, mostly for debugging
    heap->ksp_size = _ksp_size;
    // Array of samples to store the heap
    heap->arr = (Sample *) xmalloc(_ksp_size * sizeof(Sample));
    heap->arr_size = _ksp_size;
    // Initialize mapping from (k,t) to index in the array
    heap->kt2idx = (int *) xmalloc(_ksp_size * sizeof(int));
    for( long i = 0 ; i < _ksp_size ; ++i ){
        heap->kt2idx[i] = -1;
    }

    // Heap size = 0
    heap->size = 0;
    // Push in random order to break ties randomly
    long *perm = (long *) xmalloc(_ksp_size*sizeof(long));
    for( int i = 0 ; i < _ksp_size ; ++i )
        perm[i] = i;
    randomperm(_ksp_size, perm);
    
    for( int kt_ind = 0 ; kt_ind < _ksp_size ; kt_ind++ ){
        assert(perm[kt_ind] < _ksp_size && perm[kt_ind] >= 0);
        push(heap, perm[kt_ind], deltaJ[perm[kt_ind]]);
    }

    free(perm);
    debug_printf(DP_DEBUG3, "heap size: %d\n", heap->size);

    return heap;
}

static int validateHeapAux(Sample *h, const int i, const int sz){
    int il = 2*(i+1);
    int ir = 2*(i+1) - 1;
    if( il >= sz && ir >= sz ){
        return 1;
    }
    // Check parent / left child
    if( il < sz && h[il].dJ < h[i].dJ){
        return 0;
    }
    // Check parent / right child
    if( ir < sz && h[ir].dJ < h[i].dJ){
        return 0;
    }
    // Check left heap
    int leftIsHeap =  (il >= sz) || validateHeapAux(h, il, sz);
    // Check right heap
    int rightIsHeap = (ir >= sz) || validateHeapAux(h, ir, sz);
    return leftIsHeap && rightIsHeap;
}
int validateHeap(SampleHeap *heap){
    return validateHeapAux(heap->arr, 0, heap->size);
}

int validateHeapIndexArray(SampleHeap *heap){
    // Check heap kt2idx array
    for( unsigned long kt_ind = 0 ; kt_ind < heap->size ; kt_ind++ ){
        if( heap->kt2idx[kt_ind] != -1 && heap->arr[heap->kt2idx[kt_ind]].kt != (signed) kt_ind){
            debug_printf(DP_ERROR, "kt_ind: %lu, kt: %ld, size: %ld\n", kt_ind, heap->arr[heap->kt2idx[kt_ind]].kt, heap->size);
            debug_printf(DP_INFO, "kt2idx: ");
            for( unsigned long i = 0 ; i < heap->size ; i++ ){
                debug_printf(DP_INFO, "%d ", heap->kt2idx[i]);
            }
            debug_printf(DP_INFO, "\n");
            
            debug_printf(DP_INFO, "kt: ");
            for( unsigned long i = 0 ; i < heap->size ; i++ ){
                debug_printf(DP_INFO, "%d ", heap->arr[i].kt);
            }
            debug_printf(DP_INFO, "\n");
            assert(0);
            return 0;
        }
    }
    return 1;
}


void printHeap(SampleHeap *heap, const long dims[]){
    /*
    debug_printf(DP_INFO, "Indices = [\n");
    for( unsigned long i = 0 ; i < heap->size ; i++ ){
        debug_printf(DP_INFO, "%d -> %d, ... \n", 1+i, 1+heap.kt2idx[i]);
    }
    debug_printf(DP_INFO, "];\n");
    */
    
    for( unsigned long i = 0 ; i < heap->size ; i++ ){
        long kt_sub[3];
        idx2sub(3, dims, kt_sub, heap->arr[i].kt);
        debug_printf(DP_INFO, "deltaJ(%d,%d,%d,%d) = %f;\n", 
                    kt_sub[0]+1, kt_sub[1]+1, kt_sub[2]+1, deltaJPrintCount,
                    heap->arr[i].dJ);
    }
    deltaJPrintCount++;

    debug_printf(DP_INFO, "SampleHeap = [\n");
    for( unsigned long i = 0 ; i < heap->size ; i++ ){
        long kt_sub[3];
        idx2sub(3, dims, kt_sub, heap->arr[i].kt);
        debug_printf(DP_INFO, "%d %d %f; ...\n", 1+kt_sub[0], 1+kt_sub[1], heap->arr[i].dJ);
    }
    debug_printf(DP_INFO, "];\n");
    
    // Check heap index array
    validateHeapIndexArray(heap);
}

