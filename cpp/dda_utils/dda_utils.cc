#include "../misc/misc.h"
#include "../misc/misc.hpp"
#include "../misc/multind.h"
#include "../misc/mdarray.h"
#include "../debug/debug.h"

#include "dda_utils.h"
#include "sampleHeap.h"
//#include "sample.h"

#include <algorithm>
#include <vector>
#include <limits>

//#define DEBUG

using std::vector;
            
const static double Inf = std::numeric_limits<double>::max();


void printWsp(const SparseW *wsp){
    debug_printf(DP_INFO, "sparse w: \n");
    long D = wsp->D;
    int nt = wsp->dims[D];
    int nr = wsp->dims[D+2];
    for( int r = 0  ; r < nr ; r++ )
    for( int t2 = 0 ; t2 < nt ; t2++ )
    for( int t1 = 0 ; t1 < nt ; t1++ )
    for( int j = 0 ; j < wsp->len_w[t1][t2][r] ; j++ ){
        debug_printf(DP_INFO, "w(");
        for( int i = 0 ; i < D ; i++ ){
            debug_printf(DP_INFO, "%d, ", wsp->kw[t1][t2][r][j*D + i]);
        }
        debug_printf(DP_INFO,"%d, %d, %d) = %f\n", t1, t2, r, wsp->w[t1][t2][r][j]);
    }
}

void sparsifyWToK(const long D, SparseW *wsp, MDArray<double> *wmd, const long k){
    debug_print_dims(DP_INFO, D+1, wmd->dims);
    double T = kthLargestDouble(md_calc_size(D, wmd->dims) * wmd->dims[D] * wmd->dims[D], wmd->data, k);
    sparsifyW(D, wsp, wmd, T);
}

void sparsifyW(const long D, SparseW *wsp, MDArray<double> *wmd, const double T){

    debug_printf(DP_DEBUG3, "sparsifying w\n");
    // dims
    long ksize = md_calc_size(D, wmd->dims);
    long nt = wmd->dims[D];
    long nr = wmd->dims[D+2];
    long strs[D+3];
    md_calc_strides(D+3, strs, wmd->dims, 1);
    wsp->D = D;
    myAssert(D < MAX_DIMS, "D < max_dims");
    memcpy(wsp->dims, wmd->dims, wmd->D*sizeof(long));

    // set wsp->w00tt
    for( long t = 0 ; t < nt ; t++ ){
        wsp->w00tt[t] = 0;
        for( long r = 0 ; r < nr ; r++ )
            wsp->w00tt[t] = wsp->w00tt[t] + wmd->data[t * strs[D] + t*strs[D+1] + r*strs[D+2]];
    }

    // Initialize len w
    for( long r = 0 ; r < nr ; r++ )
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        wsp->len_w[t1][t2][r] = 0;
    }
    // Allocate w k-v pairs
    //debug_printf(DP_INFO, "allocating..\n");
    for( long r = 0 ; r < nr ; r++ )
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        // TODO: don't need to allocate this
        wsp->w[t1][t2][r]  = (double *) xmalloc(ksize*sizeof(double));
        wsp->kw[t1][t2][r] = (long *) xmalloc(D*ksize*sizeof(long));
    }
    // Assign k-v pairs for w
    long supp = 0;
    for( long r = 0 ; r < nr ; r++ )
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        for( long k = 0 ; k < ksize ; k++ ){
            long ksub[D];
            ind2sub<long>(D, wmd->dims, ksub, k);
            long kt1t2r_ind = k + t1*strs[D] + t2*strs[D+1] + r*strs[D+2];
            if( wmd->data[kt1t2r_ind] > T ){
                // append (k) to w_sp[t1][t2]
                /*
                debug_printf(DP_INFO, "Appending (%d %d %d %d), %f, len w: %d\n",
                            ksub[0], ksub[1], t1, t2, wmd->data[kt1t2r_ind],
                            wsp->len_w[t1][t2]);
                */
                supp++;
                wsp->w[t1][t2][r][wsp->len_w[t1][t2][r]] = wmd->data[kt1t2r_ind];
                memcpy(&wsp->kw[t1][t2][r][D*wsp->len_w[t1][t2][r]], ksub, D*sizeof(long));
                wsp->len_w[t1][t2][r]++;
            }
        }
    }
    debug_printf(DP_INFO, "w support: %d/%d = %f\n", supp, nt*nt*ksize, (float) supp / (float) (nt*nt*ksize));
}

void printPat(const int *mask, const long D, const long dims[]){
    static int maskPrintCount = 1;

    debug_printf(DP_INFO, "mask(:,%d) = [", maskPrintCount++);
    long N  = md_calc_size(D+1, dims);
    for( long i = 0 ; i < N ; i++ ){
        debug_printf(DP_INFO, "%d ", mask[i]);
    }
    debug_printf(DP_INFO, "];\n");
}

/*
 * Best candidate sampling but use a heap to store J
 */
void approxBestCandidate(const int D, double *cost, double *deltaJ, 
                                 int *mask, const SparseW *wsp,
                                 const long maxSamps[], 
                                 const long totSamps){
    //debug_level = DP_ALL;
    debug_printf(DP_INFO, "Best candidate version 2\n");
    
    long nt = wsp->dims[D];
    long nr = wsp->dims[D+2];
    long strs[D+2];
    md_calc_strides(D+2, strs, wsp->dims, 1);
    long ksize = md_calc_size(D, wsp->dims);
    long csize = nt*ksize;

    // Total cost
    *cost = 0;
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = wsp->w00tt[t];
        md_fill(D, wsp->dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

#if 0
    // DEBUG
    for( int kt_ind = 0 ; kt_ind < csize ; kt_ind++ ){
        deltaJ[kt_ind] = rand() % 100;
    }
#endif

    // Clear out the mask
    memset(mask, 0, csize*sizeof(int));

    // Current number of samples
    long nSampsCurrent[nt];
    memset(nSampsCurrent, 0, nt*sizeof(long));

    // Optional constraint on max samples in each pattern
    int n_full = 0;
    int is_full[nt];
    for( long t = 0 ; t < nt ; t++ )
        is_full[t] = maxSamps[t] == 0;

    // Heapify delta J, later can use the fact that deltaJ is constant to speed this up
    debug_printf(DP_INFO, "Construct heap...\n");
    SampleHeap *heap = new SampleHeap(csize, deltaJ);  

    // Best candidate selection loop
    debug_printf(DP_INFO, "Main loop...\n");
    for( long i = 0 ; i < totSamps ; i++ ){
        //printHeap(heap, wsp->dims);

        // TODO: pop or have an option if we can reacquire
        Sample *snew = heap->getArr(0);
        
        long k_new_sub[D+1];
        ind2sub<long>(D+1, wsp->dims, k_new_sub, snew->getKT());
        long t_new = k_new_sub[D];

        debug_printf(DP_DEBUG4, "Sampling (%d %d)\n", 1+k_new_sub[0], 1+k_new_sub[1]);

        mask[snew->getKT()]++;
        nSampsCurrent[t_new]++;

        // --- Update total cost
        *cost += deltaJ[snew->getKT()];

        // --- Update insert cost (deltaJ), even for frames that are full
        for( long r = 0 ; r < nr ; r++ )
        for( long t = 0 ; t < nt ; t++ ){

            // w( Delta k, t, t' )
            for( long j = 0 ; j < wsp->len_w[t][t_new][r] ; j++ ){
                // tmp = k' + (Delta k)_j
                long tmp[D+1];
                add<long>(D, tmp, k_new_sub, &wsp->kw[t][t_new][r][j*D]);
                mod<long>(D, tmp, tmp, wsp->dims);
                tmp[D] = t;
                long i1 = sub2ind(D+1,strs,tmp);
                heap->increaseKey(heap->getKt2idx(i1), wsp->w[t][t_new][r][j]);
            }

            // w( Delta k, t', t)
             for( long j = 0 ; j < wsp->len_w[t_new][t][r] ; j++ ){
                // tmp = k' - (Delta k)_j
                long tmp[D+1];
                sub<long>(D, tmp, k_new_sub, &wsp->kw[t_new][t][r][j*D]);
                mod<long>(D, tmp, tmp, wsp->dims);
                tmp[D] = t;
                long i1 = sub2ind(D+1,strs,tmp);
                // deltaJ[k' + (Delta k)_j] += w[t'][t][j]
                heap->increaseKey(heap->getKt2idx(i1), wsp->w[t_new][t][r][j]);
            }
        } // end t loop
    
        //printPat(mask, D, wsp->dims);
        //printHeap(heap, wsp->dims);
        //debug_printf(DP_DEBUG4, "Done updating heap...\n");

        // --- Check if this pattern is full
        if( nSampsCurrent[t_new] == maxSamps[t_new] ){
            is_full[t_new] = 1;
            n_full++;
            debug_printf(DP_DEBUG3, "%d is full\n", t_new);
            // TODO: do something about this
            //md_fill(D, wsp->dims, &deltaJ[t_new*ksize], &Inf, sizeof(double));
        }

        // --- Break if all patterns are full
        if( n_full == nt ){
            debug_printf(DP_DEBUG3, "All frames full, breaking at %d/%d\n", i+1, totSamps);
            break;
        }

        // --- Progress
        if( totSamps > 10 && totSamps % (totSamps / 10)  == 0 ){
            debug_printf(DP_INFO, "\n%d%", (int) (100.0 * (float) i / (float) totSamps));
        }
    } // end i loop
    
    delete heap;
    debug_printf(DP_INFO, "Done!\n");
}

/*
 * Best candidate sampling
 * compute deltaJ, cost change associated with adding a sample 
 *
 * INPUTS:
 *   D        = number of dimensions of k-space (usually 1 or 2)
 *   w        = weighting cost function
 *   maxSamps = max samples per pattern
 *   totSamps = total number of samples
 *
 * OUTPUTS:
 *   cost    = final cost
 *   deltaJ   = final insertion cost
 *   mask    = final sampling pattern
 *   samples = sample list
 */
void exactBestCandidate(const int D, double *cost, double *deltaJ, 
            int *mask, const MDArray<double> *w, 
            long *samples[], const long maxSamps[], 
            const long totSamps){
    //assert(D == 2);
    const long REIM_DIM = 4;
    long nt = w->dims[D];
    long ksize = md_calc_size(D, w->dims);
    long csize = nt*ksize;

    // Total cost
    *cost = 0;
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = w->data[t*w->strs[D] + t*w->strs[D+1]];
        md_fill(D, w->dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

    // Clear out the mask
    memset(mask, 0, csize*sizeof(int));

    // Current number of samples
    long nSampsCurrent[nt];
    memset(nSampsCurrent, 0, nt*sizeof(long));

    // Optional constraint on max samples in each pattern
    int n_full = 0;
    int is_full[nt];
    for( long t = 0 ; t < nt ; t++ )
        is_full[t] = maxSamps[t] == 0;

    // Best candidate selection loop
    for( long i = 0 ; i < totSamps ; i++ ){
        // k_new = argmin deltaJ(:,:,m) 
        long kt_new_ind = argmin<double>(csize, deltaJ);
        long k_new_sub[D+1];
        ind2sub<long>(D+1, w->dims, k_new_sub, kt_new_ind);
        long t_new = k_new_sub[D];
        myAssert(!is_full[t_new], "argmin should not have selected a full frame!");

        // --- Add sample: TODO one function
        // Increment mask
        mask[kt_new_ind]++;
        // Insert into sample list
        memcpy(&samples[t_new][D*nSampsCurrent[t_new]], k_new_sub, D*sizeof(long));
        nSampsCurrent[t_new]++;

        // --- Update total cost
        *cost += deltaJ[kt_new_ind];

        // --- Update insert cost (deltaJ), even for frames that are full
        for( long t = 0 ; t < nt ; t++ ){
        for( long j = t*ksize ; j < (t+1)*ksize ; j++ ){

            // diff = (k_new_sub - sj, t_new, t)
            long sj[D+1];
            ind2sub<long>(D+1, w->dims, sj, j);

            long diff[D+3]; 
            sub<long>(D, diff, k_new_sub, sj);
            mod<long>(D, diff, diff, w->dims);
            diff[D+0] = t_new;
            diff[D+1] = t;
           
            // deltaJ[j] += w[snew - sj, t_new, tj]
            deltaJ[j] += w->data[sub2ind(D+2, w->strs, diff)];

            // diff = (sj - k_new_sub, t, t_new)
            sub<long>(D, diff, sj, k_new_sub);
            mod<long>(D, diff, diff, w->dims);
            diff[D+0] = t;
            diff[D+1] = t_new;

            // deltaJ[j] += w[sj - snew, tj, t_new]
            deltaJ[j] += w->data[sub2ind(D+2, w->strs, diff)];
    
            // Part 2
            if( w->dims[REIM_DIM] == 2){
                // diff = (sj + k_new_sub, t, t_new)
                add<long>(D, diff, sj, k_new_sub);
                mod<long>(D, diff, diff, w->dims);
                diff[D+0] = t;
                diff[D+1] = t_new;
                diff[D+2] = 1;

                // deltaJ[j] += w[sj + snew, tj, t_new, 1]
                long i1 = sub2ind(D+3, w->strs, diff);
                deltaJ[j] += w->data[i1];
                // deltaJ[j] += w[sj + snew, t_new, tj, 1]
                diff[D+0] = t_new;
                diff[D+1] = t;
                i1 = sub2ind(D+3, w->strs, diff);
                deltaJ[j] += w->data[i1];
            }
        } // end j loop
        } // end m loop

        // --- Check if this pattern is full
        if( nSampsCurrent[t_new] == maxSamps[t_new] ){
            is_full[t_new] = 1;
            n_full++;
            debug_printf(DP_DEBUG3, "%d is full\n", t_new);
            md_fill(D, w->dims, &deltaJ[t_new*ksize], &Inf, sizeof(double));
        }

        // --- Break if all patterns are full
        if( n_full == nt ){
            debug_printf(DP_DEBUG3, "All frames full, breaking at %d/%d\n", i, totSamps);
            break;
        }

        // --- Progress
        if( totSamps % (totSamps / 10)  == 0 ){
            debug_printf(DP_INFO, "\r%d%", (int) (100.0 * (float) i / (float) totSamps));
        }
    } // end i loop
}

void computeDeltaJ(const int D, MDArray<double> *deltaJ, const MDArray<double> *w, long *samples[], const long Nsamps[]){
    debug_printf(DP_DEBUG1, "Computing DeltaJ...\n");

    computeDeltaJ2(D, deltaJ->data, w, samples, Nsamps, 1);
    int Nw = w->dims[D+2];
    if( Nw == 2){
        debug_printf(DP_DEBUG1, "Computing second term in DeltaJ for real-valued images...\n");
        // Compute the second component
        double *deltaJ2 = new double[deltaJ->len];
        computeDeltaJ2(D, deltaJ2, w, samples, Nsamps, 2);

        // Add to the cost
        add<double>(deltaJ->len, deltaJ->data, deltaJ->data, deltaJ2);

        // Average
        smul(deltaJ->len, deltaJ->data, 0.5, deltaJ->data);
        delete[] deltaJ2;
    }
    debug_printf(DP_DEBUG1, "Done computing DeltaJ\n");
}

/*
 * compute deltaJ, cost change associated with adding a sample 
 */
void computeDeltaJ2(const int D, double *deltaJ, const MDArray<double> *w, long *samples[], const long Nsamps[], const int w_type){
        
    debug_printf(DP_DEBUG1, "Computing first term in DeltaJ...\n");
    assert(D == 2);
    long nt = w->dims[D];
    long ksize = md_calc_size(D, w->dims);
    
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = w->data[t*w->strs[D] + t*w->strs[D+1]];
        md_fill(D, w->dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

    long csize = nt*ksize;
    for( long kt_new_ind = 0 ; kt_new_ind < csize ; kt_new_ind++ ){
        for( long t = 0 ; t < nt ; t++ ){
        for( long i = 0 ; i < Nsamps[t] ; i++ ){
            long k_new_sub[D+1];
            ind2sub<long>(D+1, w->dims, k_new_sub, kt_new_ind);

            // diff = (k_new - si, t_new, t)
            long diff[D+2]; 
            if( w_type == 1 ){
                //debug_printf(DP_DEBUG3, "[%ld %ld] - [%ld %ld] ", k_new_sub[0], k_new_sub[1], samples[t][D*i], samples[t][D*i+1]);
                sub<long>(D, diff, k_new_sub, &samples[t][D*i]);
            }else{
                //debug_printf(DP_DEBUG3, "[%ld %ld] + [%ld %ld] ", samples[t][D*i], samples[t][D*i+1], k_new_sub[0], k_new_sub[1] );
                add<long>(D, diff, k_new_sub, &samples[t][D*i]);
            }
            mod<long>(D, diff, diff, w->dims);
            diff[D+0] = k_new_sub[D];
            diff[D+1] = t;
            //debug_printf(DP_DEBUG3, " = [%ld %ld %ld %ld]", diff[0], diff[1], diff[D+0], diff[D+1]);
            //debug_printf(DP_DEBUG3, ";   ");
            //debug_printf(DP_DEBUG3, ", w(%ld,%ld) = %f ", diff[0], diff[1], w->data[sub2ind(D+2, w->strs, diff)]);
           
            // deltaJ[k_new] += w[diff]
            deltaJ[kt_new_ind] += w->data[sub2ind(D+2, w->strs, diff)];

            // diff = (si - k_new_sub, t, t_new);
            if( w_type == 1 ){
                //debug_printf(DP_DEBUG3, "[%ld %ld] - [%ld %ld] ", samples[t][D*i], samples[t][D*i+1], k_new_sub[0], k_new_sub[1] );
                sub<long>(D, diff, &samples[t][D*i], k_new_sub);
            }else{
                //debug_printf(DP_DEBUG3, "[%ld %ld] + [%ld %ld] ", samples[t][D*i], samples[t][D*i+1], k_new_sub[0], k_new_sub[1] );
                add<long>(D, diff, &samples[t][D*i], k_new_sub);
            }
            mod<long>(D, diff, diff, w->dims);
            diff[D+0] = t;
            diff[D+1] = k_new_sub[D];
            //debug_printf(DP_DEBUG3, " = [%ld %ld %ld %ld]", diff[0], diff[1], diff[D+0], diff[D+1]);

            // deltaJ[k_new] += w[diff]
            //debug_printf(DP_DEBUG3, ", w(%ld,%ld) = %f ", diff[0], diff[1], w->data[sub2ind(D+2, w->strs, diff)]);
            //debug_printf(DP_DEBUG3, "\n");
            deltaJ[kt_new_ind] += w->data[sub2ind(D+2, w->strs, diff)];

        } // end i loop
        } // end m loop
        debug_printf(DP_INFO, "\r%d%", 100.0 * (float) kt_new_ind / (float) csize);
    } // end k new loop
    debug_printf(DP_DEBUG1, "Done computing first term in DeltaJ.\n");
}

/* 
 * D         = dimension
 * p         = differential domain matrix of size [dims, Nbins, Nbins], dims is size vector of length D
 * mask_dims = sampling mask dims [dims Nbins]
 * samples   = gridded samples
 * N         = number of samples in each bin
*/
void computeDiffDist(const int D, MDArray<double> *p, 
        const long *samples[], 
        const long N[]){
    int i,j,bi,bj;

    //md_clear(p->D, p->dims, p->data, p->el_size);
    memset(p->data, 0, p->el_size*p->len);
    assert(p->dims[D+0] == p->dims[D+1]);

    for( bi=0 ; bi < p->dims[D+0] ; bi++ )
    for( bj=0 ; bj < p->dims[D+1] ; bj++ )
    for( i =0 ; i < N[bi] ; i++ )
    for( j =0 ; j < N[bj] ; j++ ){
        long diff[D+2]; /* si - sj + nbins */
        sub<long>(D, diff, &samples[bi][D*i], &samples[bj][D*j]);
        mod<long>(D, diff, diff, p->dims);
        diff[D+0] = bi;
        diff[D+1] = bj;
#ifdef DEBUG
        assert_in_bounds(D, diff, p->dims, "diff out of bounds of mask");
#endif
        p->data[sub2ind(D+2, p->strs, diff)]++;
    }
}

