#include "misc.h"
#include "misc.hpp"
#include "multind.h"
#include "mdarray.h"
#include "debug.h"
#include "dda_utils.h"
#include "sampleHeap.h"

#include <algorithm>
#include <vector>
#include <limits>

//#define DEBUG

using std::vector;
            
const static double Inf = std::numeric_limits<double>::max();


void SparseKernel::Print() const {
    debug_printf(DP_INFO, "sparse w: \n");
    int nt = dims[D];
    for( int t2 = 0 ; t2 < nt ; t2++ )
    for( int t1 = 0 ; t1 < nt ; t1++ )
    for( int j = 0 ; j < len_w[t1][t2] ; j++ ){
        debug_printf(DP_INFO, "w(");
        for( int i = 0 ; i < D ; i++ ){
            debug_printf(DP_INFO, "%d, ", kw[t1][t2][j*D + i]);
        }
        debug_printf(DP_INFO,"%d, %d) = %f\n", t1, t2, w[t1][t2][j]);
    }
}

SparseKernel sparsifyWToK(const MDArray<4, double> &wmd, const long k){
    double T = kthLargestDouble(wmd.len, wmd.data, k);
    return sparsifyW(wmd, T);
}

SparseKernel sparsifyW(const MDArray<4, double> &wmd, const double T){

    // dims
    long ksize = md_calc_size(kPhaseEncodeDims, wmd.dims);
    long nt = wmd.dims[kPhaseEncodeDims];
    SparseKernel wsp;

    debug_printf(DP_DEBUG3, "sparsifying w, nt: %d\n", nt);
    long strs[kPhaseEncodeDims+2];
    md_calc_strides(kPhaseEncodeDims+2, strs, wmd.dims, 1);
    memcpy(wsp.dims, wmd.dims, (kPhaseEncodeDims+2)*sizeof(long));


    // set wsp.w00tt
    for( long t = 0 ; t < nt ; t++ ){
        wsp.w00tt[t] = 0;
        wsp.w00tt[t] = wsp.w00tt[t] + wmd.data[t * strs[kPhaseEncodeDims] + t*strs[kPhaseEncodeDims+1]];
    }

    // Initialize len w
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        wsp.len_w[t1][t2] = 0;
    }
    // Allocate w k-v pairs
    //debug_printf(DP_INFO, "allocating..\n");
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        // TODO: don't need to allocate this
        wsp.w[t1][t2]  = (double *) xmalloc(ksize*sizeof(double));
        wsp.kw[t1][t2] = (long *) xmalloc(kPhaseEncodeDims*ksize*sizeof(long));
    }
    // Assign k-v pairs for w
    long supp = 0;
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        for( long k = 0 ; k < ksize ; k++ ){
            long ksub[kPhaseEncodeDims];
            ind2sub<long>(kPhaseEncodeDims, wmd.dims, ksub, k);
            long kt1t2r_ind = k + t1*strs[kPhaseEncodeDims] + t2*strs[kPhaseEncodeDims+1];
            if( wmd.data[kt1t2r_ind] > T ){
                // append (k) to w_sp[t1][t2]
                /*
                debug_printf(DP_INFO, "Appending (%d %d %d %d), %f, len w: %d\n",
                            ksub[0], ksub[1], t1, t2, wmd.data[kt1t2r_ind],
                            wsp.len_w[t1][t2]);
                */
                supp++;
                wsp.w[t1][t2][wsp.len_w[t1][t2]] = wmd.data[kt1t2r_ind];
                memcpy(&wsp.kw[t1][t2][kPhaseEncodeDims*wsp.len_w[t1][t2]], ksub, kPhaseEncodeDims*sizeof(long));
                wsp.len_w[t1][t2]++;
            }
        }
    }
    debug_printf(DP_INFO, "w support: %d/%d = %f\n", supp, nt*nt*ksize, (float) supp / (float) (nt*nt*ksize));
    return wsp;
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
void approxBestCandidate(const int D, double &cost, MDArray<3, double> &deltaJ, 
                                 int *mask, const SparseKernel &wsp,
                                 const long maxSamps[], 
                                 const long totSamps){
    //debug_level = DP_ALL;
    
    long nt = wsp.dims[D];

    debug_printf(DP_INFO, "Approximate best candidate sampling, nt: %d\n", nt);

    long strs[D+2];
    md_calc_strides(D+2, strs, wsp.dims, 1);
    long ksize = md_calc_size(D, wsp.dims);
    long csize = nt*ksize;

    // Total cost
    cost = 0;
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = wsp.w00tt[t];
        md_fill(D, wsp.dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
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
    SampleHeap heap(deltaJ);

    // Best candidate selection loop
    debug_printf(DP_INFO, "Main loop...\n");
    for( long i = 0 ; i < totSamps ; i++ ){
        heap.Print(wsp.dims);

        // TODO: pop or have an option if we can reacquire
        const Sample &snew = heap.getArr(0);
        
        long k_new_sub[D+1];
        ind2sub<long>(D+1, wsp.dims, k_new_sub, snew.getKTIndex());
        long t_new = k_new_sub[D];

        debug_printf(DP_DEBUG4, "Sampling (%d %d)\n", 1+k_new_sub[0], 1+k_new_sub[1]);

        mask[snew.getKTIndex()]++;
        nSampsCurrent[t_new]++;

        // --- Update total cost
        cost += deltaJ[snew.getKTIndex()];

        // --- Update insert cost (deltaJ), even for frames that are full
        for( long t = 0 ; t < nt ; t++ ){

            // w( Delta k, t, t' )
            for( long j = 0 ; j < wsp.len_w[t][t_new] ; j++ ){
                // tmp = k' + (Delta k)_j
                long tmp[D+1];
                add<long>(D, tmp, k_new_sub, &wsp.kw[t][t_new][j*D]);
                mod<long>(D, tmp, tmp, wsp.dims);
                tmp[D] = t;
                long i1 = sub2ind(D+1,strs,tmp);
                heap.increaseKey(heap.getKt2idx(i1), wsp.w[t][t_new][j]);
            }

            // w( Delta k, t', t)
             for( long j = 0 ; j < wsp.len_w[t_new][t] ; j++ ){
                // tmp = k' - (Delta k)_j
                long tmp[D+1];
                sub<long>(D, tmp, k_new_sub, &wsp.kw[t_new][t][j*D]);
                mod<long>(D, tmp, tmp, wsp.dims);
                tmp[D] = t;
                long i1 = sub2ind(D+1,strs,tmp);
                // deltaJ[k' + (Delta k)_j] += w[t'][t][j]
                heap.increaseKey(heap.getKt2idx(i1), wsp.w[t_new][t][j]);
            }
        } // end t loop
    
        //printPat(mask, D, wsp.dims);
        heap.Print(wsp.dims);
        //debug_printf(DP_DEBUG4, "Done updating heap...\n");

        // --- Check if this pattern is full
        if( nSampsCurrent[t_new] == maxSamps[t_new] ){
            is_full[t_new] = 1;
            n_full++;
            debug_printf(DP_DEBUG3, "%d is full\n", t_new);
            // TODO: do something about this
            //md_fill(D, wsp.dims, &deltaJ[t_new*ksize], &Inf, sizeof(double));
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
void exactBestCandidate(const int D, double &cost, MDArray<3, double> &deltaJ, 
            int *mask, const MDArray<4, double> &w, 
            long *samples[], const long maxSamps[], 
            const long totSamps){
    //assert(D == 2);
    const long nt = w.dims[D];
    long ksize = md_calc_size(D, w.dims);
    long csize = nt*ksize;

    // Total cost
    cost = 0;
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = w.data[t*w.strs[D] + t*w.strs[D+1]];
        md_fill(D, w.dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
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
        long kt_new_ind = argmin<double>(csize, deltaJ.data);
        long k_new_sub[D+1];
        ind2sub<long>(D+1, w.dims, k_new_sub, kt_new_ind);
        long t_new = k_new_sub[D];
        myAssert(!is_full[t_new], "argmin should not have selected a full frame!");

        // --- Add sample: TODO one function
        // Increment mask
        mask[kt_new_ind]++;
        // Insert into sample list
        memcpy(&samples[t_new][D*nSampsCurrent[t_new]], k_new_sub, D*sizeof(long));
        nSampsCurrent[t_new]++;

        // --- Update total cost
        cost += deltaJ[kt_new_ind];

        // --- Update insert cost (deltaJ), even for frames that are full
        for( long t = 0 ; t < nt ; t++ ){
        for( long j = t*ksize ; j < (t+1)*ksize ; j++ ){

            // diff = (k_new_sub - sj, t_new, t)
            long sj[D+1];
            ind2sub<long>(D+1, w.dims, sj, j);

            long diff[D+3]; 
            sub<long>(D, diff, k_new_sub, sj);
            mod<long>(D, diff, diff, w.dims);
            diff[D+0] = t_new;
            diff[D+1] = t;
           
            // deltaJ[j] += w[snew - sj, t_new, tj]
            deltaJ[j] += w.data[sub2ind(D+2, w.strs, diff)];

            // diff = (sj - k_new_sub, t, t_new)
            sub<long>(D, diff, sj, k_new_sub);
            mod<long>(D, diff, diff, w.dims);
            diff[D+0] = t;
            diff[D+1] = t_new;

            // deltaJ[j] += w[sj - snew, tj, t_new]
            deltaJ[j] += w.data[sub2ind(D+2, w.strs, diff)];
        } // end j loop
        } // end m loop

        // --- Check if this pattern is full
        if( nSampsCurrent[t_new] == maxSamps[t_new] ){
            is_full[t_new] = 1;
            n_full++;
            debug_printf(DP_DEBUG3, "%d is full\n", t_new);
            md_fill(D, w.dims, &deltaJ[t_new*ksize], &Inf, sizeof(double));
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

MDArray<3, double> computeDeltaJ(const MDArray<4, double> &w, const vector<vector<long> > &samples){

    // Create output
    const int kPhaseEncodeDims = 2u;

    MDArray<3, double> deltaJ(w.dims);

    debug_printf(DP_DEBUG1, "Computing first term in DeltaJ...\n");
    long nt = w.dims[kPhaseEncodeDims];
    long ksize = md_calc_size(kPhaseEncodeDims, w.dims);
    
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = w.data[t*w.strs[kPhaseEncodeDims] + t*w.strs[kPhaseEncodeDims+1]];
        md_fill(kPhaseEncodeDims, w.dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

    long csize = nt*ksize;
    for( long kt_new_ind = 0 ; kt_new_ind < csize ; kt_new_ind++ ){
        for( long t = 0 ; t < nt ; t++ ){
        for( long i = 0 ; i < samples[t].size() ; i++ ){
            long k_new_sub[kPhaseEncodeDims+1];
            ind2sub<long>(kPhaseEncodeDims+1, w.dims, k_new_sub, kt_new_ind);

            // diff = (k_new - si, t_new, t)
            long diff[kPhaseEncodeDims+2]; 
            //debug_printf(DP_DEBUG3, "[%ld %ld] - [%ld %ld] ", k_new_sub[0], k_new_sub[1], samples[t][D*i], samples[t][D*i+1]);
            sub<long>(kPhaseEncodeDims, diff, k_new_sub, &samples[t][kPhaseEncodeDims*i]);
            mod<long>(kPhaseEncodeDims, diff, diff, w.dims);
            diff[kPhaseEncodeDims+0] = k_new_sub[kPhaseEncodeDims];
            diff[kPhaseEncodeDims+1] = t;
            //debug_printf(DP_DEBUG3, " = [%ld %ld %ld %ld]", diff[0], diff[1], diff[kPhaseEncodeDims+0], diff[kPhaseEncodeDims+1]);
            //debug_printf(DP_DEBUG3, ";   ");
            //debug_printf(DP_DEBUG3, ", w(%ld,%ld) = %f ", diff[0], diff[1], w.data[sub2ind(D+2, w.strs, diff)]);
           
            // deltaJ[k_new] += w[diff]
            deltaJ[kt_new_ind] += w.data[sub2ind(kPhaseEncodeDims+2, w.strs, diff)];

            // diff = (si - k_new_sub, t, t_new);
            //debug_printf(DP_DEBUG3, "[%ld %ld] - [%ld %ld] ", samples[t][kPhaseEncodeDims*i], samples[t][kPhaseEncodeDims*i+1], k_new_sub[0], k_new_sub[1] );
            sub<long>(kPhaseEncodeDims, diff, &samples[t][kPhaseEncodeDims*i], k_new_sub);
            mod<long>(kPhaseEncodeDims, diff, diff, w.dims);
            diff[kPhaseEncodeDims+0] = t;
            diff[kPhaseEncodeDims+1] = k_new_sub[kPhaseEncodeDims];
            //debug_printf(DP_DEBUG3, " = [%ld %ld %ld %ld]", diff[0], diff[1], diff[kPhaseEncodeDims+0], diff[kPhaseEncodeDims+1]);

            // deltaJ[k_new] += w[diff]
            //debug_printf(DP_DEBUG3, ", w(%ld,%ld) = %f ", diff[0], diff[1], w.data[sub2ind(D+2, w.strs, diff)]);
            //debug_printf(DP_DEBUG3, "\n");
            deltaJ[kt_new_ind] += w.data[sub2ind(kPhaseEncodeDims+2, w.strs, diff)];

        } // end i loop
        } // end m loop
        debug_printf(DP_INFO, "\r%d%", 100.0 * (float) kt_new_ind / (float) csize);
    } // end k new loop
    debug_printf(DP_DEBUG1, "Done computing first term in DeltaJ.\n");

    return deltaJ;
}

/* 
 * D         = dimension
 * p         = differential domain matrix of size [dims, Nbins, Nbins], dims is size vector of length D
 * mask_dims = sampling mask dims [dims Nbins]
 * samples   = gridded samples
 * N         = number of samples in each bin
*/
void computeDiffDist(const int D, 
        const long *samples[], 
        const long N[],
        MDArray<4, double> &p){
    int i,j,bi,bj;

    //md_clear(p.D, p.dims, p.data, p.el_size);
    p.Clear();
    assert(p.dims[D+0] == p.dims[D+1]);

    for( bi=0 ; bi < p.dims[D+0] ; bi++ )
    for( bj=0 ; bj < p.dims[D+1] ; bj++ )
    for( i =0 ; i < N[bi] ; i++ )
    for( j =0 ; j < N[bj] ; j++ ){
        long diff[D+2]; /* si - sj + nbins */
        sub<long>(D, diff, &samples[bi][D*i], &samples[bj][D*j]);
        mod<long>(D, diff, diff, p.dims);
        diff[D+0] = bi;
        diff[D+1] = bj;
#ifdef DEBUG
        assert_in_bounds(D, diff, p.dims, "diff out of bounds of mask");
#endif
        p.data[sub2ind(D+2, p.strs, diff)]++;
    }
}

