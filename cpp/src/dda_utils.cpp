#include "misc.h"
#include "misc.hpp"
#include "multind.h"
#include "mdarray.h"
#include "debug.h"
#include "dda_utils.h"
#include "sampleHeap.h"
#include "config.h"

#include <algorithm>
#include <vector>
#include <random>
#include <limits>

//#define DEBUG

using std::vector;
            
const static double Inf = std::numeric_limits<double>::max();


void SparseKernel::Print() const {
    debug_printf(DP_INFO, "sparse w: \n");
    int nt = dims[kPhaseEncodeDims];
    for( int t2 = 0 ; t2 < nt ; t2++ )
    for( int t1 = 0 ; t1 < nt ; t1++ )
    for( int j = 0 ; j < len_w[t1][t2] ; j++ ){
        debug_printf(DP_INFO, "w(");
        for( int i = 0 ; i < kPhaseEncodeDims ; i++ ){
            debug_printf(DP_INFO, "%d, ", kw[t1][t2][j*D + i]);
        }
        debug_printf(DP_INFO,"%d, %d) = %f\n", t1, t2, w[t1][t2][j]);
    }
}

SparseKernel sparsifyWToK(const MDArray<4, double> &kernel, const long k){
    double T = kthLargestDouble(kernel.length(), kernel.data, k);
    return sparsifyW(kernel, T);
}

SparseKernel sparsifyW(const MDArray<4, double> &kernel, const double T){

    // dims
    long ksize = md_calc_size(kPhaseEncodeDims, kernel.dims);
    long nt = kernel.dims[kPhaseEncodeDims];
    SparseKernel sparse_kernel;

    debug_printf(DP_DEBUG3, "sparsifying w, nt: %d\n", nt);
    long strs[kPhaseEncodeDims+2];
    md_calc_strides(kPhaseEncodeDims+2, strs, kernel.dims, 1);
    memcpy(sparse_kernel.dims, kernel.dims, (kPhaseEncodeDims+2)*sizeof(long));


    // set sparse_kernel.w00tt
    for( long t = 0 ; t < nt ; t++ ){
        sparse_kernel.w00tt[t] = 0;
        sparse_kernel.w00tt[t] = sparse_kernel.w00tt[t] + kernel.data[t * strs[kPhaseEncodeDims] + t*strs[kPhaseEncodeDims+1]];
    }

    // Initialize len w
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        sparse_kernel.len_w[t1][t2] = 0;
    }
    // Allocate w k-v pairs
    //debug_printf(DP_INFO, "allocating..\n");
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        // TODO: don't need to allocate this
        sparse_kernel.w[t1][t2]  = (double *) xmalloc(ksize*sizeof(double));
        sparse_kernel.kw[t1][t2] = (long *) xmalloc(kPhaseEncodeDims*ksize*sizeof(long));
    }
    // Assign k-v pairs for w
    long supp = 0;
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        for( long k = 0 ; k < ksize ; k++ ){
            long ksub[kPhaseEncodeDims];
            ind2sub(kPhaseEncodeDims, kernel.dims, ksub, k);
            long kt1t2r_ind = k + t1*strs[kPhaseEncodeDims] + t2*strs[kPhaseEncodeDims+1];
            if( kernel.data[kt1t2r_ind] > T ){
                // append (k) to w_sp[t1][t2]
                /*
                debug_printf(DP_INFO, "Appending (%d %d %d %d), %f, len w: %d\n",
                            ksub[0], ksub[1], t1, t2, kernel.data[kt1t2r_ind],
                            sparse_kernel.len_w[t1][t2]);
                */
                supp++;
                sparse_kernel.w[t1][t2][sparse_kernel.len_w[t1][t2]] = kernel.data[kt1t2r_ind];
                memcpy(&sparse_kernel.kw[t1][t2][kPhaseEncodeDims*sparse_kernel.len_w[t1][t2]], ksub, kPhaseEncodeDims*sizeof(long));
                sparse_kernel.len_w[t1][t2]++;
            }
        }
    }
    debug_printf(DP_INFO, "w support: %d/%d = %f\n", supp, nt*nt*ksize, (float) supp / (float) (nt*nt*ksize));
    return sparse_kernel;
}

void printPat(const MDArray<3, int> &mask){
    static int maskPrintCount = 1;

    debug_printf(DP_INFO, "mask(:,%d) = [", maskPrintCount++);
    for( long i = 0 ; i < mask.length() ; i++ ){
        debug_printf(DP_INFO, "%d ", mask[i]);
    }
    debug_printf(DP_INFO, "];\n");
}

/*
 * Best candidate sampling but use a heap to store J
 */
void approxBestCandidate(const SparseKernel &wsp,
                         const vector<long> &max_samples_per_frame, 
                         const long totSamps, MDArray<3, int> &mask,
                         double &cost){
    //debug_level = DP_ALL;
    
    const long nt = wsp.dims[kPhaseEncodeDims];
    const long ksize = md_calc_size(kPhaseEncodeDims, wsp.dims);
    const long csize = nt*ksize;

    debug_printf(DP_INFO, "Approximate best candidate sampling, nt: %d\n", nt);

    long strs[kPhaseEncodeDims+2];
    md_calc_strides(kPhaseEncodeDims+2, strs, wsp.dims, 1);


    // Total cost
    cost = 0;
    MDArray<3, double> deltaJ(mask.dims);
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = wsp.w00tt[t];
        md_fill(kPhaseEncodeDims, wsp.dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

    // Clear out the mask
    mask.Clear();

    // Current number of samples
    vector<long> nSamplesCurrent(nt, 0);

    // Optional constraint on max samples in each pattern
    int n_full = 0;
    vector<bool> is_full(nt);
    for( long t = 0 ; t < nt ; t++ ){
        is_full[t] = max_samples_per_frame[t] == 0;
        if( is_full[t] ){
            n_full++;
        }
    }

    // Heapify delta J, later can use the fact that deltaJ is constant to speed this up
    debug_printf(DP_INFO, "Construct heap...\n");
    SampleHeap heap(deltaJ);

    // Best candidate selection loop
    debug_printf(DP_INFO, "Main loop...\n");
    for( long i = 0 ; i < totSamps ; i++ ){
        //heap.Print(wsp.dims);

        // TODO: pop or have an option if we can reacquire
        const Sample &snew = heap.getArr(0);
        
        vector<long> kt_new_sub(kPhaseEncodeDims+1);
        ind2sub(kPhaseEncodeDims+1, wsp.dims, &kt_new_sub[0], snew.getKTIndex());
        long t_new = kt_new_sub[kPhaseEncodeDims];

        debug_printf(DP_DEBUG4, "Sampling (%d %d)\n", 1+kt_new_sub[0], 1+kt_new_sub[1]);

        mask[snew.getKTIndex()]++;
        nSamplesCurrent[t_new]++;

        // --- Update total cost
        cost += deltaJ[snew.getKTIndex()];

        // --- Update insert cost (deltaJ), even for frames that are full
        for( long t = 0 ; t < nt ; t++ ){
            // w( Delta k, t, t' )
            for( long j = 0 ; j < wsp.len_w[t][t_new] ; j++ ){
                // tmp = k' + (Delta k)_j
                long tmp[kPhaseEncodeDims+1];
                add(kPhaseEncodeDims, tmp, &kt_new_sub[0], &wsp.kw[t][t_new][j*kPhaseEncodeDims]);
                mod(kPhaseEncodeDims, tmp, tmp, wsp.dims);
                tmp[kPhaseEncodeDims] = t;
                long i1 = sub2ind(kPhaseEncodeDims+1,strs,tmp);
                heap.increaseKey(heap.getKt2idx(i1), wsp.w[t][t_new][j]);
            }

            // w( Delta k, t', t)
            for( long j = 0 ; j < wsp.len_w[t_new][t] ; j++ ){
                // tmp = k' - (Delta k)_j
                long tmp[kPhaseEncodeDims+1];
                sub(kPhaseEncodeDims, tmp, &kt_new_sub[0], &wsp.kw[t_new][t][j*kPhaseEncodeDims]);
                mod(kPhaseEncodeDims, tmp, tmp, wsp.dims);
                tmp[kPhaseEncodeDims] = t;
                long i1 = sub2ind(kPhaseEncodeDims+1,strs,tmp);
                // deltaJ[k' + (Delta k)_j] += w[t'][t][j]
                heap.increaseKey(heap.getKt2idx(i1), wsp.w[t_new][t][j]);
            }
        } // end t loop
    
        //printPat(mask, D, wsp.dims);
        //heap.Print(wsp.dims);
        //debug_printf(DP_DEBUG4, "Done updating heap...\n");

        // --- Check if this pattern is full
        if( nSamplesCurrent[t_new] == max_samples_per_frame[t_new] ){
            is_full[t_new] = true;
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
            debug_printf(DP_DEBUG1, "\n%d%", (int) (100.0 * (float) i / (float) totSamps));
        }
    } // end i loop
    
    debug_printf(DP_INFO, "Done!\n");
}

// Functor to get the minimum element of an array breaking ties randomly. This does bias toward choosing elements
static long GetArgmin(const MDArray<3, double> &array) {
    assert(array.length() > 0);

    std::vector<long> perm(array.length());
    for( long i = 0 ; i < array.length() ; i++ )
        perm[i] = i;
    std::random_shuffle(perm.begin(), perm.end());

    long argmin = 0;
    double min = array[argmin];
    for( long j = 0 ; j < array.length() ; j++ ){
        //long i = perm[j];
        long i = j;
        if( array[i] < min || (array[i] == min && (rand() % 2 < 1)) ){
            argmin = i;
            min = array[i];
        }
    }
    return argmin;
}

/*
 * Best candidate sampling
 * compute deltaJ, cost change associated with adding a sample 
 *
 * INPUTS:
 *   w        = weighting cost function
 *   max_samples_per_frame = max samples per frame
 *   totSamps = total number of samples
 *
 * OUTPUTS:
 *   cost    = final cost
 *   deltaJ   = final insertion cost
 *   mask    = final sampling pattern
 *   samples = sample list
 */
void exactBestCandidate(const MDArray<4, double> &kernel, 
            const vector<long> &max_samples_per_frame, 
            const long totSamps, MDArray<3, int> &mask, double &cost,
            MDArray<3, double> &deltaJ){
    const long nt = kernel.dims[kPhaseEncodeDims];
    const long ksize = md_calc_size(kPhaseEncodeDims, kernel.dims);
    const long csize = nt*ksize;

    // Total cost
    cost = 0;
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = kernel.data[t*kernel.strs[kPhaseEncodeDims] + t*kernel.strs[kPhaseEncodeDims+1]];
        md_fill(kPhaseEncodeDims, kernel.dims, &deltaJ.data[t*ksize], &deltaJ_t, sizeof(double));
    }

    // Clear out the mask
    mask.Clear();

    // Current number of samples
    vector<long> nSamplesCurrent(nt, 0);

    // Optional constraint on max samples in each pattern
    int n_full = 0;
    vector<bool> is_full(nt);
    for( long t = 0 ; t < nt ; t++ ){
        is_full[t] = max_samples_per_frame[t] == 0;
        if( is_full[t] )
            n_full++;
    }

    // Best candidate selection loop
    for( long i = 0 ; i < totSamps ; i++ ){
        // k_new = argmin deltaJ(k,t) 
        long kt_new_ind = GetArgmin(deltaJ);
        long kt_new_sub[kPhaseEncodeDims+1];
        ind2sub(kPhaseEncodeDims+1, mask.dims, kt_new_sub, kt_new_ind);
        long t_new = kt_new_sub[kPhaseEncodeDims];
        myAssert(!is_full[t_new], "argmin should not have selected a full frame!");

        // --- Add sample
        // Increment mask
        mask[kt_new_ind]++;        
        nSamplesCurrent[t_new]++;

        // --- Update total cost
        cost += deltaJ[kt_new_ind];

        // --- Update insert cost (deltaJ), even for frames that are full
        for( long sj_ind = 0 ; sj_ind < csize ; sj_ind++ ){

            // diff = (kt_new_sub - sj, t_new, tj)
            long sj[kPhaseEncodeDims+1];
            ind2sub(kPhaseEncodeDims+1, mask.dims, sj, sj_ind);
            const long tj = sj[kPhaseEncodeDims];

            long diff[kPhaseEncodeDims+2];
            sub(kPhaseEncodeDims, diff, kt_new_sub, sj);
            mod(kPhaseEncodeDims, diff, diff, mask.dims);
            diff[kPhaseEncodeDims+0] = t_new;
            diff[kPhaseEncodeDims+1] = tj;
           
            // deltaJ[sj_ind] += w[snew - sj, t_new, tj]
            deltaJ.data[sj_ind] += kernel.data[sub2ind(kPhaseEncodeDims+2, kernel.strs, diff)];

            // diff = (sj - kt_new_sub, tj, t_new)
            sub(kPhaseEncodeDims, diff, sj, kt_new_sub);
            mod(kPhaseEncodeDims, diff, diff, mask.dims);
            diff[kPhaseEncodeDims+0] = tj;
            diff[kPhaseEncodeDims+1] = t_new;

            // deltaJ[sj_ind] += w[sj - snew, tj, t_new]
            deltaJ.data[sj_ind] += kernel.data[sub2ind(kPhaseEncodeDims+2, kernel.strs, diff)];
        } // end sj_ind loop

        // --- Check if this pattern is full
        if( nSamplesCurrent[t_new] == max_samples_per_frame[t_new] ){
            is_full[t_new] = true;
            n_full++;
            debug_printf(DP_DEBUG3, "%d is full\n", t_new);
            md_fill(kPhaseEncodeDims, kernel.dims, &deltaJ[t_new*ksize], &Inf, sizeof(double));
        }

        // --- Break if all patterns are full
        if( n_full == nt ){
            debug_printf(DP_DEBUG3, "All frames full, breaking at %d/%d\n", i, totSamps);
            break;
        }

        // --- Progress
        if( totSamps % (totSamps / 10)  == 0 ){
            debug_printf(DP_DEBUG1, "\r%d%", (int) (100.0 * (float) i / (float) totSamps));
        }
    } // end i loop
}

MDArray<3, double> computeDeltaJ(const MDArray<4, double> &kernel, const MDArray<3, int> &mask){

    MDArray<3, double> deltaJ(mask.dims);

    debug_printf(DP_DEBUG1, "Computing first term in DeltaJ...\n");
    const long nframes = kernel.dims[kPhaseEncodeDims];
    const long ksize = md_calc_size(kPhaseEncodeDims, kernel.dims);
    const long csize = nframes*ksize;
    
    for( long t = 0 ; t < nframes ; t++ ){
        double deltaJ_t = kernel.data[t*kernel.strs[kPhaseEncodeDims] + t*kernel.strs[kPhaseEncodeDims+1]];
        md_fill(kPhaseEncodeDims, kernel.dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

    for( long kt_new_ind = 0 ; kt_new_ind < csize ; kt_new_ind++ ){
        for( long kt_sample_ind = 0 ; kt_sample_ind < csize ; kt_sample_ind++ ){
            if( mask[kt_sample_ind] == 0 ){
                continue;
            }

            long kt_new_sub[kPhaseEncodeDims+1];
            ind2sub(kPhaseEncodeDims+1, mask.dims, kt_new_sub, kt_new_ind);

            long kt_sample_sub[kPhaseEncodeDims + 1];
            ind2sub(kPhaseEncodeDims+1, mask.dims, kt_sample_sub, kt_sample_ind);

            // diff = (k_new - si, t_new, t)
            long diff[kPhaseEncodeDims+2]; 
            sub(kPhaseEncodeDims, diff, kt_new_sub, kt_sample_sub);
            mod(kPhaseEncodeDims, diff, diff, mask.dims);
            diff[kPhaseEncodeDims+0] = kt_new_sub[kPhaseEncodeDims];
            diff[kPhaseEncodeDims+1] = kt_sample_sub[kPhaseEncodeDims];

            // deltaJ[k_new] += w[diff]
            deltaJ[kt_new_ind] += static_cast<double>(mask[kt_sample_ind]) * 
                kernel.data[sub2ind(kPhaseEncodeDims+2, kernel.strs, diff)];

            // diff = (si - kt_new_sub, t, t_new);
            sub(kPhaseEncodeDims, diff, kt_sample_sub, kt_new_sub);
            mod(kPhaseEncodeDims, diff, diff, mask.dims);
            diff[kPhaseEncodeDims+0] = kt_sample_sub[kPhaseEncodeDims];
            diff[kPhaseEncodeDims+1] = kt_new_sub[kPhaseEncodeDims];

            // deltaJ[k_new] += w[diff]
            deltaJ[kt_new_ind] += static_cast<double>(mask[kt_sample_ind]) * 
                kernel.data[sub2ind(kPhaseEncodeDims+2, kernel.strs, diff)];

        } // end i loop
        debug_printf(DP_DEBUG1, "\r%d%", 100.0 * (float) kt_new_ind / (float) csize);
    } // end k new loop
    debug_printf(DP_DEBUG1, "Done computing first term in DeltaJ.\n");

    return deltaJ;
}

/* 
 * Note this is slower than the method using the Fourier transform or the method that 
 * exploits the sparsity of the mask.
 */
MDArray<4, double> computeDiffDist(const MDArray<3, int> &mask){
    const long p_dims[4] = {mask.dims[0], mask.dims[1], mask.dims[2], mask.dims[2]};
    MDArray<4, double> p(p_dims);
    p.Clear();

    for( int i =0 ; i < mask.length() ; i++ )
    for( int j =0 ; j < mask.length() ; j++ ){
        if( mask[i] == 0 || mask[j] == 0 ){
            continue;
        }

        long diff[kPhaseEncodeDims+2];

        long si[kPhaseEncodeDims + 1u];
        long sj[kPhaseEncodeDims + 1u];

        ind2sub(kPhaseEncodeDims + 1u, mask.dims, si, i);
        ind2sub(kPhaseEncodeDims + 1u, mask.dims, sj, j);

        sub(kPhaseEncodeDims, diff, si, sj);
        mod(kPhaseEncodeDims, diff, diff, p.dims);
        diff[kPhaseEncodeDims+0] = si[kPhaseEncodeDims];
        diff[kPhaseEncodeDims+1] = sj[kPhaseEncodeDims];
#ifdef DEBUG
        assert_in_bounds(kPhaseEncodeDims, diff, p.dims, "diff out of bounds of mask");
#endif
        p.data[sub2ind(kPhaseEncodeDims+2, p.strs, diff)]++;
    }
    return p;
}

