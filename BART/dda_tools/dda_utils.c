#define _GNU_SOURCE
#include <math.h>
#include <complex.h>
#include <float.h>
#include "misc/debug.h"
#include "misc/misc.h"
#include "misc/mri.h"


#include "misc_dd.h"

#include "num/fft.h"
#include "num/multind.h"
#include "num/flpmath.h"

#include "dda_utils.h"
#include "sampleHeap.h"

const double Inf = DBL_MAX;

void print_wsp(const SparseW *wsp){
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

void sparsify_w_to_k(const long D, const long dims[], SparseW *wsp, double *wmd, const long k){

    double T = kthLargest(md_calc_size(D, dims) * dims[D] * dims[D], wmd, k);
    sparsify_w(D, dims, wsp, wmd, T);
}

void sparsify_w(const long D, const long dims[], SparseW *wsp, double *wmd, const double T){

    debug_printf(DP_DEBUG3, "sparsifying w\n");
    // dims
    long ksize = md_calc_size(D, dims);
    long nt = dims[D];
    long nr = dims[D+2];
    long strs1[D+3];
    md_calc_strides(D+3, strs1, dims, 1); // TODO D+3???
    wsp->D = D;
    myAssert(D < DIMS-2, "D < max_dims");
    memcpy(wsp->dims, dims, DIMS*sizeof(long));

    // set wsp->w00tt
    for( long t = 0 ; t < nt ; t++ ){
        wsp->w00tt[t] = 0;
        for( long r = 0 ; r < nr ; r++ )
            wsp->w00tt[t] = wsp->w00tt[t] + wmd[t * strs1[D] + t*strs1[D+1] + r*strs1[D+2]];
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
        wsp->w[t1][t2][r]  = (double *) xmalloc(ksize*sizeof(double));
        wsp->kw[t1][t2][r] = (long *) xmalloc(D*ksize*sizeof(long));
    }
    // Assign k-v pairs for w
    //debug_printf(DP_INFO, "assigning k-v pairs...\n");
    long supp = 0;

    for( long r = 0 ; r < nr ; r++ )
    for( long t2 = 0 ; t2 < nt ; t2++ )
    for( long t1 = 0 ; t1 < nt ; t1++ ){
        for( long k = 0 ; k < ksize ; k++ ){
            long ksub[D];
            idx2sub(D, dims, ksub, k);
            long kt1t2r_ind = k + t1*strs1[D] + t2*strs1[D+1] + r*strs1[D+2];
            if( wmd[kt1t2r_ind] > T ){
                // append (k) to w_sp[t1][t2]
                /*
                debug_printf(DP_INFO, "Appending (%d %d %d %d), %f, len w: %d\n",
                            ksub[0], ksub[1], t1, t2, wmd[kt1t2r_ind],
                            wsp->len_w[t1][t2]);
                */
                supp++;
                wsp->w[t1][t2][r][wsp->len_w[t1][t2][r]] = wmd[kt1t2r_ind];
                memcpy(&wsp->kw[t1][t2][r][D*wsp->len_w[t1][t2][r]], ksub, D*sizeof(long));
                wsp->len_w[t1][t2][r]++;
            }
        }
    }
    debug_printf(DP_INFO, "w support: %d/%d = %f\n", supp, nt*nt*ksize, (float) supp / (float) (nt*nt*ksize));
}

/*
static void printPat(const int *mask, const long D, const long dims[]){
    static int maskPrintCount = 1;

    debug_printf(DP_INFO, "mask(:,%d) = [", maskPrintCount++);
    long N  = md_calc_size(D+1, dims);
    for( long i = 0 ; i < N ; i++ ){
        debug_printf(DP_INFO, "%d ", mask[i]);
    }
    debug_printf(DP_INFO, "];\n");
}
*/

/*
 * Best candidate sampling but use a heap to store J
 */
void approxBestCandidate(const int D,
                         double *cost,
                         double *deltaJ,
                         int *mask,
                         const SparseW *wsp,
                         const long maxSamps[],
                         const long totSamps)
{
    debug_printf(DP_DEBUG1, "Approximate best candidate smapling...\n");
    
    long nt = wsp->dims[D];
    long nr = wsp->dims[D+2];
    long strs1[D+2];
    md_calc_strides(D+2, strs1, wsp->dims, 1);
    long ksize = md_calc_size(D, wsp->dims);
    long csize = nt*ksize;

    // Total cost
    *cost = 0;
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = wsp->w00tt[t];
        md_fill(D, wsp->dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

    // Clear out the mask
    memset(mask, 0, csize*sizeof(int));

    // Current number of samples
    long nSampsCurrent[nt];
    memset(nSampsCurrent, 0, nt*sizeof(long));

    // Optional constraint on max samples in each pattern
    int n_full = 0;

    // Heapify delta J, later can use the fact that deltaJ is constant to speed this up
    debug_printf(DP_INFO, "Create heap...\n");
    SampleHeap *heap = build_sample_heap(csize, deltaJ);
        
    // Best candidate selection loop
    debug_printf(DP_INFO, "Main loop...\n");
    for( long i = 0 ; i < totSamps ; i++ ){
        /*
        printHeap(heap, wsp->dims);
        if( !validateHeap(heap) ){
            debug_printf(DP_ERROR, "Invalid heap!\n");
            printHeap(heap, wsp->dims);
            exit(0);
        }
        */

        Sample snew;
        long k_new_sub[D+1];
        long t_new;
        do{
            //snew = heap->arr[0];
            snew = pop(heap);
            idx2sub(D+1, wsp->dims, k_new_sub, snew.kt);
            t_new = k_new_sub[D];
        }while( nSampsCurrent[t_new] >= maxSamps[t_new] );

        debug_printf(DP_DEBUG4, "Sampling (%d %d %d)\n", 1+k_new_sub[0], 1+k_new_sub[1], 1+t_new);

        mask[snew.kt]++;
        nSampsCurrent[t_new]++;

        // --- Update total cost
        *cost += deltaJ[snew.kt];

        // --- Update insert cost (deltaJ), even for frames that are full
        for( long r = 0 ; r < nr ; r++ )
        for( long t = 0 ; t < nt ; t++ ){

            // w( Delta k, t, t' )
            for( long j = 0 ; j < wsp->len_w[t][t_new][r] ; j++ ){
                // tmp = k' + (Delta k)_j
                long tmp[D+1];
                addl(D, tmp, k_new_sub, &wsp->kw[t][t_new][r][j*D]);
                modl(D, tmp, tmp, wsp->dims);
                tmp[D] = t;
                long i1 = sub2idx(D+1,strs1,tmp);
                debug_printf(DP_DEBUG4, "Increase DeltaJ(%d,%d) by %f\n", tmp[0]+1, tmp[1]+1, wsp->w[t][t_new][r][j]);
                increaseKey(heap, heap->kt2idx[i1], wsp->w[t][t_new][r][j]);
            }

            // w( Delta k, t', t)
             for( long j = 0 ; j < wsp->len_w[t_new][t][r] ; j++ ){
                // tmp = k' - (Delta k)_j
                long tmp[D+1];
                subl(D, tmp, k_new_sub, &wsp->kw[t_new][t][r][j*D]);
                modl(D, tmp, tmp, wsp->dims);
                tmp[D] = t;
                long i1 = sub2idx(D+1,strs1,tmp);
                // deltaJ[k' + (Delta k)_j] += w[t'][t][j]
                debug_printf(DP_DEBUG4, "Increase DeltaJ(%d,%d) by %f\n", tmp[0]+1, tmp[1]+1, wsp->w[t_new][t][r][j]);
                increaseKey(heap, heap->kt2idx[i1], wsp->w[t_new][t][r][j]);
            }
        } // end t loop
    
        //printPat(mask, D, wsp->dims);
        //printHeap(heap, wsp->dims);
        //debug_printf(DP_DEBUG4, "Done updating heap...\n");

        // --- Check if this pattern is full
        if( nSampsCurrent[t_new] == maxSamps[t_new] ){
            n_full++;
            debug_printf(DP_DEBUG1, "%d is full\n", t_new);
        }

        // --- Break if all patterns are full
        if( n_full == nt ){
            debug_printf(DP_DEBUG1, "All frames full, breaking at %d/%d\n", i+1, totSamps);
            break;
        }

        // --- Progress
        if( totSamps > 10 && totSamps % (totSamps / 10)  == 0 ){
            debug_printf(DP_INFO, "\r%d%", (int) (100.0 * (float) i / (float) totSamps));
        }
    } // end i loop
    
    free_heap(heap);
    debug_printf(DP_INFO, "Approximate BC done!\n");
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
void exactBestCandidate(const int D, const long dims[], double *cost, double *deltaJ,
                        int *mask, const double *w,
                        long *samples[], const long maxSamps[],
                        const long totSamps)
{
    //assert(D == 2);
    const long REIM_DIM = 4;
    long nt = dims[D];
    long ksize = md_calc_size(D, dims);
    long csize = nt*ksize;

    long strs1[D];
    md_calc_strides(D+2, strs1, dims, 1);

    // Total cost
    *cost = 0;
    
    // Initialize deltaJ
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = w[t*strs1[D] + t*strs1[D+1]];
        md_fill(D, dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
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
        long kt_new_ind = argmind(csize, deltaJ);
        long k_new_sub[D+1];
        idx2sub(D+1, dims, k_new_sub, kt_new_ind);
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
            idx2sub(D+1, dims, sj, j);

            long diff[D+3]; 
            subl(D, diff, k_new_sub, sj);
            modl(D, diff, diff, dims);
            diff[D+0] = t_new;
            diff[D+1] = t;
           
            // deltaJ[j] += w[snew - sj, t_new, tj]
            deltaJ[j] += w[sub2idx(D+2, strs1, diff)];

            // diff = (sj - k_new_sub, t, t_new)
            subl(D, diff, sj, k_new_sub);
            modl(D, diff, diff, dims);
            diff[D+0] = t;
            diff[D+1] = t_new;

            // deltaJ[j] += w[sj - snew, tj, t_new]
            deltaJ[j] += w[sub2idx(D+2, strs1, diff)];
    
            // Part 2
            if( dims[REIM_DIM] == 2){
                // diff = (sj + k_new_sub, t, t_new)
                addl(D, diff, sj, k_new_sub);
                modl(D, diff, diff, dims);
                diff[D+0] = t;
                diff[D+1] = t_new;
                diff[D+2] = 1;

                // deltaJ[j] += w[sj + snew, tj, t_new, 1]
                long i1 = sub2idx(D+3, strs1, diff);
                deltaJ[j] += w[i1];
                // deltaJ[j] += w[sj + snew, t_new, tj, 1]
                diff[D+0] = t_new;
                diff[D+1] = t;
                i1 = sub2idx(D+3, strs1, diff);
                deltaJ[j] += w[i1];
            }
        } // end j loop
        } // end m loop

        // --- Check if this pattern is full
        if( nSampsCurrent[t_new] == maxSamps[t_new] ){
            is_full[t_new] = 1;
            n_full++;
            debug_printf(DP_DEBUG3, "%d is full\n", t_new);
            md_fill(D, dims, &deltaJ[t_new*ksize], &Inf, sizeof(double));
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

/*
 * compute deltaJ, cost change associated with adding a sample 
 */
static void compute_deltaj2(const int D, const long dims[], double *deltaJ, const double *w, long *samples[], const long Nsamps[], const int w_type){
        
    debug_printf(DP_DEBUG1, "Computing first term in DeltaJ...\n");
    assert(D == 2);
    long nt = dims[D];
    long ksize = md_calc_size(D, dims);

    long strs1[D+2];
    md_calc_strides(D+2, strs1, dims, 1);
    
    for( long t = 0 ; t < nt ; t++ ){
        double deltaJ_t = w[t*strs1[D] + t*strs1[D+1]];
        md_fill(D, dims, &deltaJ[t*ksize], &deltaJ_t, sizeof(double));
    }

    long csize = nt*ksize;
    for( long kt_new_ind = 0 ; kt_new_ind < csize ; kt_new_ind++ ){
        for( long t = 0 ; t < nt ; t++ ){
        for( long i = 0 ; i < Nsamps[t] ; i++ ){
            long k_new_sub[D+1];
            idx2sub(D+1, dims, k_new_sub, kt_new_ind);

            // diff = (k_new - si, t_new, t)
            long diff[D+2]; 
            if( w_type == 1 ){
                //debug_printf(DP_DEBUG3, "[%ld %ld] - [%ld %ld] ", k_new_sub[0], k_new_sub[1], samples[t][D*i], samples[t][D*i+1]);
                subl(D, diff, k_new_sub, &samples[t][D*i]);
            }else{
                //debug_printf(DP_DEBUG3, "[%ld %ld] + [%ld %ld] ", samples[t][D*i], samples[t][D*i+1], k_new_sub[0], k_new_sub[1] );
                addl(D, diff, k_new_sub, &samples[t][D*i]);
            }
            modl(D, diff, diff, dims);
            diff[D+0] = k_new_sub[D];
            diff[D+1] = t;
            //debug_printf(DP_DEBUG3, " = [%ld %ld %ld %ld]", diff[0], diff[1], diff[D+0], diff[D+1]);
            //debug_printf(DP_DEBUG3, ";   ");
            //debug_printf(DP_DEBUG3, ", w(%ld,%ld) = %f ", diff[0], diff[1], w[sub2idx(D+2, w->strs1, diff)]);
           
            // deltaJ[k_new] += w[diff]
            deltaJ[kt_new_ind] += w[sub2idx(D+2, strs1, diff)];

            // diff = (si - k_new_sub, t, t_new);
            if( w_type == 1 ){
                //debug_printf(DP_DEBUG3, "[%ld %ld] - [%ld %ld] ", samples[t][D*i], samples[t][D*i+1], k_new_sub[0], k_new_sub[1] );
                subl(D, diff, &samples[t][D*i], k_new_sub);
            }else{
                //debug_printf(DP_DEBUG3, "[%ld %ld] + [%ld %ld] ", samples[t][D*i], samples[t][D*i+1], k_new_sub[0], k_new_sub[1] );
                addl(D, diff, &samples[t][D*i], k_new_sub);
            }
            modl(D, diff, diff, dims);
            diff[D+0] = t;
            diff[D+1] = k_new_sub[D];
            //debug_printf(DP_DEBUG3, " = [%ld %ld %ld %ld]", diff[0], diff[1], diff[D+0], diff[D+1]);

            // deltaJ[k_new] += w[diff]
            //debug_printf(DP_DEBUG3, ", w(%ld,%ld) = %f ", diff[0], diff[1], w[sub2idx(D+2, w->strs1, diff)]);
            //debug_printf(DP_DEBUG3, "\n");
            deltaJ[kt_new_ind] += w[sub2idx(D+2, strs1, diff)];

        } // end i loop
        } // end m loop
        if( 0 == kt_new_ind % (csize/100) ){
            debug_printf(DP_INFO, "\r%.0f%", 100.0 * (float) kt_new_ind / (float) csize);
        }
    } // end k new loop
    debug_printf(DP_DEBUG1, "Done computing first term in DeltaJ.\n");
}

void computeDeltaJ(const int D, const long dims[], double *deltaJ, const double *w, long *samples[], const long Nsamps[]){
    debug_printf(DP_DEBUG1, "Computing DeltaJ...\n");

    compute_deltaj2(D, dims, deltaJ, w, samples, Nsamps, 1);
    long deltaJ_size = md_calc_size(D+1, dims);

    int Nw = dims[D+2];
    if( Nw == 2){
        debug_printf(DP_DEBUG1, "Computing second term in DeltaJ for real-valued images...\n");
        // Compute the second component
        double *deltaJ2 = xmalloc(deltaJ_size * sizeof(double));
        compute_deltaj2(D, dims, deltaJ2, w, samples, Nsamps, 2);

        // Add to the cost
        addd(deltaJ_size, deltaJ, deltaJ, deltaJ2);

        // Average
        smuld(deltaJ_size, deltaJ, 0.5, deltaJ);
        free(deltaJ2);
    }
    debug_printf(DP_DEBUG1, "Done computing DeltaJ\n");
}

void circRev2d(const long dims[], complex float *dst, const complex float *src){
    dst[0] = src[0];
    for( long i = 1 ; i < dims[0] ; ++i ){
        dst[i] = src[dims[0] - i];
    }
    for( long i = 1 ; i < dims[1] ; ++i ){
        dst[dims[0]*i] = src[dims[0]*(dims[1] - i)];
    }
    for( long kz = 1 ; kz < dims[1] ; kz++ )
    for( long ky = 1 ; ky < dims[0] ; ky++ ){
        dst[ky + kz*dims[0]] = src[(dims[0] - ky) + dims[0] * (dims[1] - kz)];
    }
}

void computeDeltaJFftMethod(const int kDims, const long w_dims[], complex float *deltaJ, const complex float *w, const complex float *pat){
    debug_printf(DP_DEBUG1, "Computing DeltaJ using the FFT method...\n");

    assert(w_dims[kDims] == w_dims[kDims+1]);

    long pat_dims[DIMS];
    md_select_dims(DIMS, 7, pat_dims, w_dims);

    long kSize = md_calc_size(kDims, w_dims);
    //long ktSize = md_calc_size(kDims+1, w_dims);
    long nt = w_dims[kDims];
    long wstrs1[kDims+2];
    md_calc_strides(kDims+2, wstrs1, w_dims, 1);

    complex float *Fmask = md_alloc(DIMS, pat_dims, CFL_SIZE);
    complex float *Fw    = md_alloc(DIMS, w_dims, CFL_SIZE);
    long tmp_dims[DIMS];
    md_select_dims(DIMS, 3u, tmp_dims, w_dims);
    complex float *tmp   = md_alloc(DIMS, tmp_dims, CFL_SIZE);

    fft(DIMS, pat_dims, 3u, Fmask, pat);
    fft(DIMS, w_dims, 3u, Fw, w);

    for( int t1 = 0 ; t1 < nt ; t1++ ){
        const complex float *w00tt    = &w[t1*wstrs1[kDims]+t1*wstrs1[kDims+1]];
        complex float *deltaJt = &deltaJ[t1*kSize];
        for( long i = 0 ; i < kSize ; i++ ){
            deltaJt[i] = w00tt[0];
        }
        for( int t2 = 0 ; t2 < nt ; t2++ ){
            //deltaJ(:,:,t1) = deltaJ(:,:,t1) + 
            //    ifft(ifft(Fmask(:,:,t2) .* (Fw(:,:,t1,t2) + circRev2d(Fw(:,:,t2,t1))), [], 1), [], 2);
            

            // tmp = circRev2d(Fw(:,:,t2,t1))
            circRev2d(tmp_dims, tmp, &Fw[t2*wstrs1[kDims]+t1*wstrs1[kDims+1]]);
            // tmp = tmp + Fw(:,:,t1,t2)
            md_zadd(DIMS, tmp_dims, tmp, tmp, &Fw[t1*wstrs1[kDims]+t2*wstrs1[kDims+1]]);

            // tmp = tmp * Fmask(:,:,t2)
            md_zmul(DIMS, tmp_dims, tmp, tmp, &Fmask[t2*kSize]);
            
            // tmp = ifft(tmp)
            ifft(DIMS, tmp_dims, 3u, tmp, tmp); // this line is broken
            md_zsmul(DIMS, tmp_dims, tmp, tmp, (complex float) 1.0 / (float) md_calc_size(2, tmp_dims));
            // deltaJ(:,:,t1) = deltaJ(:,:,t1) + tmp
            md_zadd(DIMS, tmp_dims, deltaJt, deltaJt, tmp);
        }
    }
    
    md_free(tmp);
    md_free(Fw);
    md_free(Fmask);
    debug_printf(DP_DEBUG1, "Done computing DeltaJ\n");
}

// Check sns_dims
void check_sns_dims(const long sns_dims[]){
    long sns_dims1[DIMS];
    md_select_dims(DIMS, ~(FFT_FLAGS | TE_FLAG | COIL_FLAG | MAPS_FLAG | MAPS2_FLAG), sns_dims1, sns_dims);

    for( unsigned int i = 0 ; i < DIMS ; i++ ){
        if( sns_dims1[i] != 1 ){
            debug_printf(DP_ERROR, "Bad dimensions for coil sensitivities: ");
            debug_print_dims(DP_ERROR, DIMS, sns_dims);
            debug_print_dims(DP_ERROR, DIMS, sns_dims1);
            assert(0);
        }
    }

}
void get_w_dims(long w_dims[], const long sns_dims[]){
    md_singleton_dims(DIMS, w_dims);
    w_dims[0] = sns_dims[PHS1_DIM];
    w_dims[1] = sns_dims[PHS2_DIM];
    w_dims[2] = w_dims[3] = sns_dims[TE_DIM];
}

static void md_zfloat2doublereal(const long D, const long dims[], double *dst, const complex float *src){
    long N = md_calc_size(D, dims);
    for( long i = 0 ; i < N ; i++ ){
        dst[i] = (double) crealf(src[i]);
    }
}
/*
*/

void buildW(double *w, const long sns_dims[], const complex float *sns_maps){
    // Set up complex float output
    long w_dims[DIMS];
    get_w_dims(w_dims, sns_dims);
    complex float *wc = md_calloc(DIMS, w_dims, CFL_SIZE);

    // buildW
    buildW2(wc, sns_dims, sns_maps);
 
    // Convert to double for output
    md_zfloat2doublereal(DIMS, w_dims, w, wc);
    md_free(wc);
}

/*
 *
 */
void buildW2(complex float *wc, const long sns_dims[], 
            const complex float *sns_maps){

    // Timing
    double ts[4];
    double start_time;
    start_time = timestamp();

    // permute sns_maps to [y z map map2 x coil time] for speed
    unsigned int order[DIMS] = {1,2,4,6,0,3,5,7,8,9,10,11,12,13,14,15};
    //unsigned int order[DIMS] = {1,2,4,5,0,3,10,6,7,8,9,11,12,13,14,15};

    long sp_dims[DIMS];
    md_permute_dims(DIMS, order, sp_dims, sns_dims); // sp[order[i]] = s[i]
    complex float *sns_maps_perm = md_alloc(DIMS, sp_dims, CFL_SIZE);
    md_permute(DIMS, order, sp_dims, sns_maps_perm, sns_dims, sns_maps, CFL_SIZE);
    long sp_strs1[DIMS];
    md_calc_strides(DIMS, sp_strs1, sp_dims, 1);

    // Set up output
    long w_dims[DIMS];
    get_w_dims(w_dims, sns_dims);

    long w_strs1[DIMS];
    md_calc_strides(DIMS, w_strs1, w_dims, 1);

    long NyNz = w_dims[0]*w_dims[1];
    md_clear(DIMS, w_dims, wc, CFL_SIZE);

    // Allocate temporary array to store the sensitivty maps in the product
    long tmp_dims[DIMS];
    md_select_dims(DIMS, ~(READ_FLAG | COIL_FLAG | TE_FLAG), tmp_dims, sns_dims);
    long tmp_strs1[DIMS];
    md_calc_strides(DIMS, tmp_strs1, tmp_dims, 1);
    long tmp_strs[DIMS];
    md_calc_strides(DIMS, tmp_strs, tmp_dims, CFL_SIZE);

    complex float *tmp = md_alloc(DIMS, tmp_dims, CFL_SIZE);

    // Phase encode dimensions
    long pe_dims[DIMS];
    md_select_dims(DIMS, PHS1_FLAG | PHS2_FLAG, pe_dims, sns_dims);
    long pe_strs[DIMS];
    md_calc_strides(DIMS, pe_strs, pe_dims, CFL_SIZE);

    // [ny nz ncoils ncoils nt nt] array of all products of sensitivities
    long sns_prod_dims[DIMS];
    md_singleton_dims(DIMS, sns_prod_dims);
    sns_prod_dims[0] = pe_dims[PHS1_DIM];
    sns_prod_dims[1] = pe_dims[PHS2_DIM];
    sns_prod_dims[3] = sns_dims[READ_DIM];
    sns_prod_dims[4] = sns_prod_dims[5] = sns_dims[COIL_DIM];
    sns_prod_dims[6] = sns_prod_dims[7] = sns_dims[TE_DIM];

    long sns_prod_strs1[DIMS];
    md_calc_strides(DIMS, sns_prod_strs1, sns_prod_dims, 1);

    // All products of sensitivities
    complex float *sns_prod = md_calloc(DIMS, sns_prod_dims, CFL_SIZE);
        
    // Print progress
    long n_loop = sns_dims[TE_DIM]*sns_dims[TE_DIM]
                 *sns_dims[COIL_DIM]*sns_dims[COIL_DIM]
                 *sns_dims[READ_DIM];

    long loop_count = 0;
    
    ts[0] = timestamp() - start_time;
    start_time = timestamp();

    for( long x = 0  ; x < sns_dims[READ_DIM]  ; x++ )
    for( long t2 = 0 ; t2 < sns_dims[TE_DIM] ; t2++ )
    for( long c2 = 0 ; c2 < sns_dims[COIL_DIM] ; c2++ ){
        // sns2 = sns_maps(x,:,:,c2,:,t2,:)
        for( long t1 = 0 ; t1 < sns_dims[TE_DIM] ; t1++ )
        for( long c1 = 0 ; c1 < sns_dims[COIL_DIM] ; c1++ ){
            // Product of (x,c1,t1) and (x,c2,t2)
            complex float *sns_prod1 = &sns_prod[ x * sns_prod_strs1[3] 
                                               + c1 * sns_prod_strs1[4] 
                                               + c2 * sns_prod_strs1[5] 
                                               + t1 * sns_prod_strs1[6] 
                                               + t2 * sns_prod_strs1[7]];
            // tmp = sns_maps(x,:,:,c1,:,t1,:)
            // tmp = tmp * conj(sns2)
            md_zmulc(DIMS, tmp_dims, tmp, 
                &sns_maps_perm[x*sp_strs1[4] + c1*sp_strs1[5] + t1*sp_strs1[6]], 
                &sns_maps_perm[x*sp_strs1[4] + c2*sp_strs1[5] + t2*sp_strs1[6]]);

            // sum over MAPS and MAPS2_DIMS
            // This part is the bottleneck (3 seconds)
            for( long l = 0 ; l < sns_dims[MAPS2_DIM] ; l++ )
            for( long m = 0 ; m < sns_dims[MAPS_DIM]  ; m++ )
            for( long r = 0 ; r < NyNz ; r++ )
                sns_prod1[r] += tmp[r + l*tmp_strs1[MAPS2_DIM] + m*tmp_strs1[MAPS_DIM]];
        }
    }
    md_free(tmp);

    ts[1] = timestamp() - start_time;
    start_time = timestamp();

    // sns_prod = F(sns_prod)
    fft(DIMS, sns_prod_dims, 3u, sns_prod, sns_prod);

    ts[2] = timestamp() - start_time;
    start_time = timestamp();

    // Compute summation
    for( long x = 0  ; x < sns_dims[READ_DIM]  ; x++ )
    for( long t1 = 0 ; t1 < sns_dims[TE_DIM] ; t1++ )
    for( long t2 = 0 ; t2 < sns_dims[TE_DIM] ; t2++ )
    for( long c1 = 0 ; c1 < sns_dims[COIL_DIM] ; c1++ )
    for( long c2 = 0 ; c2 < sns_dims[COIL_DIM] ; c2++ ){
        complex float *sns_prod1 = &sns_prod[ x * sns_prod_strs1[3] 
                                           + c1 * sns_prod_strs1[4] 
                                           + c2 * sns_prod_strs1[5] 
                                           + t1 * sns_prod_strs1[6] 
                                           + t2 * sns_prod_strs1[7]];

        // sns_prod = abs(sns_prod).^2
        md_zmulc(DIMS, pe_dims, sns_prod1, sns_prod1, sns_prod1);

        // weighting(:,:,t1,t2,ri) = weighting(:,:,t1,t2,ri) + shiftdim(sns_prod, 1);
        
        complex float *wt1t2 = wc + t1*w_strs1[KSP_DIMS] + t2*w_strs1[KSP_DIMS+1];
        md_zadd(DIMS, pe_dims, wt1t2, wt1t2, sns_prod1);

        loop_count++;
        if( n_loop > 10 && loop_count % (n_loop / 10 ) == 0 ){
            debug_printf(DP_DEBUG1, "\r%.0f%%", 100*(float) loop_count / 
                            (float) n_loop);
        }
    }

    //weighting = 1/(prod(sns_dims(1:3))^2) * weighting;
    md_zsmul(DIMS, w_dims, wc, wc, (complex float) 1.0 / pow((float) md_calc_size(2, &sns_dims[1]),2));

    md_free(sns_prod);
    md_free(sns_maps_perm);

    ts[3] = timestamp() - start_time;
	debug_printf(DP_INFO, "Timing: %f, %f, %f, %f\n", 
                    ts[0], ts[1], ts[2], ts[3]);
}

void compute_dd_fft(const int kDims, const long dimsp[], int *p, const int *pat){

    long p_size = md_calc_size(kDims+2, dimsp);
    long pat_size = md_calc_size(kDims+1, dimsp);
    long k_size = md_calc_size(kDims, dimsp);
    
    complex float *p_cfl = md_alloc(kDims+2, dimsp, sizeof(complex float));
    for( long i = 0 ; i < p_size ; ++i)
        p_cfl[i] = p[i];

    complex float *pat_cfl = md_alloc(kDims+1, dimsp, sizeof(complex float));
    for( long i = 0 ; i < pat_size ; ++i )
        pat_cfl[i] = pat[i];
    
    complex float *psf = md_alloc(kDims+1, dimsp, sizeof(complex float));
    
    // PSF = ifft(pat)
    assert(kDims == 2); // 3u
    ifft(kDims+1, dimsp, 3u, psf, pat_cfl); // needs to be scaled by 1/N
    /*
    for( long i = 0 ; i < pat_size ; i++ ){
        float *psfi = (float *) &psf[i];
        debug_printf(DP_INFO, "%f + i%f\n", psfi[0], psfi[1]);
    }
    */

    // Compute products of the PSF
    complex float *psf1psf2 = md_alloc(kDims, dimsp, sizeof(complex float));
    long nt = dimsp[kDims];
    long strs1[kDims+2];
    md_calc_strides(kDims+2, strs1, dimsp, 1);

    for( long t1 = 0 ; t1 < nt ; t1++ )
    for( long t2 = 0 ; t2 < nt ; t2++ ){
        md_zmulc(kDims, dimsp, psf1psf2, &psf[t1*k_size], &psf[t2*k_size]);
        assert(kDims == 2); // 3u
        fft(kDims, dimsp, 3u, &p_cfl[t1*strs1[kDims] + t2*strs1[kDims+1]], psf1psf2);
    }
    md_zsmul(kDims+2, dimsp, p_cfl, p_cfl, (complex float) 1.0 / (float) k_size); 

    // Round because p_cfl is within numerical precision of an array of integers
    for( long i = 0 ; i < p_size ; ++i ){
        p[i] = (int) (p_cfl[i] + 0.5);
    }
    
    md_free(psf1psf2);
    md_free(p_cfl);
    md_free(pat_cfl);
    md_free(psf);
}

/* 
 * D         = dimension
 * p         = differential domain matrix of size [dims, Nbins, Nbins], dims is size vector of length D
 * mask_dims = sampling mask dims [dims Nbins]
 * samples   = gridded samples
 * N         = number of samples in each bin
*/
void compute_dd_sparse_pat(const int D, const long dims[], double *p, 
        const long *samples[], 
        const long N[]){
    int i,j,bi,bj;

    long p_size = md_calc_size(D+2, dims);

    long strs1[D+2];
    md_calc_strides(D+2, strs1, dims, 1);

    memset(p, 0, sizeof(double)*p_size);
    assert(dims[D+0] == dims[D+1]);

    for( bi=0 ; bi < dims[D+0] ; bi++ )
    for( bj=0 ; bj < dims[D+1] ; bj++ )
    for( i =0 ; i < N[bi] ; i++ )
    for( j =0 ; j < N[bj] ; j++ ){
        long diff[D+2]; /* si - sj + nbins */
        subl(D, diff, &samples[bi][D*i], &samples[bj][D*j]);
        modl(D, diff, diff, dims);
        diff[D+0] = bi;
        diff[D+1] = bj;
        p[sub2idx(D+2, strs1, diff)]++;
    }
}

