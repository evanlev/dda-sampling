#include <getopt.h>
#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include <math.h>

#include "num/multind.h"
#include "num/flpmath.h"

#include "misc/misc.h"
#include "misc/mri.h"
#include "misc/debug.h"
#include "misc/mmio.h"
#include "misc/opts.h"

#include "dda_tools/dda_utils.h"

static const char* usage_str = 	"Generate a sampling pattern for sensitivity maps. Map dims are xres,yres,zres,coils,maps,time,maps2,\n"
                                "where the image is formed as DF{ sum_{k,l} sns_maps(:,:,:,:,k,:,l) m(:,:,:,1,k,1,l)}. It is sometimes \n"
                                "convenient to use both k and l (maps and maps2) indices (dimensions). Best candidate sampling is used, \n"
                                "and can be made faster with -T or -K options.";

static void int2cfl(const long N, complex float *dst, const int *src){
    for( long i = 0 ; i < N ; ++i )
        dst[i] = (complex float) src[i];
}

int main_dda_bc(int argc, char* argv[])
{
    int maxST = 0;
    int totSamps = 0;
    bool exact = false;
    float T = 0;
    int K = 0;
    const char* pat_file = NULL;
    const char* maps_file = NULL;

    struct arg_s args[] = {
        ARG_INFILE(true, &maps_file, "sensitivity maps"),
        ARG_OUTFILE(true, &pat_file, "pattern"),
    };

    const struct opt_s opts[] = {
        OPT_INT('d', &debug_level, "level", "debug level"),
        OPT_INT('t', &maxST, "maxST", "max samples per frame"),
        OPT_INT('M', &totSamps, "tot", "max total samples (optional)"),
        OPT_SET('e', &exact, "Exact BC sampling"),
        OPT_INT('K', &K, "K", "Limit the size of the support set of w to this"),
        OPT_FLOAT('T', &T, "T", "Set w(Delta k, t, t') < T to 0"),
    };

    cmdline(&argc, argv, ARRAY_SIZE(args), args, usage_str, ARRAY_SIZE(opts), opts);
	
    // Check options: see setup later
    if( (K != 0 || T != 0) && exact ){
        debug_printf(DP_ERROR, "For exact BC sampling, do not provide thresholding parameters -T or -K\n");
        exit(0);
    }
    if( K && T > 0 ){
        debug_printf(DP_ERROR, "Provide only -K or -T\n");
        exit(0);
    }
    if( maxST == 0 && totSamps == 0){
        debug_printf(DP_ERROR, "Need to specify either total samples (-M) or max samples per frame (-t)\n");
        exit(0);
    }
    if(totSamps){
        debug_printf(DP_INFO, "Total samples:         %d\n", totSamps);
    }
    if(maxST){
        debug_printf(DP_INFO, "Max samples per phase: %d\n", maxST);
    }

    
    // Read in sensitivity maps
    long sns_dims[DIMS];
    complex float *sns_maps = load_cfl(argv[1], DIMS, sns_dims);
    check_sns_dims(sns_dims);
    long nt = sns_dims[TE_DIM];

    // Convert sensitivity maps to w
    debug_printf(DP_INFO, "Compute w from sensitivity maps...\n");
    long w_dims[DIMS];
    get_w_dims(w_dims, sns_dims);
    double *w = md_alloc(DIMS, w_dims, sizeof(double));
    buildW(w, sns_dims, sns_maps);
    
    // Free sensitivty maps
    unmap_cfl(DIMS, sns_dims, sns_maps);

    // Allocate sampling pattern
    long pat_dims[DIMS];
    md_singleton_dims(DIMS, pat_dims);
    pat_dims[0] = sns_dims[1];
    pat_dims[1] = sns_dims[2];
    pat_dims[2] = nt;
    long pat_size = md_calc_size(DIMS, pat_dims);
    int *pat_int = xmalloc(pat_size*sizeof(int));
    
    // Allocate deltaJ = change in cost associated with sampling (ky,kz,t)
    double cost;
    debug_printf(DP_INFO, "Allocating Delta J, size %d\n", pat_size);
    double *deltaJ = xmalloc(pat_size * sizeof(double));

    // Set up constraints on sampling from options
    if(maxST == 0){
        maxST    = totSamps; // totSamps = no constraint
    }else if(totSamps == 0){
        totSamps = maxST * nt;
    }

    // Max samples per frame: same for all frames
    long maxSampsPerFrame[nt];
    for( long t = 0 ; t < nt ; t++ ){
        maxSampsPerFrame[t] = maxST;
    }

    // Generate the sampling pattern
    if( !exact ){
        // Set up wsp a thresholded version of w
        SparseW wsp;
        if( K ){
            sparsify_w_to_k(KSP_DIMS, w_dims, &wsp, w, (long) K);
        }else{
            sparsify_w(KSP_DIMS, w_dims, &wsp, w, T);
        }
        // Generate it
        approxBestCandidate(KSP_DIMS, &cost, deltaJ, pat_int, &wsp, 
                                     maxSampsPerFrame, totSamps);
        /*
        */
    }else{
        // Set up array with list of samples acquired
        long *samples[nt];
        for( long t = 0 ; t < nt ; t++ ){
            if( maxSampsPerFrame[t] > 0 ){
                samples[t] = xmalloc(KSP_DIMS * maxSampsPerFrame[t] * sizeof(long));
            }
        }

        // Generate it
        exactBestCandidate(KSP_DIMS, w_dims, &cost, deltaJ, pat_int, w, samples,
                            maxSampsPerFrame, totSamps);

        // Free array with list of samples acquired
        for( long t = 0 ; t < nt ; t++ ){
            if( maxSampsPerFrame[t] > 0 ){
                free(samples[t]);
            }
        }
    }

    // Convert output for complex
    debug_printf(DP_DEBUG3, "Converting output pattern\n");
    complex float *pat_cfl = create_cfl(pat_file, DIMS, pat_dims);
    int2cfl(pat_size, pat_cfl, pat_int);
    
    // cleanup
    debug_printf(DP_DEBUG3, "Clean up ...\n");
    unmap_cfl(DIMS, pat_dims, pat_cfl);
    md_free(w);
    free(deltaJ);
    free(pat_int);

    debug_printf(DP_DEBUG3, "Done.\n");
    
    return 0;
}




