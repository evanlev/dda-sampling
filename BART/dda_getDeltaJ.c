
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
#include "dda_tools/misc_dd.h"

static const char* usage_str = "<w> <pat> <Delta J>";
static const char* help_str = 	"Compute the change in cost associated with inserting a sample DeltaJ(k,t)\n";

static void dbl2cfl(const long N, complex float *dst, const double *src){
    for( long i = 0 ; i < N ; ++i )
        dst[i] = (complex float) src[i];
}
static void cfl2dbl(const long N, double *dst, const complex float *src){
    for( long i = 0 ; i < N ; ++i )
        dst[i] = (double) crealf(src[i]);
}

static void checkWAndPat(const long w_dims[], const long pat_dims[]){
    const char *err_msg = "w dimensions and pattern dimensions do not match.\nw should be Ny x Nz x nt x nt and pat should be Ny x Nz x nt\n";
    for( int i = 0 ; i < 3 ; i++ ){
        if( w_dims[i] != pat_dims[i] ){
            debug_printf(DP_ERROR, err_msg);
            exit(0);
        }
    }
    if( w_dims[2] != w_dims[3] ){
        debug_printf(DP_ERROR, err_msg);
        exit(0);
    }
}

int main_dda_getDeltaJ(int argc, char* argv[])
{
    bool conv_method = false;
    const struct opt_s opts[] = {
        OPT_INT('d', &debug_level, "level", "debug level"),
        OPT_SET('c', &conv_method, "convolution method (usually slower)"),
    };
    cmdline(&argc, argv, 3, 3, usage_str, help_str, ARRAY_SIZE(opts), opts);
    
    // Read in w and pattern
    long w_dims[DIMS];
    long pat_dims[DIMS];
    complex float *w = load_cfl(argv[1], DIMS, w_dims);
    complex float *pat = load_cfl(argv[2], DIMS, pat_dims);
    long nt = pat_dims[KSP_DIMS];

    // Check input
    checkWAndPat(w_dims, pat_dims);

    // Allocate output
    complex float *deltaJ = create_cfl(argv[3], DIMS, pat_dims);

    // Compute Delta J
    if( conv_method ){
        // Find samples
        long *samples[nt];
        long Nsamps[nt];
        find_samples(samples, Nsamps, pat, pat_dims, KSP_DIMS);

        // Convert to double
        double *deltaJd = md_alloc(DIMS, pat_dims, sizeof(double));
        double *wd = md_alloc(DIMS, w_dims, sizeof(double));
        cfl2dbl(md_calc_size(DIMS, w_dims), wd, w);

        computeDeltaJ(KSP_DIMS, w_dims, deltaJd, wd, samples, Nsamps);

        dbl2cfl(md_calc_size(DIMS, pat_dims), deltaJ, deltaJd);

        // Free double arrays
        md_free(deltaJd);
        md_free(wd);
    }else{
        // Usually faster
        computeDeltaJFftMethod(KSP_DIMS, w_dims, deltaJ, w, pat);
    }
    
    // cleanup
    unmap_cfl(DIMS, w_dims, w);
    unmap_cfl(DIMS, pat_dims, deltaJ);
    unmap_cfl(DIMS, pat_dims, pat);
    
    return 0;
}




