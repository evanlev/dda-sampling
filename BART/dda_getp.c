
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

static const char* usage_str = 	"Get the differential distribution for a sampling pattern\n";

static void cfl2int(const long N, int *dst, const complex float *src){
    for( long i = 0 ; i < N ; ++i )
        dst[i] = (int) src[i];
}

int main_dda_getp(int argc, char* argv[])
{
    const char* pat_file = NULL;
    const char* p_file = NULL;

    struct arg_s args[] = {
        ARG_INFILE(true, &pat_file, "pattern"),
        ARG_OUTFILE(true, &p_file, "p"),
    };

    const struct opt_s opts[] = {
        OPT_INT('d', &debug_level, "level", "debug level"),
    };

    cmdline(&argc, argv, ARRAY_SIZE(args), args, usage_str, ARRAY_SIZE(opts), opts);

    const long D = 2; // 2 k-space dimensions
    
    long dims[DIMS];
    complex float* pat = load_cfl(pat_file, DIMS, dims);

    for (unsigned int d = 3; d < DIMS; ++d)
    {
        if (dims[d] > 1)
        {
            debug_printf(DP_ERROR, "pattern must be a 3D array with dimensions ky kz t\n");
            exit(0);
        }
    }

    long dimsp[DIMS];
    md_copy_dims(DIMS,dimsp, dims);
    dimsp[D] = dimsp[D+1] = dims[D];

    complex float *diff_dist = create_cfl(p_file, DIMS, dimsp);

    long p_size = md_calc_size(DIMS, dimsp);
    long pat_size = md_calc_size(DIMS, dims);

    // Convert to int
    int *diff_dist_int = xmalloc(p_size*sizeof(int));
    cfl2int(p_size, diff_dist_int, diff_dist);

    int *pat_int = xmalloc(pat_size*sizeof(int));
    cfl2int(pat_size, pat_int, pat);

    // Compute differential distribution: prefer to use int as in/out type
    compute_dd_fft(2, dimsp, diff_dist_int, pat_int);

    // Copy to output
    for( long i = 0 ; i < p_size ; i++ )
    {
        diff_dist[i] = (complex float) diff_dist_int[i];
    }
    
    // cleanup
    unmap_cfl(DIMS, dims, pat);
    unmap_cfl(DIMS, dims, diff_dist);
    free(pat_int);
    free(diff_dist_int);
    
    return 0;
}




