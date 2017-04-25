
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

static const char* usage_str = "<input maps> <output w>";
static const char* help_str = 	"Convert sensitivity maps to a weighting w\n"
                                "maps have dimensions [nx ny nz nc nt nl]\n"
                                "The forward model is \n"
                                "y(k,t,c,l,m) = D_t F{ sum_{l,m} S(r,t,c,l,m) x_{l,m}},\n"
                                "where D_t are sampling operators, S are maps, x are images\n";

int main_dda_getw(int argc, char* argv[])
{
    const struct opt_s opts[] = {
        OPT_INT('d', &debug_level, "level", "debug level"),
    };
    cmdline(&argc, argv, 2, 2, usage_str, help_str, ARRAY_SIZE(opts), opts);
    
    // Read in sensitivity maps
    long sns_dims[DIMS];
    complex float *sns_maps = load_cfl(argv[1], DIMS, sns_dims);

    // Check input
    check_sns_dims(sns_dims);

    // Convert sensitivity maps to w
    long w_dims[DIMS];
    get_w_dims(w_dims, sns_dims);
    double *w = md_alloc(DIMS, w_dims, sizeof(double));
    complex float *wc = create_cfl(argv[2], DIMS, w_dims);
    buildW2(wc, sns_dims, sns_maps);
    
    // cleanup
    md_free(w);
    unmap_cfl(DIMS, w_dims, wc);
    unmap_cfl(DIMS, sns_dims, sns_maps);
    
    return 0;
}




