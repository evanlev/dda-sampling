
#include <string.h> // memcpy
#include "misc/debug.h"
#include "misc/misc.h"
#include "num/multind.h"

#include "misc_dd.h"
#include <assert.h>

// swap two non-overlapping blocks of memory of size s
void swap(void *p1, void *p2, size_t s){
    unsigned char p3[s];
    //void *p3 = xmalloc(s);
    memcpy(p3, p2, s); // p3 = p2
    memcpy(p2, p1, s); // p2 = p1
    memcpy(p1, p3, s); // p1 = p3
    //free(p3);
}

// k-th largest but destroy input array
static double kthLargest2(const long N, double *arr, long k){
    assert( k < N );
    if( k == 0 && N == 1 ){
        return arr[0];
    }
    long left = 0;
    long right = N-1;
    long pivot = (right + left)/2;
    double pivotValue = arr[pivot];

    while( left <= right ){
        while( arr[left] < pivotValue ){
            left++;
        }
        while( arr[right] > pivotValue ){
            right--;
        }
        if( left <= right ){
            swap(&arr[left], &arr[right], sizeof(double));
            left++;
            right--;
        }
    }
    if( k < N - left ){
        return kthLargest2(N - left, &arr[left], k);
    }else{
        return kthLargest2(left, arr, k - N + left);
    }
}

// k-th largest element, k = 0 means largest element
double kthLargest(const long N, const double *arr, long k){
    if( k >= N || k < 0 ){
        debug_printf(DP_ERROR, "k: %d, N: %d\n", k, N);
        exit(0);
    }
    assert(k <  N);
    assert(k >= 0);

    double *tmp = xmalloc(N*sizeof(double));
    memcpy(tmp, arr, N*sizeof(double));

    return kthLargest2(N, tmp, k);

    free(tmp);
}

static void rpermute(const long n, long *a) {
    long k;
    for (k = 0; k < n; k++)
        a[k] = k;
    for (k = n-1; k > 0; k--) {
        long j = rand() % (k+1);
        long temp = a[j];
        a[j] = a[k];
        a[k] = temp;
    }
}

long
sub2idx(const unsigned int D, const long *strides, const long *sub){
    long res = 0;
    for( unsigned int i = 0 ; i < D ; i++ ){
        res += sub[i]*strides[i];
    }
    return res;
}


/* convert linear index to MD subscript */
void
idx2sub(const long D, const long *dims, long *subs, const long ind)
{
    long dims_prod = prodl(D-1, dims); 
    long ind_r = ind;
    for( long k = D-1 ; k >= 0 ; --k ){
        if( ind_r == 0 ){
            memset(subs, 0, sizeof(long)*k);
        }
        subs[k] = ind_r / dims_prod;
        if( subs[k] > 0 ){
            ind_r -= subs[k] * dims_prod;
        }
        if( k > 0 ){
            dims_prod /= dims[k-1];
        }
    }
}


/* randomly permutes the elements of the array perm of length n */
void randomperm( const long n, long perm[])
{
    long *ind = xmalloc(n*sizeof(long));
    long *tmp = xmalloc(n*sizeof(long));
    long k;
    rpermute(n, ind);
    for( k = 0 ; k < n ; k++ ){
        tmp[k] = perm[ind[k]];
    }
    for( k = 0 ; k < n ; k++ ){
        perm[k] = tmp[k];
    }
    free(ind);
    free(tmp);
}

/* Sampling */
static void find_samples_1frame(long *samples, const complex float *mask, const long mask_dims[], const int D){
    long nksp_pts = md_calc_size(D, mask_dims);
    long tail = 0;
    for( long ik = 0 ; ik < nksp_pts ; ik++ ){
        if( 1.0 == mask[ik] ){
            long ik_sub[D];
            idx2sub(D, mask_dims, ik_sub, ik);
            for( int d = 0 ; d < D ; d++ ){
                samples[D*tail + d] = ik_sub[d];
            }
            tail++;
        }
    }
}

/* Sampling */
void find_samples(long *samples[], long Nsamps[], const complex float *mask, const long mask_dims[], const int nksp_dims){
    long Nt = mask_dims[nksp_dims];
    long nksp_pts = md_calc_size(nksp_dims, mask_dims);
    for( int t = 0 ; t < Nt ; t++ ){
        const complex float* maskb = &mask[nksp_pts*t];
        Nsamps[t] = (long) sumcf(nksp_pts, maskb);
        if( Nsamps[t] > 0 ){
            samples[t] = xmalloc(sizeof(long)*Nsamps[t]*nksp_dims);
            find_samples_1frame(samples[t], maskb, mask_dims, nksp_dims);
        }
    }
}



// -------- vector operations

complex float sumcf(const long N, const complex float *x){
    complex float res = 0;
    for( long i = 0 ; i < N ; ++i ){
        res += x[i];
    }
    return res;
}


// product of an array
long prodl(long N, const long* a){
    long res = 1;
    for( long i = 0 ; i < N ; i++ )
        res *= a[i];
    return res;
}


/* add
 * dst = src1 + src2
 * */
void addd(long N, double *dst, const double *src1,  const double *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i] + src2[i];
    }
}


/* add
 * dst = src1 + src2
 * */
void addl(long N, long *dst, const long *src1,  const long *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i] + src2[i];
    }
}

/* res[i] = alpha * x[i]
*/
void smuld(const long N, double *alphax, const double alpha, const double *x){
    for( int i = 0 ; i < N ; i++ ){
        alphax[i] = x[i]*alpha;
    }
}

static void min_auxd(long *min_el, double *min_val, const long N, const double *x){
    *min_el = 0;
    *min_val = x[0];
    for( long i = 0 ; i < N ; i++ ){
        if( x[i] < *min_val ){
            *min_el = i;
            *min_val = x[i];
        }
    }
}

double mind(const long N, const double *x){
    long min_el;
    double min_val;
    min_auxd(&min_el, &min_val, N, x);
    return min_val;
}


long argmind(const long N, const double *x){
    long min_el;
    double min_val;
    min_auxd(&min_el, &min_val, N, x);
    return min_el;
}

/* sub 
 * dst = src1 - src2
 * */
void subl(long N, long *dst, const long *src1,  const long *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i] - src2[i];
    }
}

static long mod1l(const long x, const long y){
    long r = x % y;
    return r < 0 ? r + y : r;
}

/*
 * dst = mod(src1, src2)
 */
void modl(long N, long *dst, const long *src1, const long *src2){
    //cout << "in mod" << endl;
    for( long i =0 ; i < N ; i++ ){
        //cout << "mod(" << src1[i] << ", " << src2[i] << ") = ";
        dst[i] = mod1l(src1[i], src2[i]);
        //cout << dst[i] << endl;
    }
}




// ---- Debug
int is_in_bounds(const unsigned int D, const long *d, const long *dims){
    unsigned int i;
    for( i = 0 ; i < D ;i++)
    {
        if( d[i] < 0 )
            return 0;
        if( d[i] >= dims[i] )
            return 0;
    }
    return 1;
}
void assert_same_dims(const unsigned int D, const char *d1str, const long *d1, const char *d2str, const long *d2){
    for( unsigned int i = 0 ; i < D ; i++ ){
        if(d1[i] != d2[i] ){
            debug_printf(DP_ERROR, "%s dims: ", d1str);
            debug_print_dims(DP_ERROR, D, d1);
            debug_printf(DP_ERROR, "\n");
            debug_printf(DP_ERROR, "%s dims: ", d2str);
            debug_print_dims(DP_ERROR, D, d2);
            debug_printf(DP_ERROR, "\n");
            exit(0);
        }
    }

}

void assert_is_binary(const long N, const int *arr, const char *arrstr){
    for( long i = 0 ; i < N ; i++ ){
        if( arr[i] != 0 && arr[i] != 1){
            debug_printf(DP_ERROR, "%s is not binary", arrstr);
            exit(0);
        }
    }
}

void assert_in_bounds(const unsigned int D, const long *d, const long *dims, const char *str){
    if(!is_in_bounds(D, d, dims)){
        debug_printf(DP_ERROR, "%s=", str);
        debug_print_dims(DP_ERROR, D, d);
        debug_printf(DP_ERROR, " is out of bounds ");
        debug_print_dims(DP_ERROR, D, dims);
        exit(0);
    }
}

void myAssert( int exp, const char *err){
    if( !exp ){
        debug_printf(DP_ERROR, err);
        assert(0);
    }
}


