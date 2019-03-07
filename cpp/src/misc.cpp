
#include <stdarg.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <math.h>
#include <complex.h>

//#include "multind.h"
#include "debug.h"
#include "misc.h"
#include "misc.hpp"
#include "config.h"

//using std::complex;
    
void debug_print_arr(const int level, const long N, const double *arr){
    debug_printf(level, "[");
    for( long i = 0 ; i < N ; i++ ){
        debug_printf(level, "%f ", arr[i]);
    }
    debug_printf(level, "]\n");
}

/*
 * ----- Vector operations
 */

/*
 * res = cross product of x and y
*/
void cross3(double *res, const double *x, const double *y){
    res[0] = x[1]*y[2] - x[2]*y[1];
    res[1] = x[2]*y[0] - x[0]*y[2];
    res[2] = x[0]*y[1] - x[1]*y[0];
}

/*
 * = x'*y
 */
double dot(const int D, const double *x, const double *y){
   double res = 0;
   for( int i = 0 ; i < D ; ++i ){
        res += x[i]*y[i];
   }
   return res;
}

/*
 * res = mod(x,y) (double)
 */
double mod2(const double x, const double y){
    double r = remainder(x,y);
    return r < 0 ? r + y : r;
}

double normp(const long D, double *x, float p){
    long i = 0;
    double retval = 0;
    for( i = 0 ; i < D ; i++ ){
        retval += powf(x[i], p);
    }
    return powf(retval, 1.0 / p);
}

/*
 * res = ||x||
 */
double norm(const long D, double *x)
{
    return sqrt(dot(D, x, x));
}

/*
 * x = { x / ||x||2 , x != 0
 *     { 0          , x == 0
 */
void unit_vector(const long N, double *x){
    double len = norm(N, x);
    if( len > 0.0 ){
        int i;
        for( i =0 ; i < N ; i++ ){
            x[i] = x[i] / len;
        }
    }
}


/* res[i] = alpha x[i] + y[i]
*/
void axpy(const long N, double *res, const double alpha, const double *x, const double *y){
    for( int i = 0 ; i < N ; i++ ){
        res[i] = alpha*x[i] + y[i];
    }
}

/* res[i] = alpha * x[i]
*/
void smull(const long N, long *alphax, const double alpha, const long *x){
    for( int i = 0 ; i < N ; i++ ){
        alphax[i] = x[i]*alpha;
    }
}

/* res[i] = alpha * x[i]
*/
void smul(const long N, double *alphax, const double alpha, const double *x){
    for( int i = 0 ; i < N ; i++ ){
        alphax[i] = x[i]*alpha;
    }
}


/*
 * dst = mod(src1, src2)
 */
void modd(long N, double *dst, const double *src1, const double *src2){
    
    for( long i =0 ; i < N ; i++ ){
        dst[i] = mod2(src1[i], src2[i]);
    }
}

/* mul */
void mul(long N, double *dst, const double *src1,  const double *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i]*src2[i];
    }
}

int cmpfunc_double (const void * a, const void * b)
{
      if (*(double*)a > *(double*)b) return 1;
        else if (*(double*)a < *(double*)b) return -1;
          else return 0;  
}



void* xcalloc(size_t s)
{
    void *p = (void *) xmalloc(s);
    memset(p, 0, s);
    return p;
}


/*
// swap two non-overlapping blocks of memory of size s
void swap(void *p1, void *p2, size_t s){
    void *p3 = (void *) ::operator new(s);
    memcpy(p3, p2, s); // p3 = p2
    memcpy(p2, p1, s); // p2 = p1
    memcpy(p1, p3, s); // p1 = p3
    operator delete(p3);
}
*/

/*
 * Draw uniformly at random new_s in the range [0, dims-1]
 */
void draw_uniform(const unsigned int D, long *new_s, const long *dims){
    for( unsigned int i = 0 ; i < D ; i++)
        new_s[i] = rand() % dims[i];
}


