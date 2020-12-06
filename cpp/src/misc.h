#ifndef MISC_H
#define MISC_H 1

#include "mdarray.h"
#include "multind.h"
#include <stdarg.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <sstream>
//#include <ostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>

//#define DIMS 16u

using namespace std;

template<typename T>
T max(const long N, const T *x){
    long max_el = 0;
    T max_val = x[max_el];
    for( long i = 0 ; i < N ; i++ ){
        if( x[i] > max_val ){
            max_el = i;
            max_val = x[i];
        }
    }
    return max_val;
}

/* res = sum(src) */
template<typename T>
T sum(const long N, const T *src)
{
    T res = 0;
    for( long i = 0 ; i < N ; i++ ){
        res += src[i];
    }
    return res;
}

/* Sampling */  
std::vector<std::vector<long> > find_samples(const MDArray<3, int> &mask);

/* Vector operations */
void cross3(double *res, const double *x, const double *y);

double dot(const int D, const double *x, const double *y);

void unit_vector(const long N, double *x);

void axpy(const long N, double *res, const double alpha, const double *x, const double *y);

void smul(const long N, double *alphax, const double alpha, const double *x);
void smull(const long N, long *alphax, const double alpha, const long *x);

// product of an array
template<typename T>
T prod(long N, const T* a){
    T res = 1;
    for( long i = 0 ; i < N ; i++ )
        res *= a[i];
    return res;
}


/* add
 * dst = src1 + src2
 * */
template<typename T>
void add(long N, T *dst, const T *src1,  const T *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i] + src2[i];
    }
}


template<typename T>
void min_aux(long *min_el, T *min_val, const long N, const T *x){
    *min_el = 0;
    *min_val = x[0];
    for( long i = 0 ; i < N ; i++ ){
        if( x[i] < *min_val ){
            *min_el = i;
            *min_val = x[i];
        }
    }
}

template<typename T>
T min(const long N, const T *x){
    long min_el;
    T min_val;
    min_aux<T>(&min_el, &min_val, N, x);
    return min_val;
}

template<typename T>
long argmin(const long N, const T *x){
    long min_el;
    T min_val;
    min_aux<T>(&min_el, &min_val, N, x);
    return min_el;
}


/* sub 
 * dst = src1 - src2
 * */
template<typename T>
void sub(long N, T *dst, const T *src1,  const T *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i] - src2[i];
    }
}


/*
 * res = mod(x,y) 
 */
template<typename T>
T mod1(const T x, const T y){
    T r = x % y;
    return r < 0 ? r + y : r;
}

/*
 * dst = mod(src1, src2)
 */
template<typename T>
void mod(long N, T *dst, const T *src1, const T *src2){
    //cout << "in mod" << endl;
    for( long i =0 ; i < N ; i++ ){
        //cout << "mod(" << src1[i] << ", " << src2[i] << ") = ";
        dst[i] = mod1<T>(src1[i], src2[i]);
        //cout << dst[i] << endl;
    }
}



void modd(long N, double *dst, const double *src1, const double *src2);
double mod2(const double x, const double y);

void mul(long N, double *dst, const double *src1,  const double *src2);

void* xcalloc(size_t s);

double normp(const long D, double *x, float p);


double norm(const long D, double *x);

/* MD array operations */
/*
void md_copy_dims(const unsigned int D, long *dims, const long *dims1);
void md_clear( const int D, void *p, const long *p_dims, size_t size);

*/


int cmpfunc_double (const void * a, const void * b);

int compareToDouble(const double a, const double b);

/* Randomization */
void rpermute(const long n, long *a);
void randperm( const long n, long perm[]);
void draw_uniform(const unsigned int D, long *new_s, const long *dims);

#endif
