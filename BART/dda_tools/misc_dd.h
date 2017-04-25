#ifndef MISC_DD_H 
#define MISC_DD_H 1

#include <stdlib.h>
#include <complex.h>

void swap(void *p1, void *p2, size_t s);
double kthLargest(const long N, const double *arr, long k);

void randomperm( const long n, long perm[]);

void find_samples(long *samples[], long Nsamps[], const complex float *mask, const long mask_dims[], const int nksp_dims);

complex float sumcf(const long N, const complex float *x);

long prodl(long N, const long* a);
void addl(long N, long *dst, const long *src1,  const long *src2);
void addd(long N, double *dst, const double *src1,  const double *src2);
void smuld(const long N, double *alphax, const double alpha, const double *x);
double mind(const long N, const double *x);
long argmind(const long N, const double *x);
void subl(long N, long *dst, const long *src1,  const long *src2);
void modl(long N, long *dst, const long *src1, const long *src2);

long
sub2idx(const unsigned int D, const long *strides, const long *sub);

void
idx2sub(const long D, const long *dims, long *subs, const long ind);

int is_in_bounds(const unsigned int D, const long *d, const long *dims);

void assert_same_dims(const unsigned int D, const char *d1str, const long *d1, const char *d2str, const long *d2);

void assert_is_binary(const long N, const int *arr, const char *arrstr);

void assert_in_bounds(const unsigned int D, const long *d, const long *dims, const char *str);

void myAssert( int exp, const char *err);

#endif
