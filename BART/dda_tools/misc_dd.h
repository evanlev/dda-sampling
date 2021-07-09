#ifndef MISC_DD_H 
#define MISC_DD_H

#include <stdlib.h>
#include <complex.h>

extern void swap(void *p1, void *p2, size_t s);
extern double kthLargest(const long N, const double *arr, long k);

extern void randomperm( const long n, long perm[]);

extern void find_samples(long *samples[], long Nsamps[], const complex float *mask, const long mask_dims[], const int nksp_dims);

extern complex float sumcf(const long N, const complex float *x);

extern long prodl(long N, const long* a);
extern void addl(long N, long *dst, const long *src1,  const long *src2);
extern void addd(long N, double *dst, const double *src1,  const double *src2);
extern void smuld(const long N, double *alphax, const double alpha, const double *x);
extern double mind(const long N, const double *x);
extern long argmind(const long N, const double *x);
extern void subl(long N, long *dst, const long *src1,  const long *src2);
extern void modl(long N, long *dst, const long *src1, const long *src2);

extern long
sub2idx(const unsigned int D, const long *strides, const long *sub);

extern void
idx2sub(const long D, const long *dims, long *subs, const long ind);

extern int is_in_bounds(const unsigned int D, const long *d, const long *dims);

extern void assert_same_dims(const unsigned int D, const char *d1str, const long *d1, const char *d2str, const long *d2);

extern void assert_is_binary(const long N, const int *arr, const char *arrstr);

extern void assert_in_bounds(const unsigned int D, const long *d, const long *dims, const char *str);

extern void myAssert( int exp, const char *err);

#endif // MISC_DD_H
