#ifndef MISC_H
#define MISC_H 1

#include "mdarray.h"
#include "multind.h"
// #include <stdarg.h>
#include <stdlib.h>

#include <array>
#include <assert.h>
#include <iostream>

template <typename T>
T prod(long N, const T *a);

// dst = src1 + src2
template <typename T>
void add(long N, T *dst, const T *src1, const T *src2);

template <typename T>
T min(const long N, const T *x);

template <typename T>
long argmin(const long N, const T *x);

// dst = src1 - src2
template <typename T>
void sub(long N, T *dst, const T *src1, const T *src2);

// dst = mod(src1, src2)
template <typename T>
void mod(long N, T *dst, const T *src1, const T *src2);

// k-th largest element, k = 0 means largest element
template <typename T>
T kthLargest(const std::vector<T> &arr, long k);

template <typename T>
long sub2ind(const unsigned int D, const long *strides, const T *sub)
{
    long res = 0;
    for (unsigned int i = 0; i < D; i++)
    {
        res += sub[i] * strides[i];
    }
    return res;
}

template <typename T, size_t D>
long sub2ind(const std::array<T, D> &strides, const T *sub)
{
    return sub2ind(D, strides.data(), sub);
}

template <size_t D>
long sub2ind(const std::array<long, D> &strides, const std::array<long, D> &sub)
{
    return sub2ind(D, strides.data(), sub.data());
}

// convert linear index to MD subscript
template <typename T>
void ind2sub(const long D, const long *dims, T *subs, const long ind)
{
    long dims_prod = prod<long>(D - 1, dims);
    long ind_r = ind;
    for (long k = D - 1; k >= 0; --k)
    {
        if (ind_r == 0)
        {
            memset(subs, 0, sizeof(T) * k);
        }
        subs[k] = ind_r / dims_prod;
        if (subs[k] > 0)
        {
            ind_r -= subs[k] * dims_prod;
        }
        if (k > 0)
        {
            dims_prod /= dims[k - 1];
        }
    }
}

template <typename T, size_t D>
void ind2sub(const std::array<long, D> &dims, T *subs, const long ind)
{
    ind2sub(D, dims.data(), subs, ind);
}

#endif
