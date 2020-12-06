#ifndef MISC_HPP
#define MISC_HPP 1

#include "misc.h"
#include <array>
#include <assert.h>
#include <iostream>
#include <numeric>      // std::accumulate
#include <functional>
#include <algorithm>

// k-th largest but destroy input array
template<typename T>
T kthLargest2(std::vector<T> &arr, long k){
    assert(k <  arr.size());
    assert(k >= 0);

    std::nth_element(arr.begin(), arr.begin() + arr.size() - k, arr.end());

    return arr[arr.size() - k];
}

// k-th largest element, k = 0 means largest element
template <typename T>
T kthLargest(const std::vector<T> &arr, long k)
{
    std::vector<T> tmp = arr;
    return kthLargest2<T>(tmp, k);
}

template<typename T>
long
sub2ind(const unsigned int D, const vector<long> &strides, const vector<T> &sub){
    long res = 0;
    for( unsigned int i = 0 ; i < D ; i++ ){
        res += sub[i]*strides[i];
    }
    return res;
}


template<typename T>
long
sub2ind(const unsigned int D, const long *strides, const T *sub){
    long res = 0;
    for( unsigned int i = 0 ; i < D ; i++ ){
        res += sub[i]*strides[i];
    }
    return res;
}

/* convert linear index to MD subscript */
template<typename T>
void ind2sub(const long D, const long *dims, T *subs, const long ind)
{
    long dims_prod = prod<long>(D-1, dims); 
    // TODO redundant, but why D-1?

    long ind_r = ind;
    for( long k = D-1 ; k >= 0 ; --k ){
        if( ind_r == 0 ){
            memset(subs, 0, sizeof(T)*k);
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

template<typename T, size_t D>
void ind2sub(const std::array<long, D> &dims, T *subs, const long ind)
{
    ind2sub(D, dims.data(), subs, ind);
}

#endif // MISC_HPP

