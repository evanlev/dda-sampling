#ifndef MISC_HPP
#define MISC_HPP 1

#include "misc.h"
#include <assert.h>
#include <iostream>
#include <numeric>      // std::accumulate
#include <functional>

// k-th largest but destroy input array
template<typename T>
T kthLargest2(const long N, T *arr, long k){
    assert( k < N );
    if( k == 0 && N == 1 ){
        return arr[0];
    }
    long left = 0;
    long right = N-1;
    long pivot = (right + left)/2;
    T pivotValue = arr[pivot];

    while( left <= right ){
        while( arr[left] < pivotValue ){
            left++;
        }
        while( arr[right] > pivotValue ){
            right--;
        }
        if( left <= right ){
            swap(&arr[left], &arr[right], sizeof(T));
            left++;
            right--;
        }
    }
    if( k < N - left ){
        return kthLargest2<T>(N - left, &arr[left], k);
    }else{
        return kthLargest2<T>(left, arr, k - N + left);
    }
}

// k-th largest element, k = 0 means largest element
template <typename T>
T kthLargest(const long N, const T *arr, long k){
    if( k >= N || k < 0 ){
        std::cout << "k: " << k << ", N: " << N << std::endl;
    }
    assert(k <  N);
    assert(k >= 0);

    T *tmp = new T[N];
    memcpy(tmp, arr, N*sizeof(T));

    return kthLargest2<T>(N, tmp, k);

    delete [] tmp;
}

template<typename T>
vector<T> md_calc_strides(vector<T> &dims){
    // Calculate strides 
    vector<long> strs(dims.size());
    long old = 1;
    for( unsigned int i = 0 ; i < dims.size() ; i++ ){
        strs[i] = (1== dims[i]) ? 0 : old;
        old *= dims[i];
    }
    return strs;
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
vector<T>
ind2subv(const vector<T> &dims, const long ind)
{
    vector<T> subs(dims.size());
    const long D = (signed) dims.size();
    long dims_prod = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<T>());
    long ind_r = ind;
    for( long k = D-1 ; k >= 0 ; --k ){
        subs[k] = ind_r / dims_prod;
        if( subs[k] > 0 ){
            ind_r -= subs[k] * dims_prod;
        }
        if( k > 0 ){
            dims_prod /= dims[k-1];
        }
    }
}

/* convert linear index to MD subscript */
template<typename T>
void
ind2sub(const long D, const long *dims, T *subs, const long ind)
{
    long dims_prod = prod<long>(D-1, dims); 
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

#endif // MISC_HPP

