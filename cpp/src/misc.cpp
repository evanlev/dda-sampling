// #include "debug.h"
#include "misc.h"

#include <algorithm>
#include <functional>

namespace
{
    template <typename T>
    void min_aux(long *min_el, T *min_val, const long N, const T *x)
    {
        *min_el = 0;
        *min_val = x[0];
        for (long i = 0; i < N; i++)
        {
            if (x[i] < *min_val)
            {
                *min_el = i;
                *min_val = x[i];
            }
        }
    }

    template <typename T>
    inline T mod1(const T x, const T y)
    {
        T r = x % y;
        return r < 0 ? r + y : r;
    }

    // k-th largest but destroy input array
    template <typename T>
    T kthLargest2(std::vector<T> &arr, long k)
    {
        assert(k < arr.size());
        assert(k >= 0);

        std::nth_element(arr.begin(), arr.begin() + arr.size() - k, arr.end());

        return arr[arr.size() - k];
    }
}

template <typename T>
T min(const long N, const T *x)
{
    long min_el;
    T min_val;
    min_aux<T>(&min_el, &min_val, N, x);
    return min_val;
}

template <typename T>
long argmin(const long N, const T *x)
{
    long min_el;
    T min_val;
    min_aux<T>(&min_el, &min_val, N, x);
    return min_el;
}

template <typename T>
void sub(long N, T *dst, const T *src1, const T *src2)
{
    long i;
    for (i = 0; i < N; i++)
    {
        dst[i] = src1[i] - src2[i];
    }
}

template <typename T>
T prod(long N, const T *a)
{
    T res = 1;
    for (long i = 0; i < N; i++)
    {
        res *= a[i];
    }
    return res;
}

template <typename T>
void add(long N, T *dst, const T *src1, const T *src2)
{
    long i;
    for (i = 0; i < N; i++)
    {
        dst[i] = src1[i] + src2[i];
    }
}

template <typename T>
void mod(long N, T *dst, const T *src1, const T *src2)
{
    for (long i = 0; i < N; i++)
    {
        dst[i] = mod1<T>(src1[i], src2[i]);
    }
}

template <typename T>
T kthLargest(const std::vector<T> &arr, long k)
{
    std::vector<T> tmp = arr;
    return kthLargest2<T>(tmp, k);
}

template void mod(long N, long *dst, const long *src1, const long *src2);
template void add(long N, long *dst, const long *src1, const long *src2);
template long prod(long N, const long *a);
template void sub(long N, long *dst, const long *src1, const long *src2);
template double min(const long N, const double *x);
template long argmin(const long N, const double *x);
