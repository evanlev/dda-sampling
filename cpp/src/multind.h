#ifndef MULTIND_H
#define MULTIND_H 1

#include <string.h>
#include <array>

#include <stdbool.h>
#include <stdlib.h>

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) >= (Y) ? (X) : (Y))

typedef void (*md_nary_fun_t)(void* data, void* ptr[]);
typedef void (*md_loop_fun_t)(void* data, const long* pos);

#define MD_BIT(x) (1u << (x))
#define MD_IS_SET(x, y)	((x) & MD_BIT(y))

void md_nary(unsigned int C, unsigned int D, const long dim[], const long* str[], void* ptr[], void* data, md_nary_fun_t fun);
void md_loop(unsigned int D, const long dim[], void* data, md_loop_fun_t fun);
bool md_next(unsigned int D, const long dims[], unsigned int flags, long pos[]);
long md_calc_offset(unsigned int D, const long strides[], const long position[]);
long md_calc_size(unsigned int D, const long dim[]);
unsigned int md_calc_blockdim(unsigned int D, const long dim[], const long str[], size_t size);
void md_select_dims(unsigned int D, unsigned long flags, long odims[], const long idims[]);

void md_copy_dims(unsigned int D, long odims[], const long idims[]);
void md_copy_strides(unsigned int D, long ostrs[], const long istrs[]);
void md_set_dims(unsigned int D, long dims[], long val);
bool md_is_index(unsigned int D, const long pos[], const long dims[]);
bool md_check_dimensions(unsigned int N, const long dims[], unsigned int flags);
void md_singleton_dims(unsigned int D, long dims[]);
void md_singleton_strides(unsigned int D, long strs[]);
bool md_check_compat(unsigned int D, unsigned long flags, const long dim1[], const long dim2[]);
bool md_check_bounds(unsigned int D, unsigned long flags, const long dim1[], const long dim2[]);
void md_min_dims(unsigned int D, unsigned long flags, long odims[], const long idims1[], const long idims2[]);
void md_clear2(unsigned int D, const long dim[], const long str[], void* ptr, size_t size);
void md_clear(unsigned int D, const long dim[], void* ptr, size_t size);
void md_calc_strides(unsigned int D, long str[], const long dim[], size_t size);

void md_circular_swap2(unsigned M, unsigned int D, const long dims[], const long* strs[], void* ptr[], size_t size);
void md_circular_swap(unsigned M, unsigned int D, const long dims[], void* ptr[], size_t size);
void md_swap2(unsigned int D, const long dim[], const long ostr[], void* optr, const long istr[], void* iptr, size_t size);
void md_swap_flip2(unsigned int D, const long dims[], unsigned long flags, const long ostr[], void* optr, const long istr[], void* iptr, size_t size);
void md_copy2(unsigned int D, const long dim[], const long ostr[], void* optr, const long istr[], const void* iptr, size_t size);
void md_copy(unsigned int D, const long dim[], void* optr, const void* iptr, size_t size);
void md_fill2(unsigned int D, const long dim[], const long str[], void* ptr, const void* iptr, size_t size);
void md_fill(unsigned int D, const long dim[], void* ptr, const void* iptr, size_t size);
void md_flip2(unsigned int D, const long dims[], unsigned long flags, const long ostr[], void* optr, const long istr[], const void* iptr, size_t size);
void md_flip(unsigned int D, const long dims[], unsigned long flags, void* optr, const void* iptr, size_t size);
void* md_alloc(unsigned int D, const long dimensions[], size_t size);
void* md_calloc(unsigned int D, const long dimensions[], size_t size);
void* xmalloc(size_t s);

#define MD_INIT_ARRAY(x, y) { [ 0 ... ((x) - 1) ] = (y) } 

// C++ Wrappers
template <size_t D>
long md_calc_size(const std::array<long, D>& dims)
{
	return md_calc_size(D, dims.data());
}

template <size_t D>
std::array<long, D> md_calc_strides(const std::array<long, D>& dims, size_t size)
{
	std::array<long, D> str;
	md_calc_strides(D, str.data(), dims.data(), size);
	return str;
}

#endif // MULTIND_H
