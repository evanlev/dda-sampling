
#include <string.h>
#include <stdbool.h>
#include <alloca.h>
#include <strings.h>
#include <stdio.h>

#include "debug.h"

#include "multind.h"


void* xmalloc(size_t s)
{
    /*
    mexPrintf("Allocating %d\n", s);
    mexEvalString("drawnow");
	*/
    // 10 GB limit
    if( s > 1e10 || s == 0 ){
        debug_printf(DP_INFO, "Tried to allocate %d bytes. Quitting\n", s);
        exit(0);
    }
	
#ifdef MEX_COMPILE_FLAG
    void* p = mxMalloc(s);
#else
    void* p = malloc(s);
#endif
	if (NULL == p){
        debug_printf(DP_ERROR, "Could not allocate memory");
        exit(0);
    }

    return p;
}


/**
 * Generic functions which loops over all dimensions of a set of
 * multi-dimensional arrays and calls a given function for each position.
 */
void md_nary(unsigned int C, unsigned int D, const long dim[], const long* str[], void* ptr[], void* data, md_nary_fun_t fun)
{
	if (0 == D) {

		fun(data, ptr);
		return;
	}

	for (long i = 0; i < dim[D - 1]; i++) {

		void* moving_ptr[C];
	
		for (unsigned int j = 0; j < C; j++)
			moving_ptr[j] = (char *) ptr[j] + i * str[j][D - 1];

		md_nary(C, D - 1, dim, str, moving_ptr, data, fun);
	}
}


static void md_loop_r(unsigned int D, const long dim[], long pos[], void* data, md_loop_fun_t fun)
{
	if (0 == D) {

		fun(data, pos);
		return;
	}

	D--;

	for (pos[D] = 0; pos[D] < dim[D]; pos[D]++) 
		md_loop_r(D, dim, pos, data, fun);
}

/**
 * Generic function which loops over all dimensions and calls a given
 * function passing the current indices as argument.
 *
 * Runs fun( data, position ) for all position in dim
 *
 */
void md_loop(unsigned int D, const long dim[], void* data, md_loop_fun_t fun)
{
	long pos[D];
	md_loop_r(D, dim, pos, data, fun);
}



/**
 * Computes the next position. Returns true until last index.
 */
bool md_next(unsigned int D, const long dims[], unsigned int flags, long pos[])
{
	if (0 == D--)
		return false;

	if (md_next(D, dims, flags, pos))
		return true;

	if (MD_IS_SET(flags, D)) {

		assert((0 <= pos[D]) && (pos[D] < dims[D]));

		if (++pos[D] < dims[D])
			return true;

		pos[D] = 0;
	}

	return false;
}



/**
 * Returns offset for position in a multidimensional array
 *
 * return pos[0]*strides[0] + ... + pos[D-1]*strides[D-1]
 *
 * @param D number of dimensions
 * @param dim dimensions array
 */
long md_calc_offset(unsigned int D, const long strides[], const long position[])
{
	long pos = 0;

	for (unsigned int i = 0; i < D; i++)
		pos += strides[i] * position[i];

	return pos;
}



static long md_calc_size_r(unsigned int D, const long dim[], size_t size)
{
	if (0 == D)
		return size;

	return md_calc_size_r(D - 1, dim, size * dim[D - 1]);
}

/**
 * Returns the number of elements
 *
 * return dim[0]*dim[1]*...*dim[D-1]
 *
 * @param D number of dimensions
 * @param dim dimensions array
 */
long md_calc_size(unsigned int D, const long dim[])
{
	return md_calc_size_r(D, dim, 1);
}


/**
 * Computes the number of smallest dimensions which are stored
 * contineously, i.e. can be accessed as a block of memory.
 * 
 */
unsigned int md_calc_blockdim(unsigned int D, const long dim[], const long str[], size_t size)
{
	long dist = size;
	unsigned int i = 0;

	for (i = 0; i < D; i++) {

		if (!((str[i] == dist) || (dim[i] == 1)))
			break;

		dist *= dim[i];
	}

	return i;
}	



/**
 * Copy dimensions specified by flags and set remaining dimensions to 1
 *
 * odims = [ 1  idims[1]  idims[2]  1  1  idims[5] ]
 *
 * @param D number of dimensions
 * @param flags bitmask specifying which dimensions to copy
 * @param odims output dimensions
 * @param idims input dimensions
 */
void md_select_dims(unsigned int D, unsigned long flags, long odims[], const long idims[])
{
	md_copy_dims(D, odims, idims);
	
	for (unsigned int i = 0; i < D; i++)
		if (!MD_IS_SET(flags, i))
			odims[i] = 1;
}


/**
 * Copy dimensions
 *
 * odims[i] = idims[i]
 */
void md_copy_dims(unsigned int D, long odims[], const long idims[])
{
	memcpy(odims, idims, D  * sizeof(long));
}


/**
 * Copy strides
 *
 * ostrs[i] = istrs[i]
 */
void md_copy_strides(unsigned int D, long ostrs[], const long istrs[])
{
	memcpy(ostrs, istrs, D  * sizeof(long));
}

/**
 * Set all dimensions to value
 *
 * dims[i] = val
 */
void md_set_dims(unsigned int D, long dims[], long val)
{
	for (unsigned int i = 0; i < D; i++)
		dims[i] = val;
}


/**
 * returns whether or not @param pos is a valid index of an array of dimension @param dims
 */
bool md_is_index(unsigned int D, const long pos[], const long dims[])
{
	if (D == 0)
		return true;

	return ((pos[0] >= 0) && (pos[0] < dims[0]) && md_is_index(D - 1, pos + 1, dims + 1));
}



/**
 * return whether some other dimensions are >1
 */
bool md_check_dimensions(unsigned int N, const long dims[], unsigned int flags)
{
	long d[N];
	md_select_dims(N, ~flags, d, dims);
	return (1 != md_calc_size(N, d));
}



/**
 * Set all dimensions to one
 *
 * dims[i] = 1
 */
void md_singleton_dims(unsigned int D, long dims[])
{
	for (unsigned int i = 0; i < D; i++)
		dims[i] = 1;
}



/**
 * Set all strides to one
 *
 * dims[i] = 1
 */
void md_singleton_strides(unsigned int D, long strs[])
{
	for (unsigned int i = 0; i < D; i++)
		strs[i] = 0;
}



/**
 * Check dimensions for compatibility. Dimensions must be equal or
 * where indicated by a set bit in flags one must be equal to one
 * in atleast one of the arguments.
 */
bool md_check_compat(unsigned int D, unsigned long flags, const long dim1[], const long dim2[])
{
	if (0 == D)
		return true;

	D--;

	if ((dim1[D] == dim2[D]) || (MD_IS_SET(flags, D) && ((1 == dim1[D]) || (1 == dim2[D]))))
		return md_check_compat(D, flags, dim1, dim2);
		
	return false;
}



/**
 * dim1 must be bounded by dim2 where a bit is set
 */
bool md_check_bounds(unsigned int D, unsigned long flags, const long dim1[], const long dim2[])
{
	if (0 == D--)
		return true;

	if (!MD_IS_SET(flags, D) || (dim1[D] <= dim2[D]))
		return md_check_bounds(D, flags, dim1, dim2);

	return false;
}


/**
 * Set the output's flagged dimensions to the minimum of the two input dimensions
 *
 * odims = [ MIN(idims1[0],idims2[0] ... MIN(idims1[D-1],idims2[D-1]) ]
 *
 * @param D number of dimensions
 * @param flags bitmask specifying which dimensions to minimize
 * @param odims output dimensions
 * @param idims1 input 1 dimensions
 * @param idims2 input 2 dimensions
 */
void md_min_dims(unsigned int D, unsigned long flags, long odims[], const long idims1[], const long idims2[])
{
	for (unsigned int i = 0; i < D; i++)
		if (MD_IS_SET(flags, i))
			odims[i] = MIN(idims1[i], idims2[i]);
}



struct data_s {

	size_t size;
#ifdef USE_CUDA
	bool use_gpu;
#endif
};

static void nary_clear(void* _data, void* ptr[])
{
	struct data_s* data = (struct data_s*)_data;
#ifdef  USE_CUDA
	if (data->use_gpu) {

		cuda_clear(data->size, ptr[0]);
		return;
	}
#endif
	memset(ptr[0], 0, data->size);	
}

/**
 * Zero out array (with strides)
 *
 * ptr[i] = 0
 */
void md_clear2(unsigned int D, const long dim[], const long str[], void* ptr, size_t size)
{
	int skip = md_calc_blockdim(D, dim, str, size);
#ifdef  USE_CUDA
	struct data_s data = { md_calc_size(skip, dim) * size, cuda_ondevice(ptr) };
#else
	struct data_s data = { md_calc_size(skip, dim) * size };
#endif
    const long *tmp[1] = {str + skip};
    void *tmp2[1] = {ptr};
	md_nary(1, D - skip, dim + skip, tmp, tmp2, (void*)&data, &nary_clear);
	//md_nary(1, D - skip, dim + skip, (const long*[1]){ str + skip }, (void*[1]){ ptr }, (void*)&data, &nary_clear);
}

/**
 * Calculate strides in column-major format 
 * (smallest index is sequential)
 *
 * @param D number of dimensions
 * @param array of calculates strides
 * @param dim array of dimensions
 * @param size of a single element
 */
void md_calc_strides(unsigned int D, long str[], const long dim[], size_t size)
{
	long old = size;

	for (unsigned int i = 0; i < D; i++)
	{
		str[i] = (1 == dim[i]) ? 0 : old;
		old *= dim[i];
	}
}

/**
 * Zero out array (without strides)
 *
 * ptr[i] = 0
 *
 * @param D number of dimensions
 * @param dim dimensions array
 * @param ptr pointer to data to clear
 * @param size sizeof()
 */
void md_clear(unsigned int D, const long dim[], void* ptr, size_t size)
{
	long str[D];
	md_calc_strides(D, str, dim, size);
	md_clear2(D, dim, str, ptr, size);
}



struct strided_copy_s {

	long sizes[2];
	long ostr;
	long istr;
};

#ifdef USE_CUDA
static void nary_strided_copy(void* _data, void* ptr[])
{
	struct strided_copy_s* data = _data;

	cuda_memcpy_strided(data->sizes, data->ostr, ptr[0], data->istr, ptr[1]);
}
#endif

static void nary_copy(void* _data, void* ptr[])
{
	struct data_s* data = (struct data_s*)_data;
#ifdef  USE_CUDA
	if (data->use_gpu) {


		cuda_memcpy(data->size, ptr[0], ptr[1]);
		return;
	}
#endif
	memcpy(ptr[0], ptr[1], data->size);	
}



/**
 * Swap values between two arrays (with strides)
 *
 * iptr[i] = optr[i] and optr[i] = iptr[i]
 */
void md_swap2(unsigned int D, const long dim[], const long ostr[], void* optr, const long istr[], void* iptr, size_t size)
{
    const long *tmp1[2] = {ostr, istr};
    void *tmp2[2] = {optr, iptr};
	md_circular_swap2(2, D, dim, tmp1, tmp2, size);
	//md_circular_swap2(2, D, dim, (const long*[2]){ ostr, istr }, (void*[2]){ optr, iptr }, size);
}


/**
 * Allocate CPU memory
 *
 * return pointer to CPU memory
 */
void* md_alloc(unsigned int D, const long dimensions[], size_t size)
{
	return xmalloc(md_calc_size(D, dimensions) * size);
}



/**
 * Allocate CPU memory and clear
 *
 * return pointer to CPU memory
 */
void* md_calloc(unsigned int D, const long dimensions[], size_t size)
{
	void* ptr = md_alloc(D, dimensions, size);
	md_clear(D, dimensions, ptr, size);
	return ptr;
}



/**
 * Swap input and output while flipping selected dimensions
 * at the same time.
 */
void md_swap_flip2(unsigned int D, const long dims[], unsigned long flags, const long ostr[], void* optr, const long istr[], void* iptr, size_t size)
{
#if 1
	int i;
	for (i = D - 1; i >= 0; i--)
		if ((1 != dims[i]) && MD_IS_SET(flags, i))
			break;

	if (-1 == i) {

		md_swap2(D, dims, ostr, optr, istr, iptr, size);
		return;
	}

	assert(1 < dims[i]);
	assert(ostr[i] != 0);
	assert(istr[i] != 0);

	long dims2[D];
	md_copy_dims(D, dims2, dims);
	dims2[i] = dims[i] / 2;

	long off = (dims[i] + 1) / 2;
	assert(dims2[i] + off == dims[i]);

	md_swap_flip2(D, dims2, flags, ostr, optr, istr, (char *) iptr + off * istr[i], size);
	md_swap_flip2(D, dims2, flags, ostr, (char *) optr + off * ostr[i], istr, iptr, size);

	dims2[i] = 1;

	if (1 == dims[i] % 2)
		md_swap_flip2(D, dims2, flags, ostr, (char *) optr + (off - 1) * ostr[i], istr, (char *) iptr + (off - 1) * istr[i], size);
#else

	md_swap2(D, dims, ostr, optr, istr, iptr, size);
	md_flip_inpl2(D, dims, flags, ostr, optr, size);
	md_flip_inpl2(D, dims, flags, istr, iptr, size);
#endif
}

/**
 * Swap input and output while flipping selected dimensions
 * at the same time.
 */
void md_swap_flip(unsigned int D, const long dims[], unsigned long flags, void* optr, void* iptr, size_t size)
{
	long strs[D];
	md_calc_strides(D, strs, dims, size);
	md_swap_flip2(D, dims, flags, strs, optr, strs, iptr, size);
}




static void md_flip_inpl2(unsigned int D, const long dims[], unsigned long flags, const long str[], void* ptr, size_t size)
{
	int i;
	for (i = D - 1; i >= 0; i--)
		if ((1 != dims[i]) && MD_IS_SET(flags, i))
			break;

	if (-1 == i)
		return;

	assert(1 < dims[i]);
	assert(str[i] != 0);

	long dims2[D];
	md_copy_dims(D, dims2, dims);
	dims2[i] = dims[i] / 2;

	long off = str[i] * (0 + (dims[i] + 1) / 2);
	md_swap_flip2(D, dims2, flags, str, ptr, str, (char *) ptr + off, size);
}

#if 0

/**
 * Copy array (with strides)
 *
 * optr[i] = iptr[i]
 */
void md_copy2(unsigned int D, const long dim[], const long ostr[], void* optr, const long istr[], const void* iptr, size_t size)
{

	long tostr[D];
	long tistr[D];
	long tdims[D];

	md_copy_strides(D, tostr, ostr);
	md_copy_strides(D, tistr, istr);
	md_copy_dims(D, tdims, dim);

	long (*nstr2[2])[D] = { &tostr, &tistr };
#if 0
	size_t sizes[2] = { size, size };
	int ND = optimize_dims(2, D, tdims, nstr2);
	int skip = min_blockdim(2, ND, tdims, nstr2, sizes); 
#else
    int skip = 0;
    int ND = D;
#endif
	const long* nstr[2] = { *nstr2[0] + skip, *nstr2[1] + skip };

	void* nptr[2] = { optr, (void*)iptr };
#ifdef  USE_CUDA
	struct data_s data = { md_calc_size(skip, tdims) * size, (cuda_ondevice(optr) || cuda_ondevice(iptr)) };

#if 1
	if (data.use_gpu && (ND - skip == 1)) { 

		long sizes[2] = { md_calc_size(skip, tdims) * size, tdims[skip] };
		struct strided_copy_s data = { { sizes[0], sizes[1] } , (*nstr2[0])[skip], (*nstr2[1])[skip] };

		skip++;
		md_nary(2, ND - skip, tdims + skip , nstr, nptr, (void*)&data, &nary_strided_copy);
		return;
	}
#endif
#else
	struct data_s data = { md_calc_size(skip, tdims) * size };
#endif

	md_nary(2, ND - skip, tdims + skip, nstr, nptr, (void*)&data, &nary_copy);
}

/**
 * Copy array (without strides)
 *
 * optr[i] = iptr[i]
 */
void md_copy(unsigned int D, const long dim[], void* optr, const void* iptr, size_t size)
{
	long str[D];
	md_calc_strides(D, str, dim, size);
	md_copy2(D, dim, str, optr, str, iptr, size);
}

#endif

/**
 * Fill array with value pointed by pointer (with strides)
 *
 * ptr[i] = iptr[0]
 */
void md_fill2(unsigned int D, const long dim[], const long str[], void* ptr, const void* iptr, size_t size)
{
#ifdef USE_CUDA
	if (cuda_ondevice(ptr) && (!cuda_ondevice(iptr))) {

		void* giptr = gpu_constant(iptr, size);
		md_fill2(D, dim, str, ptr, giptr, size);
		md_free(giptr);
		return;
	}
#endif
#if 0
    // TODO Test this
	long istr[D];
	md_singleton_strides(D, istr);
	md_copy2(D, dim, str, ptr, istr, iptr, size);
#else
    long sz = md_calc_size(D, dim);
    for( long i = 0 ; i < sz ; i++ ){
        memcpy(((char *) ptr + i*size), ((char *)iptr), size);
    }
#endif

}



/**
 * Fill array with value pointed by pointer (without strides)
 *
 * ptr[i] = iptr[0]
 */
void md_fill(unsigned int D, const long dim[], void* ptr, const void* iptr, size_t size)
{
	long str[D];
	md_calc_strides(D, str, dim, size);
	md_fill2(D, dim, str, ptr, iptr, size);
}



struct swap_s {

	unsigned int M;
	size_t size;
};


static void nary_swap(void* _data, void* ptr[])
{
	const struct swap_s* data = (const struct swap_s*) _data;
	size_t size = data->size;
	unsigned int M = data->M;
#ifdef MEX_COMPILE
    char* tmp = (size < 32) ? (char *) alloca(size) : (char *) mxMalloc(size);
#else
    char* tmp = (size < 32) ? (char *) alloca(size) : (char *) malloc(size);
#endif
    if( tmp == NULL ){
        fprintf(stderr, "Could not allocate in nary_swap\n");
        exit(0);
    }
#ifdef  USE_CUDA
	assert(!cuda_ondevice(ptr[0]));
	assert(!cuda_ondevice(ptr[1]));
#endif
	memcpy(tmp, ptr[0], size);

	for (unsigned int i = 0; i < M - 1; i++)
		memcpy(ptr[i], ptr[i + 1], size);

	memcpy(ptr[M - 1], tmp, size);

	if (size >= 32)
		free(tmp);
}

/**
 * Swap values between a number of arrays (with strides)
 */
void md_circular_swap2(unsigned M, unsigned int D, const long dims[], const long* strs[], void* ptr[], size_t size)
{
	unsigned int skip = md_calc_blockdim(D, dims, strs[0], size);

	for (unsigned int i = 1; i < M; i++)
		skip = MIN(skip, md_calc_blockdim(D, dims, strs[i], size));

	const long* nstr[M];
	for (unsigned int i = 0; i < M; i++)
		nstr[i] = strs[i] + skip;

	struct swap_s data = { M, md_calc_size(skip, dims) * size };
	md_nary(M, D - skip, dims + skip, nstr, ptr, (void*)&data, &nary_swap);
}


/**
 * Swap values between a number of arrays
 */
void md_circular_swap(unsigned M, unsigned int D, const long dims[], void* ptr[], size_t size)
{
	long strs[M][D];

	md_calc_strides(D, strs[0], dims, size);

	const long* strp[M];

	for (unsigned int i = 1; i < M; i++) {

		md_copy_strides(D, strs[i], strs[0]);
		strp[i] = strs[i];
	}

	md_circular_swap2(M, D, dims, strp, ptr, size);
}



#if 0

/**
 * Flip array (with strides)
 *
 * optr[ dims[D] - 1 - i ] = iptr[ i ]
 *
 */
void md_flip2(unsigned int D, const long dims[], unsigned long flags, const long ostr[], void* optr, const long istr[], const void* iptr, size_t size)
{
	if (optr == iptr) {

		assert(ostr == istr);
		md_flip_inpl2(D, dims, flags, ostr, optr, size);
		return;
	}

	long off = 0;
	long ostr2[D];

	for (unsigned int i = 0; i < D; i++) {

		ostr2[i] = ostr[i];

		if (MD_IS_SET(flags, i)) {

			ostr2[i] = -ostr[i];
			off += (dims[i] - 1) * ostr[i];
		} 
	}

	md_copy2(D, dims, ostr2, (char *) optr + off, istr, iptr, size);
}


/**
 * Flip array (without strides)
 *
 * optr[ dims[D] - 1 - i ] = iptr[ i ]
 *
 */
void md_flip(unsigned int D, const long dims[], unsigned long flags, void* optr, const void* iptr, size_t size)
{
	long str[D];
	md_calc_strides(D, str, dims, size);
	md_flip2(D, dims, flags, str, optr, str, iptr, size);
}

#endif

