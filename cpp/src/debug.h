#ifndef __DEBUG_H
#define __DEBUG_H 1

//#include "misc.h"
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdarg.h>     /* va_list, va_start, va_arg, va_end */

extern int debug_level;
extern bool debug_logging;

enum debug_levels { DP_ERROR, DP_WARN, DP_INFO, DP_DEBUG1, DP_DEBUG2, DP_DEBUG3, DP_DEBUG4, DP_ALL };

extern void debug_printf(int level, const char* fmt, ...);
extern void debug_vprintf(int level, const char* fmt, va_list ap);

void logerror(const char *fmt, ... );

void myAssert( int exp, const char *err);

int is_in_bounds(const unsigned int D, const long *d, const long *dims);

void assert_is_binary(const long N, const int *arr, const char *arrstr);

void assert_in_bounds(const unsigned int D, const long *d, 
                        const long *dims, const char *str);

void debug_print_dims(const int level, const unsigned int D, const long *d1);

void assert_same_dims(const unsigned int D, const char *d1str, const long *d1, const char *d2str, const long *d2);

#ifdef __cplusplus
}
#endif

/*
*/

#endif // __DEBUG_H