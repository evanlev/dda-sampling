
#include <string.h> 

#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <time.h>
#include "cppmap.h"

#define STRSIZE 64


using namespace std;

#include "debug.h"

int debug_level = -1;
bool debug_logging = false;

void logerror(const char *fmt, ... )
{
    va_list ap;
    va_start(ap, fmt);

    char buff[512] = "";
    
    vsprintf(buff, fmt, ap);

    int openok = 1;
    cerr << buff << endl;
    va_end(ap);
}

static const char* level_strings[] = {
    "DP_ERROR",
    "DP_WARN",
    "DP_INFO",
    "DP_DEBUG1",
    "DP_DEBUG2",
    "DP_DEBUG3",
    "DP_DEBUG4"
};

static const char* get_level_str(int level)
{
        assert(level >= 0);
        return (level <= DP_DEBUG4) ? level_strings[level] : "ALL";
}

static void get_datetime_str(int len, char* datetime_str)
{
        time_t tv = time(NULL);
            struct tm* dt = gmtime(&tv);

                strftime(datetime_str, len, "%F %T", dt);
}

void debug_vprintf(int level, const char* fmt, va_list ap)
{
    if (-1 == debug_level) {

        char* str = getenv("DEBUG_LEVEL");
        debug_level = (NULL != str) ? atoi(str) : DP_INFO;
    }

    if (level <= debug_level) {

        FILE* ofp = (level < DP_INFO) ? stderr : stdout;

        if (true == debug_logging) {

            char dt_str[STRSIZE];
            get_datetime_str(STRSIZE, dt_str);

            fprintf(ofp, "[%s] [%s] - ", dt_str, get_level_str(level));

        } else
        if (debug_level < DP_INFO)
            fprintf(ofp, "%s: ", get_level_str(level));

        vfprintf(ofp, fmt, ap);
        
        fflush(ofp);
    }
}
 
void debug_printf(int level, const char *fmt, ... )
{
        va_list ap;
        va_start(ap, fmt);
        
        debug_vprintf(level, fmt, ap);
        va_end(ap);
}
    


void debug_print_dims(const int level, const unsigned int D, const long *d1){
    debug_printf(level, "[");
    for( unsigned int i = 0 ; i < D ; i++ ){
        debug_printf(level, "%ld ", d1[i]);
    }
    debug_printf(level, "]");
}


int is_in_bounds(const unsigned int D, const long *d, const long *dims){
    int i;
    for( i =0 ; i < D ;i++)
    {
        if( d[i] < 0 )
            return 0;
        if( d[i] >= dims[i] )
            return 0;
    }
    return 1;
}
void assert_same_dims(const unsigned int D, const char *d1str, const long *d1, const char *d2str, const long *d2){
    for( unsigned int i = 0 ; i < D ; i++ ){
        if(d1[i] != d2[i] ){
            debug_printf(DP_ERROR, "%s dims: ", d1str);
            debug_print_dims(DP_ERROR, D, d1);
            debug_printf(DP_ERROR, "\n");
            debug_printf(DP_ERROR, "%s dims: ", d2str);
            debug_print_dims(DP_ERROR, D, d2);
            debug_printf(DP_ERROR, "\n");
            logerror("");
        }
    }

}

void assert_is_binary(const long N, const int *arr, const char *arrstr){
    for( long i = 0 ; i < N ; i++ ){
        if( arr[i] != 0 && arr[i] != 1){
            debug_printf(DP_ERROR, "%s is not binary", arrstr);
            logerror("");
        }
    }
}

void assert_in_bounds(const unsigned int D, const long *d, const long *dims, const char *str){
    if(!is_in_bounds(D, d, dims)){
        debug_printf(DP_ERROR, "%s=", str);
        debug_print_dims(DP_ERROR, D, d);
        debug_printf(DP_ERROR, " is out of bounds ");
        debug_print_dims(DP_ERROR, D, dims);
        logerror("");
    }
}

void myAssert( int exp, const char *err){
    if( !exp ){
        std::cerr << err << std::endl;
        assert(0);
    }
}


