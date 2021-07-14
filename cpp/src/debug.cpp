
#include <string>
#include <array>

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

#include "debug.h"

int debug_level = -1;
bool debug_logging = false;

void logerror(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);

    char buff[512] = "";

    vsprintf(buff, fmt, ap);

    int openok = 1;
    std::cerr << buff << std::endl;
    va_end(ap);
}

static constexpr std::array<const char*, 7> s_level_strings = {
    "DP_ERROR",
    "DP_WARN",
    "DP_INFO",
    "DP_DEBUG1",
    "DP_DEBUG2",
    "DP_DEBUG3",
    "DP_DEBUG4"};

static constexpr const char* get_level_str(int level)
{
    assert(level >= 0);
    return (level <= DP_DEBUG4) ? s_level_strings[level] : "ALL";
}

static void get_datetime_str(int len, char *datetime_str)
{
    time_t tv = time(NULL);
    struct tm *dt = gmtime(&tv);

    strftime(datetime_str, len, "%F %T", dt);
}

void debug_vprintf(int level, const char *fmt, va_list ap)
{
    if (-1 == debug_level)
    {
        char *str = getenv("DEBUG_LEVEL");
        debug_level = (NULL != str) ? atoi(str) : DP_INFO;
    }

    if (level <= debug_level)
    {
        FILE *ofp = (level < DP_INFO) ? stderr : stdout;

        if (true == debug_logging)
        {

            char dt_str[STRSIZE];
            get_datetime_str(STRSIZE, dt_str);

            fprintf(ofp, "[%s] [%s] - ", dt_str, get_level_str(level));
        }
        else if (debug_level < DP_INFO)
            fprintf(ofp, "%s: ", get_level_str(level));

        vfprintf(ofp, fmt, ap);

        fflush(ofp);
    }
}

void debug_printf(int level, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);

    debug_vprintf(level, fmt, ap);
    va_end(ap);
}

void debug_print_dims(const int level, const unsigned int D, const long *d1)
{
    debug_printf(level, "[");
    for (unsigned int i = 0; i < D; i++)
    {
        debug_printf(level, "%ld ", d1[i]);
    }
    debug_printf(level, "]");
}

void assert_same_dims(const unsigned int D, const char *d1str, const long *d1, const char *d2str, const long *d2)
{
    for (unsigned int i = 0; i < D; i++)
    {
        if (d1[i] != d2[i])
        {
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

void assert_is_binary(const long N, const int *arr, const char *arrstr)
{
    for (long i = 0; i < N; i++)
    {
        if (arr[i] != 0 && arr[i] != 1)
        {
            debug_printf(DP_ERROR, "%s is not binary", arrstr);
            logerror("");
        }
    }
}

void assert_with_log(int exp, const char *err)
{
    if (!exp)
    {
        std::cerr << err << std::endl;
        assert(0);
    }
}
