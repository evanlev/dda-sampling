file( GLOB SOURCES *.cpp )

set(HEADERS
    config.h
    dda_utils.h
    sample.h
    sampleHeap.h
    debug.h
    mdarray.h
    misc.h
    misc.hpp
    multind.h
    cppmap.h
)

add_library(common STATIC ${SOURCES})
target_include_directories(common PUBLIC "${CMAKE_CURRENT_LIST_DIR}")
