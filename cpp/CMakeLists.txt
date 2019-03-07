cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
     
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

project(sampling)

set(HEADERS
    src/config.h
    src/dda_utils.h
    src/sample.h
    src/sampleHeap.h
    src/debug.h
    src/mdarray.h
    src/misc.h
    src/misc.hpp
    src/multind.h
    src/cppmap.h
)

file( GLOB SOURCES src/*.cpp )

add_library(${PROJECT_NAME} ${SOURCES} ${HEADERS})

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/src" ${Boost_INCLUDE_DIRS} )

#file( GLOB APP_SOURCES src/*.cxx )
#file( GLOB APP_SOURCES *.cxx )
#foreach( testsourcefile ${APP_SOURCES} )
    #string( REPLACE ".cxx" "" testname ${testsourcefile} )
    #add_executable( ${testname} ${testsourcefile} )
    #target_link_libraries(${testname} ${PROJECT_NAME} ${Boost_LIBRARIES} )
    #endforeach( testsourcefile ${APP_SOURCES})

add_executable( bcfast bcfast.cpp )
target_link_libraries( bcfast ${PROJECT_NAME} ${SOURCES} ${Boost_LIBRARIES} )

add_executable( getDeltaJ getDeltaJ.cpp )
target_link_libraries( getDeltaJ ${PROJECT_NAME} ${SOURCES} ${Boost_LIBRARIES} )

add_executable( testHeap testHeap.cpp )
target_link_libraries( testHeap ${PROJECT_NAME} ${SOURCES} ${Boost_LIBRARIES} )

add_executable( testGenerator testGenerator.cpp )
target_link_libraries( testGenerator ${PROJECT_NAME} ${SOURCES} ${Boost_LIBRARIES} )
