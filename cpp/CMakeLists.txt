cmake_minimum_required(VERSION 3.4)


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

#set(SOURCES
#    src/dda_utils.cpp
#    src/sampleHeap.cpp
#    src/sample.cpp
#    src/misc.cpp
#    src/multind.cpp
#    src/debug.cpp
#)

add_library(${PROJECT_NAME} ${SOURCES} ${HEADERS})

#include_directories(

#target_public_includes(${PROJECT_NAME}
#    HEADERS
#        ${HEADERS}
#    DONT_STRIP
#    PREFIX bcfast
#)
#set(Boost_USE_STATIC_LIBS OFF) 
#set(Boost_USE_MULTITHREADED ON)  
#set(Boost_USE_STATIC_RUNTIME OFF) 
#find_package(Boost 1.59.0 COMPONENTS programoptions) 
find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

#find_package( Boost 1.45.0 COMPONENTS program_options REQUIRED )

#include_directories(include)
include_directories("${CMAKE_SOURCE_DIR}/src" ${Boost_INCLUDE_DIRS} )

#include_directories(${Boost_INCLUDE_DIRS}) 
#include_directories(src)



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

#target_link_libraries( bcfast ${PROJECT_NAME} ${SOURCES} )


#add_executable(getDeltaJ.cpp DEPS ${PROJECT_NAME} TEST)
    


