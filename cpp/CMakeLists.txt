cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
     
set(CMAKE_CXX_STANDARD 17)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

project (sampling)

add_subdirectory(src)
add_subdirectory(bcfast)
add_subdirectory(getDeltaJ)
add_subdirectory(testHeap)
add_subdirectory(testGenerator)

