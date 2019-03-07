project(getDeltaJ)
cmake_minimum_required(VERSION 3.1)

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
include_directories("${CMAKE_SOURCE_DIR}/src" ${Boost_INCLUDE_DIRS} )

add_executable(${PROJECT_NAME} getDeltaJ.cpp)
target_link_libraries(${PROJECT_NAME} common ${Boost_LIBRARIES})
