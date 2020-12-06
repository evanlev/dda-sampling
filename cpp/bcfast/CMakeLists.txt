add_executable (bcfast bcfast.cpp)

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries (bcfast LINK_PUBLIC common ${Boost_LIBRARIES})

