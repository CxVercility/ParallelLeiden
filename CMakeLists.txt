cmake_minimum_required(VERSION 3.16)
project(ParallelLeiden)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(networkit)
include_directories(networkit/include/)
include_directories(networkit/extlibs/tlx)
set(CMAKE_CXX_FLAGS "-fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -fopenmp -Wall -pedantic -Wextra -Wformat -Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wno-unused-parameter -Wpedantic -Wnull-dereference -fno-omit-frame-pointer -fdiagnostics-color -ftemplate-depth=64 -fconstexpr-depth=64 -ftemplate-backtrace-limit=8 -Wreorder -Woverloaded-virtual")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-fopenmp")
add_executable(Testrun Test.cpp)
add_library(PL SHARED ParallelLeiden.cpp PLExperimental.cpp PLSequential.cpp )
target_link_libraries(Testrun networkit PL)


