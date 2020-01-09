cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(coalescence_afterburner)

set(SOURCE_FILES
    src/coalescence.cc
    src/fourvector.cc
)
add_executable(coalescence src/coalescence_main.cc ${SOURCE_FILES})
include_directories(include)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -std=c++11 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Werror=switch")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
