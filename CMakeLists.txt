cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(coalescence_afterburner)

find_package(Git REQUIRED)

execute_process(COMMAND git clone  https://github.com/smash-transport/smash.git
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set(SOURCE_FILES
    src/coalescence.cc

    smash/src/pdgcode.cc
    smash/src/threevector.cc
    smash/src/fourvector.cc
    smash/src/stringfunctions.cc
)
add_executable(coalescence src/coalescence_main.cc ${SOURCE_FILES})

include_directories(include)
include_directories(smash/src/include)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -std=c++11 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Werror=switch")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
