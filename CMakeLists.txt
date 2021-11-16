cmake_minimum_required(VERSION 3.11)

set( CMAKE_CXX_COMPILER "g++" )
set( DEBUG "1")
project(dmc)

set(PROJECT_DIR .)

add_library(dmc dmcpp.cpp)

target_include_directories(dmc PUBLIC ${PROJECT_DIR})

#add_executable(debug dmcpp.cpp)
