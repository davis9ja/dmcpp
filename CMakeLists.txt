cmake_minimum_required(VERSION 3.11)

set( CMAKE_CXX_COMPILER "g++" )
project(dmc)

set(PROJECT_DIR .)

add_library(dmc dmcpp.cpp)

target_include_directories(dmc PUBLIC ${PROJECT_DIR})
