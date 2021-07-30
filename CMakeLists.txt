cmake_minimum_required(VERSION 3.17)
project(hexarena_backend CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)
include_directories(engine)

add_executable(hexarena_backend main.cpp)
