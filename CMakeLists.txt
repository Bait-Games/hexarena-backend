cmake_minimum_required(VERSION 3.17)
project(hexarena_backend CXX)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)
include_directories(engine)

add_executable(hexarena_backend main.cpp engine/user.cpp engine/component.h engine/utils.h engine/engine.cpp)
