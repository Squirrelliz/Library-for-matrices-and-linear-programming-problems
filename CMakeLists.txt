cmake_minimum_required(VERSION 3.24)
project(iop)

set(CMAKE_CXX_STANDARD 23)

add_executable(iop main.cpp Fraction.cpp Fraction.h Matrix.cpp Matrix.h)
