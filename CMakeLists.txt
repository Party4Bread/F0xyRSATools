cmake_minimum_required(VERSION 3.10)
project(F0xyRSATools C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "-lgmp")

add_executable(F0xyRSATools main.c pollard.c stdafx.h kraitchik.c fermat.c)
find_path(GMP_INCLUDE_DIRS NAMES gmp.h)
find_library(GMP_LIBRARY NAMES gmp libgmp)
target_link_libraries(F0xyRSATools gmp libgmp)