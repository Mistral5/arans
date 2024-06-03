cmake_minimum_required(VERSION 3.22)
project(arans C)

set(CMAKE_C_STANDARD 11)

add_compile_options(-mavx -march=native)

add_executable(arans
        main.c
        arans_8.h
        arans_2x6.h
        arans_3x5.h
        arans_4x4.h
        arans_2x2x4.h
        arans_2x3x3.h
        arans_2x2x2x2.h
        arans_3x5_clear.h
        arans_SIMD.h
)
