cmake_minimum_required(VERSION 3.22)
project(arans C)

set(CMAKE_C_STANDARD 11)

add_compile_options(-mavx -march=native)

add_executable(arans
        main.c
        arans.h
        arans_SIMD.h
        arans_16x16.h
        arans_8x32.h
        arans_4x64.h
)
