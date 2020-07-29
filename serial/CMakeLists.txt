cmake_minimum_required(VERSION 3.15)
project(serial C)

set(CMAKE_C_STANDARD 11)

add_executable(serial main.c graph.h graph.c coarsen.h initpart.h refine.h refine.c initpart.c coarsen.c multilevel.h multilevel.c utils.h utils.c)