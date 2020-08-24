cmake_minimum_required(VERSION 3.15)
project(parallel C)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})


add_executable(parallel main.c graph.h graph.c coarsen.h initpart.h refine.h refine.c initpart.c coarsen.c multilevel.h multilevel.c utils.h utils.c)
target_link_libraries(parallel ${MPI_C_LIBRARIES})
