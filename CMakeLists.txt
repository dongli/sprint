cmake_minimum_required(VERSION 3.6)

project(sprint LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
endif ()
get_directory_property(parent_dir PARENT_DIRECTORY)

if (NOT parent_dir)
  add_subdirectory(lib/proj)
  add_subdirectory(lib/container)
  add_subdirectory(lib/kdtree)
  add_subdirectory(lib/string)
  add_subdirectory(lib/fiona)
  add_subdirectory(lib/flogger)
endif ()

set(sources
  src/bilinear_interpolator_mod.F90
)

add_library(sprint ${sources})
target_link_libraries(sprint fortran_proj fortran_container fortran_string kdtree fiona flogger fortran_proj)

if (NOT parent_dir)
  add_executable(bilinear_interpolator_test.exe src/bilinear_interpolator_test.F90)
  target_link_libraries(bilinear_interpolator_test.exe sprint)
endif ()
