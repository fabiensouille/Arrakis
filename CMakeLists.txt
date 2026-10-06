cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(
	ARRAKIS 
	VERSION 2.6
	LANGUAGES Fortran C
)

# ENABLE OPENMP
option(ENABLE_OPENMP "Enable OpenMP support" ON)

# Specify the output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Recipe
set(CMAKE_CONFIGURATION_TYPES Release Optim Debug Profiling)
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_OPTIM "-O3 -DNDEBUG -flto -march=native -funroll-loops -ftree-vectorize")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -fbacktrace -fbounds-check -finit-integer=-1 -finit-real=nan -ffpe-trap=invalid,zero,overflow")
set(CMAKE_Fortran_FLAGS_PROFILING "-g -pg")

# Set Release as default if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Set the flags based on the build type
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")

# Sources code
add_subdirectory(sources)

# Parallel using OpenMP
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(ARRAKIS PUBLIC OpenMP::OpenMP_Fortran)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  else()
    message(FATAL_ERROR "OpenMP not found")
  endif()
else()
  string(REPLACE "-fopenmp" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# External dependency
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
if(WIN32)
  find_package(LAPACK COMPONENTS OpenBLAS)
else()
  find_package(LAPACK)
endif()
if (LAPACK_FOUND)
  target_include_directories(ARRAKIS PUBLIC ${LAPACK_INCLUDE_DIR})
  target_link_libraries (ARRAKIS PUBLIC LAPACK::LAPACK)
  target_compile_definitions(ARRAKIS PUBLIC LAPACK::LAPACK)
else ()
  message(FATAL_ERROR "LAPACK not found")
endif ()
