cmake_minimum_required(VERSION 3.13)
project(CLBM)
# set compiler
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

set ( CMAKE_C_FLAGS "-O3 -march=native -mavx -ftree-vectorize -DLIKWID_PERFMON")
set ( CMAKE_CXX_FLAGS "-O3 -march=native -mavx -ftree-vectorize -DLIKWID_PERFMON")

# Find GSL package
find_package(GSL REQUIRED)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Find likwid package
find_package(likwid REQUIRED)

# Add your target executable
add_executable(CLBM src/main.cpp)

# Link with GSL and OpenMP
# Link with GSL and OpenMP
include_directories("$ENV{LIKWID_INCDIR}")
target_link_libraries(CLBM PRIVATE GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX "$ENV{LIKWID_LIBDIR}/liblikwid.so")

# Include directories
target_include_directories(CLBM PRIVATE ${GSL_INCLUDE_DIRS})

# Set compiler flags for OpenMP and optimization
if(OpenMP_CXX_FOUND)
    target_compile_options(CLBM PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(CLBM PRIVATE OpenMP::OpenMP_CXX)
endif()

# cmake_minimum_required(VERSION 3.10)
# project(CLBM)
# set compiler
# set(CMAKE_C_COMPILER gcc)
# set(CMAKE_CXX_COMPILER g++)
# 
# set ( CMAKE_C_FLAGS "-O3 -march=native")
# set ( CMAKE_CXX_FLAGS "-O3 -march=native")
# Find GSL package
# find_package(GSL REQUIRED)
# 
# Find OpenMP package
# find_package(OpenMP REQUIRED)
# 
# Add your target executable
# add_executable(CLBM src/main.cpp)
# add_executable(GRID grid/gridGen.cpp)
# Link with GSL and OpenMP
# target_link_libraries(CLBM PRIVATE GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX)
# 
# Include directories
# target_include_directories(CLBM PRIVATE ${GSL_INCLUDE_DIRS})
# 
# Set compiler flags for OpenMP and optimization
# if(OpenMP_CXX_FOUND)
#     target_compile_options(CLBM PRIVATE ${OpenMP_CXX_FLAGS}) 
#     target_link_libraries(CLBM PRIVATE OpenMP::OpenMP_CXX)
# endif()
# 
