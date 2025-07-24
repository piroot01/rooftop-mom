cmake_minimum_required(VERSION 3.10)
project(fastZeroMaskMult LANGUAGES C CXX)

# 1) Turn on compile-command export for clangd, etc.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 2) Require at least C11 and C++11 (MATLAB's API is C but we use std::vector, etc.)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 3) Locate your MATLAB installation via environment variable
if(NOT DEFINED ENV{MATLAB_ROOT})
    message(FATAL_ERROR
      "Please set the MATLAB_ROOT environment variable to your MATLAB install path,\n"
      "e.g. export MATLAB_ROOT=/usr/local/MATLAB/R2024b")
endif()
set(MATLAB_ROOT $ENV{MATLAB_ROOT})

# 4) Include headers and library paths
include_directories("${MATLAB_ROOT}/extern/include")
link_directories("${MATLAB_ROOT}/bin/glnxa64")

# 5) Treat .c file as C++ so we can use <vector>
set_source_files_properties(fastZeroMaskMult.c PROPERTIES LANGUAGE CXX)

# 6) Define the MEX shared library (no lib prefix, .mexa64 suffix)
add_library(fastZeroMaskMult SHARED fastZeroMaskMult.c)
set_target_properties(fastZeroMaskMult PROPERTIES
    PREFIX ""
    SUFFIX ".mexa64"
)

# 7) Link against the MEX and MATLAB libraries
# Order matters on some linkers: mex mx mat ut
target_link_libraries(fastZeroMaskMult
    PRIVATE mex mx mat ut
)

# 8) (Optional) Install rule to drop the .mexa64 into MATLAB's working dir
install(TARGETS fastZeroMaskMult
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
