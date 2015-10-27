# Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu>
# Licensed under the MIT license. See the license file LICENSE.

project(optRot)
cmake_minimum_required(VERSION 2.8.1)
set(POD_NAME optRot)
include(cmake/pods.cmake)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_COLOR_MAKEFILE ON)
IF(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) # Release; debug; ...
  set(CMAKE_BUILD_TYPE Debug) # Release; debug; ...
ENDIF()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif()

find_package (Eigen3 REQUIRED)
find_package (CUDA REQUIRED)

# add flags depending on build type
set(CMAKE_CXX_FLAGS_RELEASE "-g -Wall -Wno-deprecated -std=c++11 -O2 -finline-functions -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE "-g -Wall -Wno-deprecated -O2 -finline-functions")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wno-deprecated -std=c++11")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\")

include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  )

add_library(optRot SHARED
  ./src/tetrahedron.cpp
  ./src/node.cpp
  ./src/node_S3.cpp
  ./src/node_R3.cpp
  ./src/box.cpp
  ./src/combinations.cpp
  ./src/s3_tessellation.cpp
  ./src/numeric_helpers.cpp
  ./src/lower_bound_S3.cpp
  ./src/upper_bound_indep_S3.cpp
  ./src/upper_bound_convex_S3.cpp
  ./src/lower_bound_R3.cpp
  ./src/upper_bound_indep_R3.cpp
  ./src/upper_bound_convex_R3.cpp
#  ./src/branch_and_bound.cpp
  )

add_executable(testCombinations ./test/combinations.cpp)
target_link_libraries(testCombinations optRot)

add_executable(testS3Tessellation ./test/s3_tessellation.cpp)
target_link_libraries(testS3Tessellation optRot)

add_executable(testBoundsS3 ./test/bounds_S3.cpp)
target_link_libraries(testBoundsS3 optRot)

add_executable(testBoundsR3 ./test/bounds_R3.cpp)
target_link_libraries(testBoundsR3 optRot)

add_executable(testBranchAndBoundS3 ./test/branch_and_bound_S3.cpp)
target_link_libraries(testBranchAndBoundS3 optRot)

add_executable(testBranchAndBoundR3 ./test/branch_and_bound_R3.cpp)
target_link_libraries(testBranchAndBoundR3 optRot)

add_executable(testBranchAndBoundT3 ./test/branch_and_bound_T3.cpp)
target_link_libraries(testBranchAndBoundT3 optRot)

add_executable(testMinMaxInBox ./test/min_max_in_box.cpp)
target_link_libraries(testMinMaxInBox optRot)

pods_install_libraries(optRot)
pods_install_headers(
  ./include/optRot/s3_tessellation.h
  ./include/optRot/tetrahedron.h
  ./include/optRot/box.h
  ./include/optRot/combinations.h
  ./include/optRot/bound.h
  ./include/optRot/node.h
  ./include/optRot/node_impl.h
  ./include/optRot/node_S3.h
  ./include/optRot/node_R3.h
  ./include/optRot/branch_and_bound.h
  ./include/optRot/branch_and_bound_impl.h
  ./include/optRot/numeric_helpers.h
  ./include/optRot/vmf.h
  ./include/optRot/vmf_impl.h
  ./include/optRot/vmf_mm.h
  ./include/optRot/vmf_mm_impl.h
  ./include/optRot/normal.h
  ./include/optRot/normal_impl.h
  ./include/optRot/lower_bound_R3.h
  ./include/optRot/upper_bound_indep_R3.h
  ./include/optRot/upper_bound_convex_R3.h
  ./include/optRot/lower_bound_S3.h
  ./include/optRot/upper_bound_indep_S3.h
  ./include/optRot/upper_bound_convex_S3.h
  DESTINATION optRot)
